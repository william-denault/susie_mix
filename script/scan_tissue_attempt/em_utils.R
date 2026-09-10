# Base-R helpers shared by EM preparation, workers, and tests.
em_validate_priors <- function(priors) {
  cols <- c("pi_add", "pi_rec", "pi_dom")
  if (!is.data.frame(priors) || !all(c("tissue", cols) %in% names(priors)) ||
      !nrow(priors) || anyNA(priors$tissue) || any(!nzchar(priors$tissue)) ||
      anyDuplicated(priors$tissue) || !all(vapply(priors[cols], is.numeric, logical(1)))) {
    stop("Priors must contain unique tissue names and numeric pi_add, pi_rec, pi_dom.")
  }
  p <- as.matrix(priors[cols])
  if (any(!is.finite(p)) || any(p < 0) || any(abs(rowSums(p) - 1) > 1e-8)) {
    stop("Each tissue's priors must be finite, non-negative, and sum to one.")
  }
  invisible(priors)
}

em_prior_for_tissue <- function(priors, tissue) {
  i <- match(tissue, priors$tissue)
  if (is.na(i)) stop("No estimated coding prior for tissue: ", tissue)
  setNames(as.numeric(priors[i, c("pi_add", "pi_rec", "pi_dom")]),
           c("additive", "recessive", "dominant"))
}

em_predictor_weights <- function(coding, prior, predictor_names = NULL) {
  classes <- c("additive", "recessive", "dominant")
  if (!length(coding) || anyNA(coding) || any(!coding %in% classes) ||
      !is.numeric(prior) || !setequal(names(prior), classes) ||
      anyDuplicated(names(prior)) || any(!is.finite(prior)) || any(prior < 0)) {
    stop("Invalid coding labels or coding prior.")
  }
  # Preserve class mass despite different predictor counts after QC.
  # If a whole class is absent, redistribute across the remaining classes.
  counts <- table(coding)
  weights <- unname(prior[coding] / as.numeric(counts[coding]))
  if (sum(weights) <= 0) stop("Retained predictors have zero total prior mass.")
  weights <- weights / sum(weights)
  if (!is.null(predictor_names)) {
    if (length(predictor_names) != length(weights)) stop("Predictor name mismatch.")
    names(weights) <- predictor_names
  }
  weights
}

em_atomic_write <- function(object, path, csv = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".em_write_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  if (csv) write.csv(object, tmp, row.names = FALSE, na = "") else saveRDS(object, tmp)
  # Same-directory rename is atomic on the Linux cluster. Windows cannot
  # replace an existing target with file.rename, so use a copy there.
  if (.Platform$OS.type == "windows" && file.exists(path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) stop("Cannot write: ", path)
  } else if (!file.rename(tmp, path)) stop("Cannot rename output to: ", path)
  invisible(path)
}

em_read_manifest <- function(index_dir) {
  files <- sort(list.files(index_dir, "^chunk_[0-9]+_genes\\.txt$", full.names = TRUE))
  if (!length(files)) stop("No chunk gene lists found in ", index_dir)
  manifest <- do.call(rbind, lapply(seq_along(files), function(i) {
    genes <- trimws(readLines(files[i], warn = FALSE))
    if (!length(genes) || any(!nzchar(genes))) stop("Empty gene/chunk in ", files[i])
    data.frame(chunk = i, gene = genes, stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(manifest$gene) || any(grepl("[/\\\\]", manifest$gene)) ||
      any(manifest$gene %in% c(".", ".."))) {
    stop("Chunk gene names must be unique and valid file names.")
  }
  manifest
}

em_check_result_files <- function(results_dir, manifest) {
  files <- sort(list.files(results_dir, "\\.rds$", full.names = TRUE))
  genes <- sub("\\.rds$", "", basename(files))
  missing <- setdiff(manifest$gene, genes)
  extra <- setdiff(genes, manifest$gene)
  if (length(missing) || length(extra)) {
    stop("Result files do not match the gene chunks in ", results_dir,
         ". Missing: ", length(missing), "; extra: ", length(extra),
         ". Complete the source scan or reconcile the chunk lists first.")
  }
  files
}

em_pending_chunks <- function(iteration_dir, manifest) {
  chunks <- sort(unique(manifest$chunk))
  chunks[!file.exists(file.path(iteration_dir, "completed",
                               sprintf("chunk_%03d.done", chunks)))]
}

em_sum_pips <- function(files, fit_name, previous_priors = NULL) {
  classes <- c("additive", "recessive", "dominant")
  sums <- list()
  n_fits <- n_nonconverged <- integer(0)
  audit <- data.frame(gene = character(), tissue = character(),
                      issue = character(), message = character())
  add_issue <- function(gene, tissue, issue, message) {
    audit[nrow(audit) + 1L, ] <<- list(gene, tissue, issue, paste(message, collapse = "; "))
  }
  for (file in files) {
    gene <- sub("\\.rds$", "", basename(file))
    # Unreadable files are not silently treated as failed fits.
    out <- readRDS(file)
    if (!is.list(out)) stop("Invalid gene result: ", file)
    if (!is.null(out[["error"]])) {
      add_issue(gene, "", "gene_error", out[["error"]])
      next
    }
    errors <- attr(out, "tissue_errors", exact = TRUE)
    for (tissue in names(errors)) add_issue(gene, tissue, "tissue_error", errors[[tissue]])
    if (!length(out)) {
      add_issue(gene, "", "no_successful_tissues", "No tissue fits in result")
      next
    }
    if (is.null(names(out)) || any(!nzchar(names(out))) || anyDuplicated(names(out))) {
      stop("Invalid tissue names in ", file)
    }
    for (tissue in names(out)) {
      x <- out[[tissue]]
      fit <- x[[fit_name]]
      pip <- fit[["pip"]]
      map <- x[["mix_predictor_map"]]
      coding <- x[["mix_coding"]]
      if (!is.null(map)) {
        if (!all(c("predictor_index", "predictor_name", "coding") %in% names(map)) ||
            nrow(map) != length(pip) ||
            !identical(as.integer(map$predictor_index), seq_along(pip))) {
          stop("Invalid predictor map in ", file, " / ", tissue)
        }
        if (!is.null(coding) && !identical(as.character(coding), as.character(map$coding))) {
          stop("Coding/map mismatch in ", file, " / ", tissue)
        }
        coding <- as.character(map$coding)
        if (!is.null(names(pip)) && !identical(names(pip), as.character(map$predictor_name))) {
          stop("PIP/map name mismatch in ", file, " / ", tissue)
        }
      }
      if (!is.numeric(pip) || !length(pip) || any(!is.finite(pip)) ||
          any(pip < 0 | pip > 1) || length(coding) != length(pip) ||
          anyNA(coding) || any(!coding %in% classes)) {
        stop("Missing/invalid ", fit_name, " PIPs or coding labels in ", file, " / ", tissue)
      }
      if (is.null(sums[[tissue]])) {
        sums[[tissue]] <- setNames(numeric(3), classes)
        n_fits[tissue] <- n_nonconverged[tissue] <- 0L
      }
      # Include every predictor PIP, including those outside credible sets.
      sums[[tissue]] <- sums[[tissue]] + vapply(classes, function(c) sum(pip[coding == c]), numeric(1))
      n_fits[tissue] <- n_fits[tissue] + 1L
      if (identical(fit[["converged"]], FALSE)) {
        n_nonconverged[tissue] <- n_nonconverged[tissue] + 1L
        add_issue(gene, tissue, "nonconverged", "PIPs included; SuSiE reported converged=FALSE")
      }
    }
  }
  if (!length(sums)) stop("No usable tissue PIPs found in source results.")
  tissues <- sort(unique(c(names(sums), previous_priors$tissue)))
  priors <- do.call(rbind, lapply(tissues, function(tissue) {
    mass <- sums[[tissue]]
    if (is.null(mass)) mass <- setNames(numeric(3), classes)
    total <- sum(mass)
    status <- "estimated"
    if (total > 0) {
      p <- mass / total
    } else {
      if (is.null(previous_priors) || !tissue %in% previous_priors$tissue) {
        stop("Zero total PIP for ", tissue, "; initial proportions are undefined.")
      }
      p <- em_prior_for_tissue(previous_priors, tissue)
      status <- "carried_forward_zero_pip"
    }
    data.frame(tissue = tissue, pi_add = unname(p[1]), pi_rec = unname(p[2]),
               pi_dom = unname(p[3]), pip_add = unname(mass[1]),
               pip_rec = unname(mass[2]), pip_dom = unname(mass[3]), pip_total = total,
               n_fits = if (tissue %in% names(n_fits)) n_fits[[tissue]] else 0L,
               n_nonconverged = if (tissue %in% names(n_nonconverged)) n_nonconverged[[tissue]] else 0L,
               prior_status = status, stringsAsFactors = FALSE)
  }))
  em_validate_priors(priors)
  list(priors = priors, audit = audit)
}
