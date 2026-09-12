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

# The categorical latent variables in SuSiE are component assignments, not
# marginal predictor-inclusion indicators. All alpha rows belong in this
# objective, including rows with V=0 and rows without a reported credible set.
em_coding_availability <- function() {
  out <- sapply(0:2, function(j) bitwAnd(1:7, bitwShiftL(1L, j)) != 0L)
  colnames(out) <- c("additive", "recessive", "dominant")
  out
}

em_coding_q <- function(prior, counts) {
  available <- em_coding_availability()
  value <- 0
  for (i in which(rowSums(counts) > 0)) {
    probability <- prior[available[i, ]] / sum(prior[available[i, ]])
    mass <- counts[i, available[i, ]]
    if (any(!is.finite(probability)) || any(probability[mass > 0] <= 0)) return(-Inf)
    value <- value + sum(mass[mass > 0] * log(probability[mass > 0]))
  }
  value
}

em_coding_mstep <- function(counts, previous = rep(1/3, 3)) {
  available <- em_coding_availability()
  stopifnot(is.matrix(counts), identical(dim(counts), c(7L, 3L)),
            all(is.finite(counts)), all(counts >= 0), all(counts[!available] == 0),
            length(previous) == 3L, all(is.finite(previous)), all(previous >= 0),
            abs(sum(previous) - 1) < 1e-8)
  total <- sum(counts)
  if (total == 0) return(list(prior = previous, gain = 0, status = "carried_forward_no_components"))
  before <- em_coding_q(previous, counts)
  if (!is.finite(before)) stop("Component posterior has support outside the starting coding prior.")
  if (all(rowSums(counts)[1:6] == 0)) {
    # With all three classes present in every fit, the exact M-step is the
    # normalized sum of alpha, independent of unequal numbers of predictors.
    prior <- colSums(counts) / total
  } else {
    # Missing coding classes introduce -L_g log(sum_{c in A_g} pi_c) terms.
    # Optimize the concave objective in log weights, preserving old total
    # mass between disconnected sets of classes (which the data cannot identify).
    connected <- diag(TRUE, 3)
    for (i in which(rowSums(counts) > 0)) {
      j <- which(available[i, ])
      connected[j, j] <- TRUE
    }
    for (j in 1:3) connected <- connected | outer(connected[, j], connected[j, ], "&")
    prior <- previous
    remaining <- 1:3
    while (length(remaining)) {
      group <- which(connected[remaining[1], ])
      remaining <- setdiff(remaining, group)
      if (length(group) == 1L || sum(previous[group]) == 0) next
      rows <- which(rowSums(counts[, group, drop = FALSE]) > 0)
      a <- available[rows, group, drop = FALSE]
      z <- counts[rows, group, drop = FALSE]
      totals <- rowSums(z)
      unpack <- function(theta) c(theta, 0)
      objective <- function(theta, gradient = FALSE) {
        eta <- unpack(theta)
        q <- 0
        score <- colSums(z)
        for (r in seq_along(rows)) {
          eta_r <- eta[a[r, ]]
          log_denominator <- max(eta_r) + log(sum(exp(eta_r - max(eta_r))))
          q <- q + sum(z[r, a[r, ]] * (eta_r - log_denominator))
          score[a[r, ]] <- score[a[r, ]] - totals[r] * exp(eta_r - log_denominator)
        }
        if (gradient) -head(score, -1) / sum(totals) else -q / sum(totals)
      }
      start <- log(pmax(previous[group], 1e-12))
      opt <- optim(head(start - tail(start, 1), -1), objective,
                   gr = function(theta) objective(theta, TRUE), method = "BFGS",
                   control = list(reltol = 1e-13, maxit = 2000))
      if (opt$convergence != 0 || max(abs(objective(opt$par, TRUE))) > 1e-6) {
        stop("Coding-prior M-step did not converge for incomplete coding availability.")
      }
      eta <- unpack(opt$par)
      probability <- exp(eta - max(eta))
      prior[group] <- sum(previous[group]) * probability / sum(probability)
    }
  }
  gain <- em_coding_q(prior, counts) - before
  if (!is.finite(gain) || gain < -1e-8 * max(1, abs(before))) {
    stop("Coding-prior update decreased the expected complete log likelihood.")
  }
  list(prior = prior, gain = gain, status = "estimated")
}

em_bind_history <- function(history, priors) {
  if (is.null(history)) return(priors)
  # Preserve the meaning of old rows when extending a legacy PIP-share run.
  if (!"update_method" %in% names(history)) history$update_method <- "legacy_pip_share"
  columns <- union(names(history), names(priors))
  for (column in setdiff(columns, names(history))) history[[column]] <- NA
  for (column in setdiff(columns, names(priors))) priors[[column]] <- NA
  rbind(history[columns], priors[columns])
}

em_estimate_coding_priors <- function(files, fit_name, previous_priors = NULL) {
  classes <- c("additive", "recessive", "dominant")
  sums <- counts <- list()
  n_fits <- n_components <- n_zero_variance <- integer(0)
  elbo_sums <- n_elbo <- numeric(0)
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
      alpha <- fit[["alpha"]]
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
          any(pip < -1e-10 | pip > 1 + 1e-10) || length(coding) != length(pip) ||
          anyNA(coding) || any(!coding %in% classes)) {
        stop("Missing/invalid ", fit_name, " PIPs or coding labels in ", file, " / ", tissue)
      }
      if (!is.matrix(alpha) || !is.numeric(alpha) || nrow(alpha) < 1L ||
          ncol(alpha) != length(pip) || any(!is.finite(alpha)) ||
          any(alpha < -1e-10 | alpha > 1 + 1e-10) ||
          any(abs(rowSums(alpha) - 1) > 1e-8) ||
          (!is.null(fit$null_index) && fit$null_index != 0)) {
        stop("Missing/invalid component alpha in ", file, " / ", tissue,
             ". Full SuSiE fits without an explicit null column are required; PIPs alone cannot supply this EM update.")
      }
      expected_names <- if (!is.null(map)) as.character(map$predictor_name) else names(pip)
      if (!is.null(colnames(alpha)) && !is.null(expected_names) &&
          !identical(colnames(alpha), expected_names)) {
        stop("Alpha/map name mismatch in ", file, " / ", tissue)
      }
      if (!identical(fit$converged, TRUE)) {
        stop("Source fit is not confirmed converged in ", file, " / ", tissue,
             "; finish the E-step before updating coding priors.")
      }
      if (!is.numeric(fit$V) || !length(fit$V) %in% c(1L, nrow(alpha)) ||
          any(!is.finite(fit$V)) || any(fit$V < 0)) {
        stop("Missing/invalid component prior variances in ", file, " / ", tissue)
      }
      alpha[] <- pmin(1, pmax(0, alpha))
      alpha <- alpha / rowSums(alpha)
      pip <- pmin(1, pmax(0, pip))
      if (is.null(sums[[tissue]])) {
        sums[[tissue]] <- setNames(numeric(3), classes)
        counts[[tissue]] <- matrix(0, 7, 3, dimnames = list(NULL, classes))
        n_fits[tissue] <- n_components[tissue] <- n_zero_variance[tissue] <- 0L
        elbo_sums[tissue] <- n_elbo[tissue] <- 0
      }
      # PIP sums remain descriptive diagnostics only. The M-step uses alpha.
      sums[[tissue]] <- sums[[tissue]] + vapply(classes, function(c) sum(pip[coding == c]), numeric(1))
      availability <- sum(2^(0:2) * (classes %in% coding))
      counts[[tissue]][availability, ] <- counts[[tissue]][availability, ] +
        vapply(classes, function(c) sum(alpha[, coding == c, drop = FALSE]), numeric(1))
      n_fits[tissue] <- n_fits[tissue] + 1L
      n_components[tissue] <- n_components[tissue] + nrow(alpha)
      n_zero_variance[tissue] <- n_zero_variance[tissue] + sum(rep_len(fit$V, nrow(alpha)) == 0)
      final_elbo <- tail(fit$elbo, 1)
      if (is.numeric(final_elbo) && length(final_elbo) == 1L && is.finite(final_elbo)) {
        elbo_sums[tissue] <- elbo_sums[tissue] + final_elbo
        n_elbo[tissue] <- n_elbo[tissue] + 1L
      }
    }
  }
  if (!length(sums)) stop("No usable tissue component posteriors found in source results.")
  tissues <- sort(unique(c(names(sums), previous_priors$tissue)))
  priors <- do.call(rbind, lapply(tissues, function(tissue) {
    mass <- sums[[tissue]]
    if (is.null(mass)) mass <- setNames(numeric(3), classes)
    total <- sum(mass)
    z <- counts[[tissue]]
    if (is.null(z)) z <- matrix(0, 7, 3)
    previous <- if (!is.null(previous_priors) && tissue %in% previous_priors$tissue)
      em_prior_for_tissue(previous_priors, tissue) else rep(1/3, 3)
    update <- em_coding_mstep(z, previous)
    p <- update$prior
    alpha_mass <- colSums(z)
    data.frame(tissue = tissue, pi_add = unname(p[1]), pi_rec = unname(p[2]),
               pi_dom = unname(p[3]), pip_add = unname(mass[1]),
               pip_rec = unname(mass[2]), pip_dom = unname(mass[3]), pip_total = total,
               n_fits = if (tissue %in% names(n_fits)) n_fits[[tissue]] else 0L,
               n_nonconverged = 0L,
               prior_status = update$status,
               update_method = "susie_component_alpha_v1",
               alpha_add = alpha_mass[1], alpha_rec = alpha_mass[2], alpha_dom = alpha_mass[3],
               alpha_total = sum(alpha_mass),
               n_components = if (tissue %in% names(n_components)) n_components[[tissue]] else 0L,
               n_zero_variance_components = if (tissue %in% names(n_zero_variance)) n_zero_variance[[tissue]] else 0L,
               mstep_q_gain = update$gain,
               source_elbo_sum = if (tissue %in% names(n_elbo) && n_elbo[tissue] == n_fits[tissue])
                 unname(elbo_sums[tissue]) else NA_real_,
               n_source_elbo = if (tissue %in% names(n_elbo)) unname(n_elbo[tissue]) else 0L,
               max_abs_prior_change = max(abs(p - previous)), stringsAsFactors = FALSE,
               row.names = NULL)
  }))
  em_validate_priors(priors)
  list(priors = priors, audit = audit, component_counts = counts)
}

# Do not copy the old fit's pi into s_init: some SuSiE versions let it override
# the newly supplied prior_weights. Also keep V outside s_init so that older
# versions do not prune/reinitialize the zero-variance components on warm start.
em_stop_fit <- function(message) {
  stop(structure(list(message = message, call = NULL), class = c("em_fit_error", "error", "condition")))
}

em_susie_initialization <- function(fit, predictor_names, L, variance_y) {
  if (is.null(fit)) return(NULL)
  fields <- c("alpha", "mu", "mu2")
  expected <- c(min(as.integer(L), length(predictor_names)), length(predictor_names))
  if (!all(vapply(fields, function(field) {
    is.matrix(fit[[field]]) && is.numeric(fit[[field]]) &&
      identical(dim(fit[[field]]), as.integer(expected)) && all(is.finite(fit[[field]])) &&
      identical(colnames(fit[[field]]), predictor_names)
  }, logical(1))) || any(fit$alpha < 0 | fit$alpha > 1) ||
      any(abs(rowSums(fit$alpha) - 1) > 1e-8) ||
      !is.numeric(fit$V) || !length(fit$V) %in% c(1L, expected[1]) ||
      any(!is.finite(fit$V)) || any(fit$V < 0) ||
      !is.numeric(fit$sigma2) || length(fit$sigma2) != 1L ||
      !is.finite(fit$sigma2) || fit$sigma2 <= 0 ||
      !is.finite(variance_y) || variance_y <= 0) {
    em_stop_fit("Previous fit cannot initialize this model: predictor order, L, or posterior parameters differ.")
  }
  s_init <- fit[fields]
  class(s_init) <- "susie"
  list(s_init = s_init, scaled_prior_variance = rep_len(fit$V, expected[1]) / variance_y,
       residual_variance = fit$sigma2)
}
