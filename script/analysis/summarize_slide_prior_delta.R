# Source this file, then call summarize_slide_prior_delta(path, workers = 4).
# Each worker reads one RDS at a time and returns only the small CS summary.
# Uses base R + parallel; data.table speeds up binding when installed.

.slide_delta_bind <- function(parts) {
  parts <- Filter(function(x) !is.null(x) && length(x) > 0L, parts)
  if (!length(parts)) return(data.frame())
  if (requireNamespace("data.table", quietly = TRUE))
    return(as.data.frame(data.table::rbindlist(parts, use.names = TRUE)))
  if (is.data.frame(parts[[1L]])) {
    ans <- do.call(rbind, parts)
    rownames(ans) <- NULL
    return(ans)
  }
  # Plain named lists of scalars: assemble columns, not one data.frame per CS.
  columns <- lapply(names(parts[[1L]]), function(nm)
    unlist(lapply(parts, `[[`, nm), use.names = FALSE))
  names(columns) <- names(parts[[1L]])
  as.data.frame(columns, stringsAsFactors = FALSE)
}

.slide_delta_summary <- function(weights, grid, prefix) {
  total <- sum(weights)
  if (!is.finite(total) || total <= 0 || any(!is.finite(weights)) ||
      any(weights < 0)) stop("Invalid posterior weights.")
  prob <- weights / total
  ans <- c(delta_mean = sum(prob * grid),
           p_recessive_like = sum(prob[grid < 0]),
           p_recessive = sum(prob[grid == -1]),
           p_additive = sum(prob[grid == 0]),
           p_dominant_like = sum(prob[grid > 0]))
  as.list(setNames(ans, paste0(prefix, names(ans))))
}

.slide_delta_file <- function(file, min_entries) {
  # A failed file is recorded explicitly; no partial rows from that file escape.
  tryCatch({
    out <- readRDS(file)
    if (length(out) < min_entries)
      return(list(rows = NULL, error = NULL, status = "filtered"))
    n_cs <- vapply(out, function(x) length(x$fit_slide_prior$sets$cs), integer(1))
    rows <- vector("list", sum(n_cs))
    i <- 0L
    for (o in which(n_cs > 0L)) {
      fit <- out[[o]]$fit_slide_prior
      tissue <- names(out)[o]
      if (!length(tissue) || is.na(tissue) || !nzchar(tissue)) tissue <- as.character(o)
      grid <- fit$delta_grid
      dims <- c(dim(fit$alpha), length(grid))
      if (!is.matrix(fit$alpha) || !length(grid) || any(!is.finite(grid)) ||
          !identical(dim(fit$alpha_delta), as.integer(dims)))
        stop("Missing or incompatible alpha, alpha_delta or delta_grid in ", tissue)
      p <- ncol(fit$alpha)
      indices <- fit$sets$cs_index
      if (length(indices) != n_cs[o] || anyNA(indices) ||
          any(indices < 1 | indices > nrow(fit$alpha) | indices != floor(indices)))
        stop("Invalid cs_index in ", tissue)
      if (length(fit$delta_forced) != p || !is.matrix(fit$genotype_counts) ||
          nrow(fit$genotype_counts) != p ||
          !all(c("0", "1", "2") %in% colnames(fit$genotype_counts)) ||
          length(fit$min_obs) != 1L)
        stop("Missing or incompatible forced/count metadata in ", tissue)

      for (s in seq_len(n_cs[o])) {
        cs <- fit$sets$cs[[s]]
        ell <- indices[s]
        if (!length(cs) || anyNA(cs) || any(cs < 1 | cs > p | cs != floor(cs)) ||
            anyDuplicated(cs)) stop("Invalid CS members in ", tissue)
        alpha_cs <- fit$alpha[ell, cs]
        coverage <- sum(alpha_cs)
        if (any(!is.finite(alpha_cs)) || any(alpha_cs < 0) || coverage <= 0)
          stop("Invalid alpha in ", tissue)
        lead <- cs[which.max(alpha_cs)]
        cs_name <- names(fit$sets$cs)[s]
        if (!length(cs_name) || is.na(cs_name) || !nzchar(cs_name))
          cs_name <- paste0("L", ell)

        # drop = FALSE protects single-SNP CSs and single-point delta grids.
        # R stores this slice in SNP-by-grid order after removing the size-1 axis.
        joint <- matrix(fit$alpha_delta[ell, cs, , drop = FALSE],
                        nrow = length(cs), ncol = length(grid))
        i <- i + 1L
        rows[[i]] <- c(list(
          file = basename(file), tissue = tissue, cs = cs_name,
          component = ell, n_snps = length(cs), cs_coverage = coverage,
          lead_index = unname(lead), lead_alpha = unname(fit$alpha[ell, lead]),
          lead_delta_forced = unname(fit$delta_forced[lead]),
          lead_n0 = unname(fit$genotype_counts[lead, "0"]),
          lead_n1 = unname(fit$genotype_counts[lead, "1"]),
          lead_n2 = unname(fit$genotype_counts[lead, "2"]), min_obs = fit$min_obs,
          cs_forced_mass = sum(alpha_cs * fit$delta_forced[cs]) / coverage),
          .slide_delta_summary(colSums(joint), grid, "cs_"),
          .slide_delta_summary(fit$alpha_delta[ell, lead, ], grid, "lead_"))
      }
    }
    list(rows = .slide_delta_bind(rows), error = NULL, status = "ok")
  }, error = function(e) list(rows = NULL, status = "error",
                             error = list(file = file, message = conditionMessage(e))))
}

summarize_slide_prior_delta <- function(path, workers = 1L, batch_size = 100L,
                                       min_entries = 3L, files = NULL,
                                       batch_dir = NULL, verbose = TRUE) {
  # min_entries = 3 preserves the user's length(out) > 2 filter.
  # Use min_entries = 0 to include files with only one or two tissue entries.
  # batch_dir optionally saves each completed batch (summary + errors), so an
  # interruption need not lose completed work. Existing batch files are never
  # silently reused: this function is not an automatic resume mechanism.
  for (name in c("workers", "batch_size", "min_entries")) {
    value <- get(name)
    lower <- if (name == "min_entries") 0 else 1
    if (length(value) != 1L || !is.finite(value) || value < lower ||
        value != floor(value)) stop("Invalid ", name)
  }
  allocation <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (!is.na(allocation) && workers > allocation)
    stop("workers exceeds SLURM_CPUS_PER_TASK (", allocation, ").")
  if (is.null(files)) files <- list.files(path, pattern = "\\.rds$",
                                         full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No RDS files found.")
  files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  batches <- split(seq_along(files), ceiling(seq_along(files) / batch_size))
  if (!is.null(batch_dir)) {
    dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(batch_dir)) stop("Cannot create batch_dir.")
    if (length(list.files(batch_dir, pattern = "^batch_[0-9]+\\.rds$")))
      stop("batch_dir already contains batches; use a new directory.")
  }
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makePSOCKcluster(min(workers, length(files)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c(".slide_delta_bind", ".slide_delta_summary",
                                 ".slide_delta_file"),
                            envir = environment(summarize_slide_prior_delta))
    parallel::clusterEvalQ(cl, {
      if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
      NULL
    })
  }
  summaries <- errors <- vector("list", length(batches))
  counts <- c(ok = 0L, filtered = 0L, error = 0L)
  started <- proc.time()[["elapsed"]]
  for (b in seq_along(batches)) {
    idx <- batches[[b]]
    ans <- if (is.null(cl)) lapply(files[idx], .slide_delta_file, min_entries) else
      parallel::parLapplyLB(cl, files[idx], .slide_delta_file,
                           min_entries = min_entries, chunk.size = 1L)
    summaries[[b]] <- .slide_delta_bind(lapply(ans, `[[`, "rows"))
    errors[[b]] <- .slide_delta_bind(lapply(ans, `[[`, "error"))
    status <- vapply(ans, `[[`, character(1), "status")
    counts <- counts + table(factor(status, levels = names(counts)))
    if (!is.null(batch_dir))
      saveRDS(list(files = files[idx], status = status,
                   results = summaries[[b]], errors = errors[[b]]),
              file.path(batch_dir, sprintf("batch_%05d.rds", b)))
    if (verbose) message(max(idx), "/", length(files), " files; ", counts["error"],
                         " failed; ", round(proc.time()[["elapsed"]] - started, 1), " s")
  }
  result <- list(results = .slide_delta_bind(summaries),
                 errors = .slide_delta_bind(errors), counts = counts)
  if (counts["error"] > 0L)
    warning(counts["error"], " file(s) failed. Inspect $errors before using $results.",
            call. = FALSE)
  result
}
