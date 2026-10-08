# Global, tissue-specific EM for the discrete slider prior. Base R helpers.
# Source em_utils.R first for atomic writes, manifests and em_stop_fit.
spe_grid <- function() seq(-1, 1, length.out = 17L)
spe_columns <- function() sprintf("w_%02d", seq_along(spe_grid()))
spe_settings <- function() list(
  schema = 1L, delta_grid = spe_grid(), L = 10L, standardize = FALSE,
  estimate_prior_method = "EM", min_obs = 5L, min_abs_corr = 0.5,
  max_iter = 1000L, tol = 1e-5, min_samples = 50L, min_n_rec = 5L,
  min_maf = 0.05, min_maf_plink = 0, hwe_thresh = 1e-8,
  cis_window = 5e5, seed = 1L)

# A larger numerical iteration budget is allowed on retry; model parameters,
# convergence tolerance, cohort and the frozen prior must remain unchanged.
spe_fit_budget <- function(saved_max_iter) {
  value <- Sys.getenv("SUSIE_SLIDE_MAX_ITER", "")
  if (!nzchar(value)) return(saved_max_iter)
  budget <- suppressWarnings(as.numeric(value))
  if (length(budget) != 1L || !is.finite(budget) || budget != floor(budget) ||
      budget < saved_max_iter || budget > .Machine$integer.max)
    stop("SUSIE_SLIDE_MAX_ITER must be an integer at least ", saved_max_iter, ".")
  as.integer(budget)
}

spe_validate_priors <- function(priors) {
  cols <- spe_columns()
  if (!is.data.frame(priors) || !nrow(priors) ||
      !all(c("tissue", cols) %in% names(priors)) ||
      anyNA(priors$tissue) || any(!nzchar(priors$tissue)) ||
      anyDuplicated(priors$tissue) ||
      !all(vapply(priors[cols], is.numeric, logical(1))))
    stop("Invalid slider tissue prior table.")
  weights <- as.matrix(priors[cols])
  if (any(!is.finite(weights)) || any(weights < 0) ||
      any(abs(rowSums(weights) - 1) > 1e-8))
    stop("Slider weights must be finite, nonnegative and sum to one per tissue.")
  invisible(priors)
}

spe_prior <- function(priors, tissue) {
  i <- match(tissue, priors$tissue)
  if (is.na(i)) stop("No slider prior for tissue: ", tissue)
  as.numeric(priors[i, spe_columns()])
}

spe_fit_counts <- function(fit, prior, require_converged = TRUE) {
  if (!inherits(fit, "susie_slide") ||
      !isTRUE(all.equal(as.numeric(fit$delta_grid), spe_grid(), tolerance = 1e-12)) ||
      !isTRUE(all.equal(as.numeric(fit$delta_prior), prior, tolerance = 1e-10)))
    stop("Expected a discrete slider fit with the frozen grid and prior.")
  if (require_converged && !isTRUE(fit$converged)) stop("Unconverged slider fit.")
  a <- fit$alpha
  counts <- fit$delta_prior_counts
  if (!is.matrix(a) || !is.numeric(a) || any(!is.finite(a)) || any(a < 0) ||
      any(abs(rowSums(a) - 1) > 1e-7)) stop("Invalid slider alpha matrix.")
  L <- nrow(a); p <- ncol(a)
  V <- fit$V
  if (length(V) == 1L) V <- rep(V, L)
  if (length(V) != L || any(!is.finite(V)) || any(V < 0)) stop("Invalid slider V.")
  if (!is.matrix(counts) || !is.numeric(counts) ||
      !identical(dim(counts), c(L, length(spe_grid()))) ||
      any(!is.finite(counts)) || any(counts < 0))
    stop("Missing or invalid delta_prior_counts; use susieRSlidePrior >= 0.3.0.")
  if (!is.logical(fit$delta_forced) || length(fit$delta_forced) != p ||
      anyNA(fit$delta_forced) || length(fit$input_p) != 1L ||
      !is.finite(fit$input_p) || fit$input_p < 1 || fit$input_p > p ||
      fit$input_p != floor(fit$input_p)) stop("Invalid forced-SNP or null metadata.")
  if (!identical(dim(fit$alpha_delta), c(L, p, length(spe_grid()))))
    stop("Missing joint SNP/slider posterior for warm starts.")
  free <- !fit$delta_forced & seq_len(p) <= fit$input_p
  # A zero-variance component has no dependence on the slider prior after
  # integrating its assignment out. Forced-additive SNPs and null columns
  # likewise supply no slider-prior counts. No CS eligibility is consulted.
  counts[V == 0, ] <- 0
  expected <- rowSums(a[, free, drop = FALSE])
  expected[V == 0] <- 0
  if (any(abs(rowSums(counts) - expected) > 1e-7))
    stop("Slider counts do not match free-SNP posterior mass.")
  if (any(counts[, prior == 0, drop = FALSE] > 1e-12))
    stop("Posterior mass lies outside the frozen slider prior support.")
  counts
}

spe_update_weights <- function(counts, previous) {
  counts <- as.numeric(counts)
  if (length(counts) != length(previous) || any(!is.finite(counts)) || any(counts < 0))
    stop("Invalid pooled slider counts.")
  if (sum(counts) == 0) return(previous)
  counts / sum(counts)
}

spe_pool <- function(files, previous) {
  spe_validate_priors(previous)
  totals <- matrix(0, nrow(previous), length(spe_grid()),
                   dimnames = list(previous$tissue, spe_columns()))
  audit <- list(); n <- 0L
  add_audit <- function(gene, tissue = "", issue, active = 0L, mass = 0, elbo = NA_real_, detail = "") {
    n <<- n + 1L
    audit[[n]] <<- data.frame(gene, tissue, issue, active, free_mass = mass, elbo, detail)
  }
  for (file in files) {
    out <- readRDS(file)
    gene <- sub("\\.rds$", "", basename(file))
    if (!is.list(out)) stop("Invalid gene result: ", file)
    if (!is.null(out$error)) {
      add_audit(gene, issue = "gene_error", detail = as.character(out$error))
      next
    }
    errors <- attr(out, "tissue_errors", exact = TRUE)
    for (tissue in names(errors)) add_audit(gene, tissue, "tissue_error", detail = errors[[tissue]])
    if (length(out) && (is.null(names(out)) || anyDuplicated(names(out)) ||
                        any(!names(out) %in% previous$tissue))) stop("Invalid tissue names: ", file)
    for (tissue in names(out)) {
      fit <- out[[tissue]]$fit_slide_prior
      counts <- tryCatch(spe_fit_counts(fit, spe_prior(previous, tissue)),
                         error = function(e) stop(gene, " / ", tissue, ": ", conditionMessage(e)))
      totals[tissue, ] <- totals[tissue, ] + colSums(counts)
      if (!length(fit$elbo) || !is.finite(tail(fit$elbo, 1))) stop("Missing finite ELBO: ", file)
      active <- sum(rep(fit$V, length.out = nrow(counts)) > 0)
      add_audit(gene, tissue, "fit", active, sum(counts), tail(fit$elbo, 1))
    }
    if (!length(out) && !length(errors)) add_audit(gene, issue = "no_eligible_tissues")
  }
  audit <- do.call(rbind, audit)
  if (is.null(audit) || !any(audit$issue == "fit")) stop("No successful slider fits to learn from.")
  weights <- t(vapply(seq_len(nrow(previous)), function(i)
    spe_update_weights(totals[i, ], spe_prior(previous, previous$tissue[i])), numeric(17)))
  colnames(weights) <- spe_columns()
  diagnostics <- lapply(seq_len(nrow(previous)), function(i) {
    old <- spe_prior(previous, previous$tissue[i]); new <- weights[i, ]
    used <- totals[i, ] > 0
    rows <- audit$tissue == previous$tissue[i] & audit$issue == "fit"
    gain <- sum(totals[i, used] * (log(new[used]) - log(old[used])))
    if (!is.finite(gain) || gain < -1e-7) stop("Slider M-step failed its Q-increase check.")
    data.frame(n_fits = sum(rows), n_active = sum(audit$active[rows]),
               free_mass = sum(totals[i, ]), prior_max_change = max(abs(new - old)),
               m_step_q_gain = gain, source_elbo = sum(audit$elbo[rows]))
  })
  list(priors = data.frame(tissue = previous$tissue, weights, do.call(rbind, diagnostics)),
       audit = audit, counts = data.frame(tissue = previous$tissue, totals))
}

# Diagnostics for the fits at iteration k, before applying the k -> k+1 M-step.
# The total is a sum of per-gene/tissue variational objectives (a composite
# objective for overlapping genes/tissues), not the exact joint log evidence.
spe_record_objective <- function(project, iteration_dir, pooled) {
  root <- file.path(project, "results_slide_prior_em")
  priors <- read.csv(file.path(iteration_dir, "priors.csv"), stringsAsFactors = FALSE)
  iteration <- unique(priors$iteration)
  if (length(iteration) != 1L || !is.finite(iteration) || iteration < 0)
    stop("Invalid objective iteration.")
  fit_rows <- function(audit) {
    x <- audit[audit$issue == "fit", c("gene", "tissue", "elbo"), drop = FALSE]
    if (!nrow(x) || anyNA(x) || anyDuplicated(x[c("gene", "tissue")]) ||
        any(!is.finite(x$elbo))) stop("Invalid objective fit audit.")
    x[order(x$gene, x$tissue, method = "radix"), , drop = FALSE]
  }
  fits <- fit_rows(pooled$audit)
  previous_fits <- previous_priors <- NULL
  if (iteration > 0L) {
    previous_dir <- file.path(root, sprintf("iteration_%03d", iteration - 1L))
    previous_file <- file.path(previous_dir, "objective_fit_audit.csv")
    if (!file.exists(previous_file))
      stop("Verify iteration ", iteration - 1L, " first to record its objective baseline.")
    previous_fits <- fit_rows(read.csv(previous_file, stringsAsFactors = FALSE))
    if (!identical(unname(as.matrix(fits[c("gene", "tissue")])),
                   unname(as.matrix(previous_fits[c("gene", "tissue")]))))
      stop("Objective gene/tissue cohort changed; iteration totals are not comparable.")
    previous_priors <- read.csv(file.path(previous_dir, "priors.csv"), stringsAsFactors = FALSE)
    if (!setequal(previous_priors$tissue, priors$tissue)) stop("Objective tissue universe changed.")
  }
  changes <- vapply(priors$tissue, function(tissue) {
    if (is.null(previous_priors)) return(NA_real_)
    max(abs(spe_prior(priors, tissue) - spe_prior(previous_priors, tissue)))
  }, numeric(1))
  summarize <- function(tissue = NULL) {
    rows <- if (is.null(tissue)) rep(TRUE, nrow(fits)) else fits$tissue == tissue
    dx <- if (is.null(tissue)) seq_len(nrow(pooled$priors)) else match(tissue, pooled$priors$tissue)
    nf <- sum(rows)
    value <- sum(fits$elbo[rows])
    previous <- if (is.null(previous_fits)) NA_real_ else sum(previous_fits$elbo[rows])
    change <- value - previous
    # This is a floating-point warning threshold, not a convergence criterion.
    roundoff <- if (is.na(previous)) NA_real_ else max(1e-6, 100 * .Machine$double.eps * abs(previous))
    data.frame(iteration = iteration, level = if (is.null(tissue)) "overall" else "tissue",
      tissue = if (is.null(tissue)) "" else tissue, n_genes = length(unique(fits$gene[rows])),
      n_fits = nf, elbo = value, delta_elbo = change,
      relative_delta_elbo = change / max(1, abs(previous)),
      delta_elbo_per_fit = if (nf) change / nf else NA_real_,
      prior_max_change = if (is.null(tissue)) max(changes) else unname(changes[tissue]),
      next_m_step_q_gain = sum(pooled$priors$m_step_q_gain[dx]),
      next_prior_max_change = max(pooled$priors$prior_max_change[dx]),
      elbo_decreased = if (is.na(change)) NA else change < -roundoff,
      recorded_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  }
  current <- do.call(rbind, c(list(summarize()), lapply(priors$tissue, summarize)))
  history_file <- file.path(root, "objective_history.csv")
  history <- if (file.exists(history_file)) read.csv(history_file, stringsAsFactors = FALSE) else NULL
  if (!is.null(history) && any(history$iteration > iteration)) {
    old <- history[history$iteration == iteration, , drop = FALSE]
    cols <- setdiff(names(current), "recorded_at")
    old <- old[order(old$level, old$tissue), cols, drop = FALSE]
    check <- current[order(current$level, current$tissue), cols, drop = FALSE]
    if (!isTRUE(all.equal(old, check, check.attributes = FALSE, tolerance = 1e-12)))
      stop("Cannot change an earlier objective after later iterations were recorded.")
  }
  if (!is.null(history)) history <- history[history$iteration != iteration, , drop = FALSE]
  history <- rbind(history, current)
  history <- history[order(history$iteration, history$level, history$tissue), , drop = FALSE]
  # Reverification replaces an iteration's rows, so restarts never double-count.
  em_atomic_write(pooled$audit, file.path(iteration_dir, "objective_fit_audit.csv"), csv = TRUE)
  em_atomic_write(current, file.path(iteration_dir, "objective_summary.csv"), csv = TRUE)
  em_atomic_write(history, history_file, csv = TRUE)
  if (any(current$elbo_decreased %in% TRUE))
    warning("ELBO decreased for iteration ", iteration, "; inspect objective_history.csv.", call. = FALSE)
  message("Recorded iteration ", iteration, " objective: ", format(current$elbo[1], digits = 15),
          " across ", current$n_fits[1], " fits.")
  invisible(current)
}

spe_check_packages <- function() {
  for (package in c("susieRSlidePrior", "data.table", "matrixStats"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package is missing: ", package)
  if (utils::packageVersion("susieRSlidePrior") < "0.3.0" ||
      !all(c("delta_grid", "delta_prior") %in% names(formals(susieRSlidePrior::susie))))
    stop("Install susieRSlidePrior >= 0.3.0 with discrete-prior counts and warm starts.")
  as.character(utils::packageVersion("susieRSlidePrior"))
}
