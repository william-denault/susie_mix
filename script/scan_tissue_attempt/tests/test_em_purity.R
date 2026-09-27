# Rscript --vanilla script/scan_tissue_attempt/tests/test_em_purity.R
# Deterministic CS geometry plus the production EB pooling path.
source("script/scan_tissue_attempt/em_utils.R")
stopifnot(requireNamespace("susieR", quietly = TRUE))
equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

run_tests <- function() {
  set.seed(61)
  orthogonal <- qr.Q(qr(cbind(1, matrix(rnorm(320), 80, 4))))[, -1]
  X <- cbind(orthogonal[, 1], .8*orthogonal[, 1] + .6*orthogonal[, 2],
             .8*orthogonal[, 1] + .6*orthogonal[, 3], orthogonal[, 4])
  colnames(X) <- paste0("snp", 1:4)
  alpha <- rbind(c(.5, .49, .005, .005), c(.5, .49, .005, .005),
                 c(.005, .5, .49, .005), c(.5, .005, .005, .49),
                 c(.005, .005, .985, .005), c(.99, .005, .0025, .0025))
  colnames(alpha) <- colnames(X)
  fit <- structure(list(alpha = alpha, V = c(1, 1, 1, 1, 1e-12, 0),
    pip = setNames(1 - apply(1 - alpha[1:5, ], 2, prod), colnames(X)),
    mu = alpha*0, mu2 = alpha*0 + .1, sigma2 = 1,
    pi = rep(.25, 4), converged = TRUE), class = "susie")
  fit$sets <- susieR::susie_get_cs(fit, X = X, min_abs_corr = .5)
  original <- fit
  fit <- em_attach_cs_eligibility(fit, X)
  # L1/L2 have identical sets, L3 overlaps them, L4 fails purity, L5 is a
  # singleton with tiny positive V, L6 has zero V. Purity sorting reorders Ls.
  stopifnot(identical(em_eligible_components(fit), c(5L, 1L, 2L, 3L)),
            length(fit$sets$cs) == 2L)
  equal(fit$coding_prior_cs$components$min_abs_corr, c(1, .8, .8, .64))
  restored <- fit
  restored$coding_prior_cs <- NULL
  stopifnot(identical(restored, original))
  null <- fit
  null$alpha[,] <- .25
  null$V <- 1
  null$sets <- susieR::susie_get_cs(null, X = X, min_abs_corr = .5)
  null <- em_attach_cs_eligibility(null, X)
  stopifnot(is.null(null$sets$cs), length(em_eligible_components(null)) == 0L,
            all(null$V > 0), nrow(null$alpha) == 6L)

  dir.create("tmp", showWarnings = FALSE)
  test_base <- normalizePath("tmp", winslash = "/", mustWork = TRUE)
  project <- tempfile("em_purity_", tmpdir = test_base)
  dir.create(project)
  on.exit({
    resolved <- normalizePath(project, winslash = "/", mustWork = TRUE)
    stopifnot(startsWith(resolved, paste0(test_base, "/em_purity_")))
    unlink(resolved, recursive = TRUE)
  }, add = TRUE)
  coding <- c("additive", "recessive", "dominant", "additive")
  save_fit <- function(fit, name) {
    file <- file.path(project, paste0(name, ".rds"))
    saveRDS(list(Tissue = list(weighted_fit_mix = fit, mix_coding = coding)), file)
    file
  }
  signal_file <- save_fit(fit, "signal")
  null_file <- save_fit(null, "null")
  previous <- data.frame(tissue = "Tissue", pi_add = .6, pi_rec = .3, pi_dom = .1)
  pooled <- em_estimate_coding_priors(c(signal_file, null_file), "weighted_fit_mix", previous)
  expected <- colSums(alpha[c(1, 2, 3, 5), ])/4
  equal(em_prior_for_tissue(pooled$priors, "Tissue"),
        c(sum(expected[c(1, 4)]), expected[2], expected[3]))
  stopifnot(pooled$priors$n_fits == 2L, pooled$priors$n_active_components == 11L,
            pooled$priors$n_eligible_components == 4L, pooled$priors$alpha_total == 4,
            pooled$priors$n_eligible_fits == 1L, pooled$priors$n_zero_cs_fits == 1L,
            nrow(readRDS(signal_file)$Tissue$weighted_fit_mix$alpha) == 6L)
  carried <- em_estimate_coding_priors(null_file, "weighted_fit_mix", previous)$priors
  equal(em_prior_for_tissue(carried, "Tissue"), c(.6, .3, .1))
  stopifnot(carried$alpha_total == 0, carried$mstep_q_gain == 0,
            carried$prior_status == "carried_forward_no_eligible_components")
  initial <- em_estimate_coding_priors(null_file, "weighted_fit_mix")$priors
  equal(em_prior_for_tissue(initial, "Tissue"), rep(1/3, 3))
  stopifnot(initial$prior_status == "initialized_uniform_no_eligible_components")

  old <- fit
  old$coding_prior_cs <- NULL
  expect_error(em_eligible_components(old), "Missing per-component")
  stale <- fit
  stale$V[1] <- 0
  expect_error(em_eligible_components(stale), "Invalid per-component")
  bad <- fit
  bad$coding_prior_cs$dedup <- TRUE
  expect_error(em_eligible_components(bad), "Invalid per-component")
  expect_error(em_attach_cs_eligibility(fit, X[, 4:1], .5), "predictor order mismatch")
  expect_error(em_attach_cs_eligibility(fit, X, 0), "required min_abs_corr")
  cat("Purity eligibility, duplicate/overlapping CSs, zero-CS pooling, and prior fallback passed.\n")
}
run_tests()
