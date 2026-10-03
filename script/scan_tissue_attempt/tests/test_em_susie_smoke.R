# Small real SuSiE fits, independent of GTEx. Uses an installed susieR package,
# or accepts a local source checkout: Rscript .../test_em_susie_smoke.R ../susieR
# For source-only dense-matrix checks, base R supplies colSds (normally from
# matrixStats); all fitting, variance updates, PIPs, and ELBO code are SuSiE's.
source("script/scan_tissue_attempt/em_utils.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  core <- new.env(parent = globalenv())
  core$colSds <- function(x) apply(x, 2, sd)
  for (name in c("initialize.R", "single_effect_regression.R", "sparse_multiplication.R",
                  "update_each_effect.R", "elbo.R", "estimate_residual_variance.R",
                  "susie_utils.R", "susie.R")) {
    sys.source(file.path(args[1], "R", name), core)
  }
  fit_susie <- core$susie
  cat("Testing SuSiE numerical code from:", normalizePath(args[1]), "\n")
} else {
  if (!requireNamespace("susieR", quietly = TRUE)) stop("Install susieR or provide its source checkout path.")
  fit_susie <- susieR::susie
}
set.seed(782)
coding <- rep(c("additive", "recessive", "dominant"), each = 5)
predictors <- paste0("snp", rep(1:5, 3), "__", coding)
regions <- lapply(1:4, function(i) {
  genotype <- matrix(rbinom(160*5, 2, .3), 160, 5)
  X <- cbind(genotype, genotype == 2, genotype >= 1)
  storage.mode(X) <- "double"
  colnames(X) <- predictors
  list(X = X, y = 1.8*X[, c(1, 7, 13, 2)[i]] + .6*X[, 5] + rnorm(160))
})
# A fifth region is orthogonal to every predictor and has no association.
# It must fit all-zero component variances and contribute no coding counts.
null_X <- regions[[1]]$X
regions[[5]] <- list(X = null_X, y = qr.resid(qr(cbind(1, null_X)), rnorm(nrow(null_X))))
common <- list(L = 4L, standardize = FALSE, estimate_prior_method = "EM",
               min_abs_corr = 0, max_iter = 1000, tol = 1e-7)
prior <- c(.4, .25, .35)
fits <- lapply(regions, function(region) do.call(fit_susie,
  c(region, common, list(prior_weights = em_predictor_weights(coding, setNames(prior, unique(coding)))))))
stopifnot(all(vapply(fits, function(fit) isTRUE(fit$converged), logical(1))))
stopifnot(all(fits[[5]]$V == 0), is.null(fits[[5]]$sets$cs))

# Exercise the production result-pooling path, not a separate test-only alpha
# summation. Persist full fitted posteriors just as cluster workers do.
test_dir <- tempfile("em_smoke_", tmpdir = "tmp")
dir.create(test_dir, recursive = TRUE)
files <- file.path(test_dir, paste0("gene", seq_along(fits), ".rds"))
pool <- function(fits, prior) {
  for (i in seq_along(fits)) saveRDS(list(Tissue = list(weighted_fit_mix = fits[[i]],
                                                       mix_coding = coding)), files[i])
  previous <- data.frame(tissue = "Tissue", pi_add = prior[1], pi_rec = prior[2], pi_dom = prior[3])
  pooled <- em_estimate_coding_priors(files, "weighted_fit_mix", previous)
  expected_active <- sum(vapply(fits, function(fit) sum(fit$V > 0), numeric(1)))
  stopifnot(pooled$priors$n_fits == 5L, pooled$priors$n_active_fits <= 4L,
            pooled$priors$n_active_components == expected_active,
            abs(pooled$priors$alpha_total - expected_active) < 1e-8,
            pooled$priors$n_zero_variance_components >= 4L,
            pooled$priors$n_source_elbo == 5L,
            abs(pooled$priors$source_elbo_sum -
                sum(vapply(fits, function(fit) tail(fit$elbo, 1), numeric(1)))) < 1e-8)
  pooled
}
objective <- sum(vapply(fits, function(fit) tail(fit$elbo, 1), numeric(1)))
initial <- objective
for (iteration in 1:10) {
  pooled <- pool(fits, prior)
  prior <- unname(em_prior_for_tissue(pooled$priors, "Tissue"))
  weights <- em_predictor_weights(coding, setNames(prior, unique(coding)), predictors)
  next_fits <- lapply(seq_along(fits), function(i) {
    initialization <- em_susie_initialization(fits[[i]], predictors, 4L, var(regions[[i]]$y))
    do.call(fit_susie, c(regions[[i]], common, list(prior_weights = weights), initialization))
  })
  stopifnot(all(vapply(next_fits, function(fit) isTRUE(fit$converged), logical(1))),
            all(vapply(next_fits, function(fit) max(abs(fit$pi - weights)) < 1e-12, logical(1))),
            all(vapply(next_fits, function(fit) nrow(fit$alpha) == 4L, logical(1))))
  stopifnot(all(next_fits[[5]]$V == 0))
  next_objective <- sum(vapply(next_fits, function(fit) tail(fit$elbo, 1), numeric(1)))
  stopifnot(next_objective >= objective - 1e-5)
  fits <- next_fits
  objective <- next_objective
}
stopifnot(objective > initial)
cat(sprintf("SuSiE variational EM: total ELBO %.6f -> %.6f; 10 nondecreasing updates.\n", initial, objective))
cat("Warm starts retained new coding weights and all four component rows.\n")
cat("Only positive-variance rows supplied prior counts; the all-null gene supplied none.\n")
unlink(files)
unlink(test_dir)
