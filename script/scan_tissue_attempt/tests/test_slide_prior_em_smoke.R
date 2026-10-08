# Small real fits; needs the same susieRSlidePrior package as RCC.
extra_lib <- Sys.getenv("SLIDER_PRIOR_TEST_LIB")
if (nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))
for (file in c("em_utils.R", "slide_prior_em_utils.R")) source(file.path("script/scan_tissue_attempt", file))
spe_check_packages()
set.seed(374)
datasets <- lapply(1:4, function(i) {
  X <- matrix(rbinom(240 * 12, 2, .35), 240, 12)
  X[, 12] <- 0; X[1:8, 12] <- 1 # No homozygous alternate: forced additive.
  colnames(X) <- paste0("snp", 1:12)
  delta <- c(-1, -.5, .5, 1)[i]
  list(X = X, y = 1.3 * (X[, 1] + delta * (X[, 1] == 1)) + rnorm(nrow(X)))
})
prior <- rep(1 / 17, 17)
fits <- vector("list", length(datasets))
objectives <- numeric()
repo <- normalizePath(".", winslash = "/")
dir.create(file.path(repo, "tmp"), showWarnings = FALSE)
project <- tempfile("slide_objective_smoke_", tmpdir = file.path(repo, "tmp")); dir.create(project)
for (iteration in 0:5) {
  fits <- lapply(seq_along(datasets), function(i) {
    d <- datasets[[i]]
    args <- list(X = d$X, y = d$y, L = 1, min_obs = 5, standardize = FALSE,
                 estimate_prior_method = "EM", delta_grid = spe_grid(), delta_prior = prior,
                 max_iter = 2000, tol = 1e-7, coverage = NULL)
    if (iteration > 0) args$model_init <- fits[[i]]
    do.call(susieRSlidePrior::susie, args)
  })
  for (f in fits) {
    spe_fit_counts(f, prior)
    stopifnot(f$delta_forced[12])
    free <- which(!f$delta_forced & seq_len(ncol(f$alpha)) <= f$input_p)
    direct <- apply(f$alpha_delta[, free, , drop = FALSE], c(1, 3), sum)
    direct[f$V == 0, ] <- 0
    stopifnot(isTRUE(all.equal(unname(direct), unname(spe_fit_counts(f, prior)), tolerance = 1e-8)))
  }
  objectives <- c(objectives, sum(vapply(fits, function(f) tail(f$elbo, 1), numeric(1))))
  d <- file.path(project, "results_slide_prior_em", sprintf("iteration_%03d", iteration))
  dir.create(file.path(d, "results"), recursive = TRUE)
  table <- data.frame(iteration = iteration, tissue = "Fixture", as.list(setNames(prior, spe_columns())))
  write.csv(table, file.path(d, "priors.csv"), row.names = FALSE)
  for (i in seq_along(fits)) saveRDS(list(Fixture = list(fit_slide_prior = fits[[i]])),
                                    file.path(d, "results", paste0("G", i, ".rds")))
  pooled <- spe_pool(list.files(file.path(d, "results"), full.names = TRUE), table)
  recorded <- spe_record_objective(project, d, pooled)
  stopifnot(abs(recorded$elbo[1] - tail(objectives, 1)) < 1e-8)
  totals <- Reduce(`+`, lapply(fits, function(f) colSums(spe_fit_counts(f, prior))))
  next_prior <- spe_update_weights(totals, prior)
  used <- totals > 0
  stopifnot(sum(totals[used] * log(next_prior[used] / prior[used])) >= -1e-8)
  prior <- next_prior
}
stopifnot(all(diff(objectives) >= -1e-5), max(abs(prior - 1 / 17)) > .001)
print(data.frame(iteration = 0:5, elbo = objectives))
history <- read.csv(file.path(project, "results_slide_prior_em/objective_history.csv"))
stopifnot(nrow(history) == 12L, identical(history$iteration[history$level == "overall"], 0:5),
          identical(dirname(normalizePath(project, winslash = "/")), paste0(repo, "/tmp")))
unlink(project, recursive = TRUE)
cat("PASS: real 17-point slider fits, forced SNPs, full warm starts, recorded objectives and five EM updates with nondecreasing ELBO.\n")
