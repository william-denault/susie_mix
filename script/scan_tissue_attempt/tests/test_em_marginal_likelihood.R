# Analytic and enumerated-model checks; base R only, no GTEx or Slurm.
# Run from the project root with Rscript --vanilla.
source("script/scan_tissue_attempt/em_utils.R")
equal <- function(x, y, tolerance = 1e-7) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
}

# Two components can select the same predictor. Marginal PIPs cap its support
# at one inclusion; alpha counts both assignments, as the categorical Q needs.
alpha <- rbind(c(.8, .1, .1), c(.8, .1, .1))
counts <- matrix(0, 7, 3)
counts[7, ] <- colSums(alpha)
update <- em_coding_mstep(counts)
equal(update$prior, c(.8, .1, .1))
pip <- 1 - apply(1 - alpha, 2, prod)
stopifnot(max(abs(pip/sum(pip) - update$prior)) > .08,
          em_coding_q(update$prior, counts) > em_coding_q(pip/sum(pip), counts))

# Conditional coding availability changes the denominator in the objective.
# A:R posterior odds are 4:1 and R:D odds are 1:4, giving A:R:D = 4:1:4.
counts <- matrix(0, 7, 3)
counts[3, ] <- c(8, 2, 0)
counts[6, ] <- c(0, 2, 8)
update <- em_coding_mstep(counts)
equal(update$prior, c(4, 1, 4)/9, 1e-5)
stopifnot(em_coding_q(update$prior, counts) > em_coding_q(colSums(counts)/sum(counts), counts))

# Singleton-only codings carry no relative-prior information. Keep old mass
# between disconnected availability groups, instead of inventing enrichment.
counts <- matrix(0, 7, 3)
counts[3, ] <- c(8, 2, 0)
counts[4, ] <- c(0, 0, 100)
update <- em_coding_mstep(counts, c(.2, .3, .5))
equal(update$prior, c(.4, .1, .5), 1e-5)
counts <- matrix(0, 7, 3)
counts[cbind(c(1, 2, 4), 1:3)] <- c(100, 1, 50)
equal(em_coding_mstep(counts, c(.2, .3, .5))$prior, c(.2, .3, .5))

# With no active counts, retain the prior without adding pseudocounts.
counts <- matrix(0, 7, 3)
null_update <- em_coding_mstep(counts, c(.2, .3, .5))
equal(null_update$prior, c(.2, .3, .5))
stopifnot(null_update$gain == 0, null_update$status == "carried_forward_no_active_components")
# For comparison, the earlier uncollapsed update also leaves a null-only
# prior unchanged when its alpha rows equal the prior.
counts <- matrix(0, 7, 3)
counts[7, ] <- 100 * c(.2, .3, .5)
equal(em_coding_mstep(counts, c(.2, .3, .5))$prior, c(.2, .3, .5))

# Test complete EM cycles against the exact marginal likelihood of a Gaussian
# two-effect model. Enumerate all p^2 assignments, including repeated indices;
# integrate both normally distributed effect sizes analytically. This checks
# the target likelihood independently of the code's Q calculation.
set.seed(91)
genotypes <- list(c(0, 0, 0, 1, 1, 1, 2, 2), c(0, 0, 1, 1, 1, 2, 2, 2),
                  c(0, 0, 0, 0, 1, 1, 2, 2), c(0, 0, 1, 1, 2, 2, 2, 2))
make_region <- function(X, y, V) {
  configurations <- as.matrix(expand.grid(rep(list(1:3), length(V))))
  log_density <- apply(configurations, 1, function(configuration) {
    covariance <- diag(.4, nrow(X))
    for (l in seq_along(V)) covariance <- covariance + V[l] * tcrossprod(X[, configuration[l]])
    upper <- chol(covariance)
    -.5 * (nrow(X) * log(2*pi) + 2 * sum(log(diag(upper))) +
             sum(forwardsolve(t(upper), y)^2))
  })
  list(configurations = configurations, log_density = log_density, V = V, X = X, y = y)
}
regions <- lapply(seq_along(genotypes), function(g) {
  genotype <- genotypes[[g]]
  X <- scale(cbind(genotype, genotype == 2, genotype >= 1), scale = FALSE)
  y <- as.numeric(X[, if (g %% 2) 1 else 3]) + rnorm(nrow(X), sd = .4)
  make_region(X, y, c(.7, .3))
})
e_step <- function(prior, fitted_regions = regions) {
  counts <- matrix(0, 7, 3)
  likelihood <- 0
  for (region in fitted_regions) {
    configuration <- region$configurations
    # Include every latent assignment, even inactive ones, in the independent
    # exact marginal-likelihood calculation. Only the M-step counts collapse.
    log_weight <- region$log_density + rowSums(matrix(log(prior[configuration]), nrow(configuration)))
    normalizer <- max(log_weight) + log(sum(exp(log_weight - max(log_weight))))
    posterior <- exp(log_weight - normalizer)
    likelihood <- likelihood + normalizer
    for (coding in 1:3) counts[7, coding] <- counts[7, coding] +
      sum(posterior * rowSums(configuration[, region$V > 0, drop = FALSE] == coding))
  }
  list(counts = counts, likelihood = likelihood)
}
prior <- c(.2, .5, .3)
state <- e_step(prior)
initial_likelihood <- state$likelihood
for (iteration in 1:50) {
  update <- em_coding_mstep(state$counts, prior)
  next_state <- e_step(update$prior)
  stopifnot(next_state$likelihood >= state$likelihood - 1e-10)
  prior <- update$prior
  state <- next_state
}
stopifnot(state$likelihood > initial_likelihood + .1)
cat(sprintf("Exact two-effect Gaussian EM: log marginal likelihood %.6f -> %.6f; 50 nondecreasing updates.\n",
            initial_likelihood, state$likelihood))

# Appending an exactly inactive slot preserves each marginal likelihood and
# its active assignment counts, for arbitrary coding priors. An all-null
# gene contributes a constant likelihood and zero M-step counts.
with_null <- lapply(regions, function(region) make_region(region$X, region$y, c(region$V, 0)))
null_gene <- make_region(regions[[1]]$X, regions[[1]]$y, c(0, 0, 0))
for (p in list(c(.2, .5, .3), c(.7, .1, .2), c(.05, .9, .05))) {
  equal(e_step(p, with_null), e_step(p))
  equal(e_step(p, list(null_gene))$counts, matrix(0, 7, 3))
  equal(e_step(p, list(null_gene))$likelihood,
        e_step(rep(1/3, 3), list(null_gene))$likelihood)
}
mixed_regions <- c(with_null, list(null_gene))
prior <- c(.2, .5, .3)
state <- e_step(prior, mixed_regions)
initial_likelihood <- state$likelihood
for (iteration in 1:50) {
  update <- em_coding_mstep(state$counts, prior)
  next_state <- e_step(update$prior, mixed_regions)
  stopifnot(next_state$likelihood >= state$likelihood - 1e-10)
  prior <- update$prior
  state <- next_state
}
stopifnot(state$likelihood > initial_likelihood + .1)
cat(sprintf("Collapsed EM with zero-variance slots and a null gene: log marginal likelihood %.6f -> %.6f; 50 nondecreasing updates.\n",
            initial_likelihood, state$likelihood))

# Warm starts preserve the posterior, residual variance, and all component
# variances, while omitting the old pi that could overwrite the new weights.
predictor <- c("snp__additive", "snp__recessive", "snp__dominant")
colnames(alpha) <- predictor
fit <- list(alpha = alpha, mu = alpha * 0, mu2 = alpha * 0 + .1,
            V = c(.5, 0), sigma2 = .8, pi = c(.9, .05, .05))
init <- em_susie_initialization(fit, predictor, 2, 2)
equal(init$s_init$alpha, alpha)
equal(init$scaled_prior_variance, c(.25, 0))
equal(init$residual_variance, .8)
stopifnot(!"pi" %in% names(init$s_init), !"V" %in% names(init$s_init))
bad <- tryCatch(em_susie_initialization(fit, rev(predictor), 2, 2), error = identity)
stopifnot(inherits(bad, "error"))
cat("All EM marginal-likelihood and prior-update checks passed.\n")
