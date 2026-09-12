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

# Null components are uninformative, not evidence for a coding. Including
# rows equal to the prior leaves it unchanged. CS membership is never used.
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
regions <- lapply(seq_along(genotypes), function(g) {
  genotype <- genotypes[[g]]
  X <- scale(cbind(genotype, genotype == 2, genotype >= 1), scale = FALSE)
  y <- as.numeric(X[, if (g %% 2) 1 else 3]) + rnorm(nrow(X), sd = .4)
  configurations <- as.matrix(expand.grid(first = 1:3, second = 1:3))
  log_density <- apply(configurations, 1, function(configuration) {
    covariance <- diag(.4, nrow(X)) + .7 * tcrossprod(X[, configuration[1]]) +
      .3 * tcrossprod(X[, configuration[2]])
    upper <- chol(covariance)
    -.5 * (nrow(X) * log(2*pi) + 2 * sum(log(diag(upper))) +
             sum(forwardsolve(t(upper), y)^2))
  })
  list(configurations = configurations, log_density = log_density)
})
e_step <- function(prior) {
  counts <- matrix(0, 7, 3)
  likelihood <- 0
  for (region in regions) {
    configuration <- region$configurations
    log_weight <- region$log_density + log(prior[configuration[, 1]]) +
      log(prior[configuration[, 2]])
    normalizer <- max(log_weight) + log(sum(exp(log_weight - max(log_weight))))
    posterior <- exp(log_weight - normalizer)
    likelihood <- likelihood + normalizer
    for (coding in 1:3) counts[7, coding] <- counts[7, coding] +
      sum(posterior * rowSums(configuration == coding))
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
