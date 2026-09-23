# Hand-calculated analytic CI checks; base R only, without fitting models.
source("script/sim/simulation_metric_helpers.R")
near <- function(x, y) stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-12)))
expect_error <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))

# Three reported CSs, with sizes 1, 3, 5, spread unevenly across three runs.
d <- data.frame(scenario = "test", pve = .2, K = 2, configuration = "one",
  seed = 1:3, method = "SuSiE", n_cs = c(2, 0, 1), covered_cs = c(1, 0, 1),
  purity_sum = c(1.2, 0, .9), cs_size_sum = c(4, 0, 5),
  cs_size_sum_sq = c(10, 0, 25), recovered = c(1, 0, 2), n_causal = 2)
set.seed(123)
rng_before <- .Random.seed
result <- summarize_metrics(d)
stopifnot(identical(rng_before, .Random.seed))
s <- result$summary
z <- qnorm(.975)
near(s$coverage, 2/3); near(s$coverage_se, sqrt((2/3) * (1/3) / 3))
near(s$coverage_lo, max(0, 2/3 - z * s$coverage_se))
near(s$coverage_hi, 1)
near(s$purity, .7); near(s$purity_se, sqrt(.7 * .3 / 3))
near(s$power, .5); near(s$power_se, sqrt(.5 * .5 / 6))
near(s$power_lo, .5 - z * s$power_se)
near(s$power_hi, .5 + z * s$power_se)
near(s$cs_size, mean(c(1, 3, 5)))
near(s$cs_size_se, sd(c(1, 3, 5)) / sqrt(3))
near(s$cs_size_lo, 3 - z * 2 / sqrt(3))
near(s$cs_size_hi, 3 + z * 2 / sqrt(3))
stopifnot(s$coverage_n == 3, s$purity_n == 3, s$power_n == 6, s$cs_size_n == 3,
          !any(grepl("boot", names(s))))
alternative <- summarize_metrics(d, proportion_n = "replicates")$summary
near(alternative$power_se, sqrt(.5 * .5 / 3))

# Confidence level is configurable; uncertainty is not silently just +/-1 SE.
s90 <- summarize_metrics(d, level = .9)$summary
near(s90$power_hi, .5 + qnorm(.95) * s$power_se)

# No CSs supply no coverage, purity or size observations; power is still zero.
empty <- summarize_metrics(d[2, ])$summary
stopifnot(is.na(empty$coverage), is.na(empty$coverage_se), is.na(empty$coverage_lo),
          is.na(empty$purity_se), is.na(empty$cs_size_se),
          empty$power == 0, empty$power_se == 0, empty$power_hi == 0)
single <- summarize_metrics(d[3, ])$summary
stopifnot(single$coverage_se == 0, single$coverage_lo == 1,
          is.na(single$cs_size_se), is.na(single$cs_size_lo))
expect_error(summarize_metrics(d[, names(d) != "cs_size_sum_sq"]))

# CS summary deduplicates biological SNPs before taking size squares.
cs <- list(c(1L, 1L, 2L), 3L)
compact <- cs_summary(list(cs = cs, purity = list(min.abs.corr = c(.8, 1))), cs, c(1, 3))
stopifnot(compact["cs_size_sum"] == 3, compact["cs_size_sum_sq"] == 5,
          compact["covered_cs"] == 2)

# Paired Gaussian contrasts use the covariance from the shared simulations.
mix <- d; mix$method <- "SuSiE-mix"; mix$recovered <- c(2, 1, 2)
paired <- rbind(d, mix)
a <- summarize_metrics(paired)
delta <- a$differences[a$differences$metric == "power", ]
near(delta$difference, mean(c(.5, .5, 0)))
near(delta$standard_error, sd(c(.5, .5, 0)) / sqrt(3))
near(delta$lower, delta$difference - z * delta$standard_error)
# Reused seed blocks must not acquire artificial precision in paired contrasts.
copy <- paired; copy$configuration <- "two"
b <- summarize_metrics(rbind(paired, copy))
stopifnot(isTRUE(all.equal(a$differences, b$differences)))
expect_error(summarize_metrics(paired[-1, ]))
cat("PASS: denominator-based normal CIs, Gaussian CS-size CIs, empty/single-CS cases, and analytic paired contrasts.\n")
