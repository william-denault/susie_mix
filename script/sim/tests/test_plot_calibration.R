# Base-R regression checks; no refitting or production output changes.
source("script/sim/simulation_plot_helpers.R")
expect_error <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))

p <- c(0, .05, .1, .19, .5, .99, 1)
z <- calibration_bin_totals(p, c(2, 5, 7))
stopifnot(identical(z[1:10], c(2, 2, 0, 0, 0, 1, 0, 0, 0, 2)),
          sum(z[11:20]) == 3, abs(sum(z[21:30]) - sum(p)) < 1e-12)
expect_error(calibration_bin_totals(c(NA, .2), 1))
expect_error(calibration_bin_totals(c(.2, 1.1), 1))
expect_error(calibration_bin_totals(c(.2, .3), 3))

# With equal SNP counts per seed, the empirical ratio SE reduces to sd/sqrt(n).
state <- new.env(parent = emptyenv())
counts <- t(vapply(c(1, 3, 5), function(k)
  calibration_bin_totals(rep(.55, 5), seq_len(k)), numeric(31)))
meta <- data.frame(scenario = "Additive only", pve = .2, K = 1, method = "SuSiE")
calibration_merge(state, meta, 1:3, counts)
s <- calibration_summary(as.list(state))
bin <- s[s$bin == 6, ]
stopifnot(bin$n_snps == 15, bin$n_causal == 9, bin$mean_pip == .55,
          abs(bin$frequency - .6) < 1e-12,
          abs(bin$empirical_se - sd(c(.2, .6, 1)) / sqrt(3)) < 1e-12,
          all(is.na(s$frequency[s$n_snps == 0])))

# Duplicate conditions sharing the same seeds do not create independent draws.
meta$K <- 2
calibration_merge(state, meta, 1:3, counts)
pooled <- calibration_summary(as.list(state), pool_K = TRUE)
bin2 <- pooled[pooled$bin == 6, ]
stopifnot(bin2$n_snps == 30, bin2$n_causal == 18,
          bin2$n_seed_blocks == 3, bin2$n_replicates == 6,
          abs(bin2$empirical_se - bin$empirical_se) < 1e-12)

# Unequal bin occupancy: pool counts, not the per-replicate frequencies.
unequal <- new.env(parent = emptyenv())
meta$K <- 1
v <- rbind(calibration_bin_totals(c(.51, .01), 1),
           calibration_bin_totals(c(rep(.59, 9), .01), 10))
calibration_merge(unequal, meta, c(1, 2), v)
u <- calibration_summary(as.list(unequal))
u <- u[u$bin == 6, ]
stopifnot(u$n_snps == 10, u$n_causal == 1, u$frequency == .1,
          abs(u$mean_pip - .582) < 1e-12)

d <- data.frame(coverage = c(.91, .98), coverage_lo = c(.2, .3),
                purity = c(.84, .94), purity_lo = c(.1, .2),
                cs_size = c(7, 21), cs_size_hi = c(10, 40))
stopifnot(identical(panel_y_limits("coverage", d), c(.91, 1)),
          identical(panel_y_limits("purity", d), c(.84, 1)),
          identical(panel_y_limits("cs_size", d), c(0, 21)),
          identical(panel_y_limits("power", d), c(0, 1)),
          identical(panel_y_limits("coverage", d[FALSE, ]), c(0, 1)),
          diff(panel_y_limits("coverage", data.frame(coverage = 1))) > 0)
curve <- data.frame(method = "SuSiE", threshold = c(Inf, .1), fdr = c(0, .5), tpr = c(0, 1))
stopifnot(identical(panel_y_limits("power_fdr", curve, c(0, .25)), c(0, .5)))
# Nonmonotone FDR must preserve threshold ordering when clipping line segments.
curve <- rbind(curve, data.frame(method = "SuSiE", threshold = 0, fdr = .1, tpr = .8))
stopifnot(abs(visible_curve_max(curve, c(0, .25)) - .875) < 1e-12)

selection <- data.frame(scenario = "Additive only", pve = .2, K = 1,
                        configuration = "add1", seed = c(2, 1), method = "SuSiE")
stopifnot(identical(calibration_selection(selection, "SuSiE")$seed, c(1, 2)))
cat("PASS: calibration bins, pooling, seed-level SE, empty bins, and panel limits.\n")
