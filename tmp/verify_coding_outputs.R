options(susie_mix.coding_plots.run = FALSE)
source("script/sim/plot_coding_mixture.R")
d <- "simulation results/coding_figures"
config <- readRDS(file.path(d, "plot_settings.rds"))
r <- readRDS(file.path(d, "replicate_coding_metrics.rds"))
a <- read.csv(file.path(d, "file_audit.csv"))
cal <- read.csv(file.path(d, "coding_pip_calibration.csv"))
calls <- read.csv(file.path(d, "coding_discovery_accuracy.csv"))
m <- read.csv(file.path(d, "mixture_recovery.csv"))
stopifnot(sum(a$included) == nrow(r), sum(cal$n_true) == sum(r$K),
          sum(calls$n_correct[calls$threshold == 0]) == sum(r$K), sum(a$invalid) == 0,
          all(abs(calls$n_calls - calls$n_correct - calls$n_wrong_coding - calls$n_wrong_snp) < 1e-8),
          isTRUE(all.equal(sum(cal$sum_pip), sum(r[c("mass_add", "mass_rec", "mass_dom")]))))
cat("Verified", nrow(r), "replicates;", sum(r$K), "true effects;",
    sum(a$roundoff_clamped_replicates), "roundoff-corrected replicates.\n")
b <- subset(m, scenario == "Additive only")
cat("Additive-only shares by coding (ranges):\n")
print(aggregate(estimated_share ~ coding, b, range))
f <- read.csv(file.path(d, "additive_only_false_call_rates.csv"))
cat("Maximum fraction with a false nonadditive call at PIP 0.9 or higher:\n")
print(sapply(f[c("any_false_rec", "any_false_dom")], max))
cat("Additive-only included replicates:", sum(r$scenario == "Additive only"), "\n")
analysis <- list(replicates = r, calibration = cal, discovery = calls, mixture = m,
                 confusion = read.csv(file.path(d, "coding_confusion.csv")),
                 classification = read.csv(file.path(d, "coding_classification.csv")))
cm_render_plots(analysis, config)
stopifnot(length(list.files(d, "[.]pdf$")) == 8L, length(list.files(d, "[.]png$")) == 8L)
cat("Verified and redrawn all eight PDF/PNG figures from the final count tables.\n")
