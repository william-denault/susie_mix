# Redraw only coverage from the user's existing summary table.
# Do not read checkpoints or rerun the full plotting analysis.
code <- parse("script/sim/plot_simulations.R")
e <- new.env()
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("dir.create"))) break
  eval(z, e)
}
e$summary_table <- read.csv(file.path(e$output_dir, "metric_summary.csv"))
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("<-")) && is.symbol(z[[2]]) &&
      as.character(z[[2]]) %in% c("draw_figure", "save_figure")) eval(z, e)
}
e$save_figure("coverage", e$pure_rows, "coverage_pure")
e$save_figure("coverage", e$mixed_rows, "coverage_mixed")
cat("Refreshed the two coverage PDFs and PNGs from metric_summary.csv.\n")

# Descriptive assessment of the already-saved ROC counts at fixed FPRs.
roc <- readRDS(file.path(e$output_dir, "roc_counts.rds"))
for (cutoff in c(.001, .01)) {
  d <- roc[roc$fpr <= cutoff, ]
  best <- aggregate(tpr ~ scenario + pve + K + method, d, max)
  paired <- merge(best[best$method == "SuSiE", ], best[best$method == "SuSiE-mix", ],
                  by = c("scenario", "pve", "K"), suffixes = c("_add", "_mix"))
  paired$gain <- paired$tpr_mix - paired$tpr_add
  cat("\nFPR <=", cutoff, "(best available saved threshold):\n")
  print(aggregate(gain ~ scenario, paired, function(x) c(mean = mean(x), min = min(x), max = max(x))))
}
