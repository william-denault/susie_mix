# Redraw only the combined-L ROC figures from existing count summaries.
code <- parse("script/sim/plot_simulations.R")
e <- new.env()
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("dir.create"))) break
  eval(z, e)
}
e$roc_table <- e$add_fdr(readRDS(file.path(e$output_dir, "roc_counts.rds")))
e$pooled_roc_table <- e$pool_curve_counts(e$roc_table)
saveRDS(e$pooled_roc_table, file.path(e$output_dir, "roc_counts_all_L.rds"))
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("<-")) && is.symbol(z[[2]]) &&
      as.character(z[[2]]) %in% c("draw_figure", "save_figure")) eval(z, e)
}
e$save_figure("roc", e$pure_rows, "roc_pure_all_L")
e$save_figure("roc", e$mixed_rows, "roc_mixed_all_L")
e$save_figure("power_fdr", e$pure_rows, "power_fdr_pure_all_L")
e$save_figure("power_fdr", e$mixed_rows, "power_fdr_mixed_all_L")
cat("Saved pooled ROC and power-FDR figures as PDF and PNG.\n")
