# Redraw ROC figures from the existing count summaries, without reading chunks.
code <- parse("script/sim/plot_simulations.R")
e <- new.env()
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("dir.create"))) break
  eval(z, e)
}
e$roc_table <- e$add_fdr(readRDS(file.path(e$output_dir, "roc_counts.rds")))
e$pooled_roc_table <- e$pool_curve_counts(e$roc_table)
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("<-")) && is.symbol(z[[2]]) &&
      as.character(z[[2]]) %in% c("draw_figure", "save_figure", "save_roc_figures")) eval(z, e)
}
e$save_roc_figures()
cat("Saved per-L ROC figures and the enabled overview/collection outputs.\n")
