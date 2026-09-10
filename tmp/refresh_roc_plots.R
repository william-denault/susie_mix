# Redraw ROC figures from the existing count summaries, without reading chunks.
code <- parse("script/sim/plot_simulations.R")
e <- new.env()
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("dir.create"))) break
  eval(z, e)
}
e$roc_table <- readRDS(file.path(e$output_dir, "roc_counts.rds"))
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("<-")) && is.symbol(z[[2]]) &&
      as.character(z[[2]]) %in% c("draw_figure", "save_figure", "save_roc_figures")) eval(z, e)
}
e$save_roc_figures()
cat("Saved nine individual ROC figures (PDF and PNG), one per L and scenario group.\n")
