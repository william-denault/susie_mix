# Redraw only the combined-L ROC figures from existing count summaries.
code <- parse("script/sim/plot_simulations.R")
e <- new.env()
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("dir.create"))) break
  eval(z, e)
}
e$roc_table <- readRDS(file.path(e$output_dir, "roc_counts.rds"))
for (z in code) {
  if (is.call(z) && identical(z[[1]], as.name("<-")) && is.symbol(z[[2]]) &&
      as.character(z[[2]]) %in% c("draw_figure", "save_figure")) eval(z, e)
}
e$save_figure("roc", e$pure_rows, "roc_pure_all_L")
e$save_figure("roc", e$mixed_rows, "roc_mixed_all_L")
cat("Saved the pure and mixed all-L ROC figures as PDF and PNG.\n")
