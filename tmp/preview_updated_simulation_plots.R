e <- new.env()
for (expr in parse('script/sim/plot_simulations.R')) {
  if (is.call(expr) && identical(expr[[1]], as.name('dir.create'))) break
  eval(expr, e)
}
for (expr in parse('script/sim/plot_simulations.R')) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      as.character(expr[[2]]) %in% c('figure_data', 'draw_figure', 'save_figure', 'save_roc_figures')) eval(expr, e)
}
old_output <- e$output_dir
e$summary_table <- read.csv(file.path(old_output, 'metric_summary.csv'))
e$roc_table <- readRDS(file.path(old_output, 'roc_counts.rds'))
e$pooled_roc_table <- readRDS(file.path(old_output, 'roc_counts_all_L.rds'))
e$output_dir <- file.path(e$project_dir, 'tmp/plot_shared_axis_preview')
dir.create(e$output_dir, recursive = TRUE, showWarnings = FALSE)
# Observe the limits passed to every data panel on both rendering devices.
e$plot <- function(x, ..., ylim) {
  e$observed_limits <- rbind(e$observed_limits, ylim)
  graphics::plot(x, ..., ylim = ylim)
}
for (metric in c('coverage', 'purity', 'cs_size', 'power', 'power_fdr')) {
  e$observed_limits <- NULL
  scenarios <- e$scenario_groups$pairs_2
  e$save_figure(metric, scenarios, paste0(metric, '_pairs_2'))
  expected <- e$figure_y_limits(metric, e$figure_data(metric, scenarios), c(0, e$fdr_max))
  stopifnot(nrow(e$observed_limits) == 50,
            all(e$observed_limits[, 1] == expected[1]),
            all(e$observed_limits[, 2] == expected[2]))
  cat(metric, ': every panel uses', expected, '\n')
}
cat('Axis previews saved in', e$output_dir, '\n')
