e <- new.env()
for (expr in parse('script/sim/plot_simulations.R')) {
  if (is.call(expr) && identical(expr[[1]], as.name('dir.create'))) break
  eval(expr, e)
}
for (expr in parse('script/sim/plot_simulations.R')) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      as.character(expr[[2]]) %in% c('draw_figure', 'save_figure', 'save_roc_figures')) eval(expr, e)
}
old_output <- e$output_dir
e$summary_table <- read.csv(file.path(old_output, 'metric_summary.csv'))
e$roc_table <- readRDS(file.path(old_output, 'roc_counts.rds'))
e$pooled_roc_table <- readRDS(file.path(old_output, 'roc_counts_all_L.rds'))
e$output_dir <- file.path(e$project_dir, 'tmp/plot_axis_preview')
dir.create(e$output_dir, recursive = TRUE, showWarnings = FALSE)
for (metric in c('coverage', 'purity', 'cs_size', 'power', 'power_fdr'))
  e$save_figure(metric, e$scenario_groups$pure, paste0(metric, '_pure'))
cat('Axis previews saved in', e$output_dir, '\n')
