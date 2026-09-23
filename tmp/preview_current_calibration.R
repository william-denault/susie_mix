e <- new.env()
for (expr in parse('script/sim/plot_simulations.R')) {
  if (is.call(expr) && identical(expr[[1]], as.name('dir.create'))) break
  eval(expr, e)
}
for (expr in parse('script/sim/plot_simulations.R')) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      as.character(expr[[2]]) %in% c('figure_data', 'draw_figure', 'save_figure', 'save_roc_figures', 'save_calibration_figures')) eval(expr, e)
}
e$summary_table <- read.csv(file.path(e$output_dir, 'metric_summary.csv'))
e$roc_table <- readRDS(file.path(e$output_dir, 'roc_counts.rds'))
e$pooled_roc_table <- readRDS(file.path(e$output_dir, 'roc_counts_all_L.rds'))
e$calibration_table <- read.csv(file.path(e$output_dir, 'pip_calibration.csv'))
e$pooled_calibration_table <- read.csv(file.path(e$output_dir, 'pip_calibration_all_L.csv'))
baseline <- e$roc_table[e$roc_table$threshold == 0, ]
bins <- aggregate(e$calibration_table[c('n_snps', 'n_causal')],
  e$calibration_table[c('scenario', 'pve', 'K', 'method')], sum)
key <- function(d) paste(d$scenario, d$pve, d$K, d$method)
rows <- match(key(baseline), key(bins))
stopifnot(!anyNA(rows), all(baseline$tp == bins$n_causal[rows]),
  all(baseline$tp + baseline$fp == bins$n_snps[rows]),
  sum(e$calibration_table$n_causal) == sum(e$pooled_calibration_table$n_causal),
  sum(e$calibration_table$n_snps) == sum(e$pooled_calibration_table$n_snps))
e$output_dir <- file.path(e$project_dir, 'tmp/calibration_by_method_preview')
dir.create(e$output_dir, recursive = TRUE, showWarnings = FALSE)
obsolete <- file.path(e$output_dir, c('pip_calibration_pure_L1.pdf',
  'pip_calibration_pure_all_L.png', 'pip_calibration_pairs_1_by_L.pdf'))
retained <- file.path(e$output_dir, c('pip_calibration_all_L.csv', 'coverage_pure.pdf'))
invisible(file.create(c(obsolete, retained)))
e$save_calibration_figures()
stopifnot(!any(file.exists(obsolete)), all(file.exists(retained)),
  length(list.files(e$output_dir, pattern = '^pip_calibration_.*_all_L_susie.*\\.pdf$')) == 15,
  length(list.files(e$output_dir, pattern = '^pip_calibration_.*_all_L_susie.*\\.png$')) == 15)
for (method in e$method_names) {
  d <- e$figure_data('pip_calibration', e$scenario_groups$pure, methods = method)
  stopifnot(identical(unique(d$method), method), all(d$K == 0))
}
cat('PASS: calibration/ROC counts agree; 15 pooled figures per format; one method per figure; obsolete plots removed safely.\n')
