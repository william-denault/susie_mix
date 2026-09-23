root <- normalizePath('.', winslash = '/')
out <- file.path(root, 'tmp/plot_pipeline_validation')
overrides <- list(project_dir = root,
  chunk_dir = file.path(root, 'tmp/slide_simulation_validation/smoke_chunks'),
  output_dir = out, write_png = TRUE)
run_pipeline <- function(refresh = FALSE) {
  e <- new.env()
  for (expr in parse('script/sim/plot_simulations.R')) {
    if (is.call(expr) && identical(expr[[1]], as.name('<-')) && is.symbol(expr[[2]])) {
      name <- as.character(expr[[2]])
      if (name %in% names(overrides)) expr[[3]] <- overrides[[name]]
      if (name == 'reuse_saved_summaries') expr[[3]] <- refresh
    }
    if (is.call(expr) && identical(expr[[1]], as.name('for')) &&
        identical(expr[[2]], as.name('metric'))) next
    if (is.call(expr) && identical(expr[[1]], as.name('save_roc_figures'))) next
    if (is.call(expr) && identical(expr[[1]], as.name('save_calibration_figures'))) next
    eval(expr, e)
  }
  e
}
e <- run_pipeline()
stopifnot(sum(e$audit$included) == 27,
  sum(e$calibration_table$n_causal) == sum(e$replicates$n_causal),
  sum(e$pooled_calibration_table$n_causal) == sum(e$calibration_table$n_causal))
baseline <- e$roc_table[e$roc_table$threshold == 0, ]
stopifnot(sum(e$calibration_table$n_snps) == sum(baseline$tp + baseline$fp))
cache <- file.path(out, 'pip_calibration_seed_counts.rds')
before <- tools::md5sum(cache)
stale <- e$summary_table
stale$coverage_lo <- -99
write.csv(stale, file.path(out, 'metric_summary.csv'), row.names = FALSE)
refresh <- run_pipeline(TRUE)
stopifnot(identical(before, tools::md5sum(cache)),
          isTRUE(all.equal(e$calibration_table, refresh$calibration_table)),
          isTRUE(all.equal(e$summary_table, refresh$summary_table)))
# Older compact caches trigger a complete rebuild, never reuse old bootstrap CIs.
old <- e$replicates
old$cs_size_sum_sq <- NULL
saveRDS(old, file.path(out, 'replicate_metrics.rds'))
upgraded <- run_pipeline(TRUE)
stopifnot('cs_size_sum_sq' %in% names(upgraded$replicates),
          isTRUE(all.equal(e$summary_table, upgraded$summary_table)))
e$save_figure('pip_calibration', e$scenario_groups$pure, 'pip_calibration_pure_all_L')
e$save_figure('coverage', e$scenario_groups$pure, 'coverage_normal_ci')
e$save_figure('cs_size', e$scenario_groups$pure, 'cs_size_gaussian_ci')
cat('PASS: checkpoint reader, analytic summary refresh, old-cache upgrade, exact calibration/ROC totals, and rendering.\n')
