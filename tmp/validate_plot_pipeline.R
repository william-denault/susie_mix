root <- normalizePath('.', winslash = '/')
out <- file.path(root, 'tmp/plot_pipeline_validation')
overrides <- list(project_dir = root,
  chunk_dir = file.path(root, 'tmp/slide_simulation_validation/smoke_chunks'),
  output_dir = out, bootstrap_reps = 50, write_png = TRUE)
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
refresh <- run_pipeline(TRUE)
stopifnot(identical(before, tools::md5sum(cache)),
          isTRUE(all.equal(e$calibration_table, refresh$calibration_table)))
e$save_figure('pip_calibration', e$scenario_groups$pure, 'pip_calibration_pure_all_L')
cat('PASS: normal reader, exact calibration/ROC totals, cache reuse, and calibration rendering.\n')
