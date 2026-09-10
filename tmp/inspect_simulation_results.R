.libPaths(c('C:/Users/willi/AppData/Local/R/win-library/4.6',
            'C:/Users/willi/AppData/Local/R/win-library/4.5', .libPaths()))
cat('Packages:', sapply(c('ggplot2', 'data.table'), requireNamespace, quietly = TRUE), '\n')
files <- c(
  'additive_add1_rec0_dom0_n500_L10_pve0.2_seed1e+06_reps200_chunk1.RData',
  'additive_add1_rec0_dom0_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData',
  'additive_recessive_dominant_add1_rec1_dom1_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData')
for (f in files) {
  e <- new.env()
  load(file.path('simulation results/chunks', f), envir = e)
  cat('\nFILE:', f, '\nobjects:', ls(e), '\nreplicates:', length(e$results), '\n')
  good <- which(vapply(e$results, function(x) is.null(x$error), logical(1)))
  cat('Failures:', length(e$results) - length(good), '\n')
  if (length(good)) {
    x <- e$results[[good[1]]]
    str(x, max.level = 1)
    print(x$settings)
    print(x$metrics)
    str(x$susie_cs)
    str(x$susie_mix_cs)
    cat('PIP lengths:', length(x$susie_pip), length(x$susie_mix_pip_snp), '\n')
  }
}
