options(susie_mix.lead_pip_plots.run = FALSE)
source('script/sim/plot_lead_pip_coding.R')
stage <- file.path(project_dir, 'tmp/manuscript_lead_pip_update')
figdir <- file.path(stage, 'simulation_figures')
tabdir <- file.path(stage, 'simulation_tables')
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
a <- readRDS('simulation results/lead_pip_figures/lead_pip_analysis.rds')
r <- a$replicates
high <- r[!is.na(r$called_coding) & r$lead_pip >= .9, ]
summarize_high <- function(d, keys) {
  out <- aggregate(data.frame(n_leads = rep(1L, nrow(d)), n_exact = as.integer(d$outcome == 'exact'),
    n_wrong_coding = as.integer(d$outcome == 'wrong_coding'), n_wrong_snp = as.integer(d$outcome == 'wrong_snp'),
    sum_pip = d$lead_pip), d[keys], sum)
  out$mean_pip <- out$sum_pip / out$n_leads
  out$exact_accuracy <- out$n_exact / out$n_leads
  out
}
write.csv(summarize_high(high, 'called_coding'), file.path(tabdir, 'lead_pip_high_confidence.csv'), row.names = FALSE)
write.csv(summarize_high(high, c('pve', 'called_coding')), file.path(tabdir, 'lead_pip_high_confidence_by_pve.csv'), row.names = FALSE)
write.csv(summarize_high(high, c('pve', 'K', 'called_coding')), file.path(tabdir, 'lead_pip_high_confidence_by_pve_K.csv'), row.names = FALSE)
write.csv(summarize_high(high[high$scenario == 'Additive only', ], c('pve', 'K', 'called_coding')),
  file.path(tabdir, 'lead_pip_additive_only_high_confidence.csv'), row.names = FALSE)
stopifnot(nrow(r) == 87988, nrow(high) == 24438, sum(r$scenario == 'Additive only') == 7996,
  sum(high$scenario == 'Additive only' & high$called_coding == 'dominant') == 27,
  sum(high$scenario == 'Additive only' & high$called_coding == 'recessive') == 3)
mix <- read.csv('simulation results/coding_figures/mixture_recovery.csv')
mix <- mix[mix$scenario == 'Additive only', ]
write.csv(mix, file.path(tabdir, 'additive_only_pip_allocation.csv'), row.names = FALSE)
cal <- a$accuracy[a$accuracy$scope == 'All simulations', ]
cal <- aggregate(cal[c('n_leads', 'sum_lead_pip', 'n_exact')], cal[c('pve', 'called_coding', 'bin')], sum)
cal$mean_pip <- cal$sum_lead_pip / cal$n_leads
cal$accuracy <- cal$n_exact / cal$n_leads
write.csv(cal, file.path(tabdir, 'lead_pip_accuracy_pooled_K.csv'), row.names = FALSE)
for (name in list.files('simulation results/lead_pip_figures', '[.]csv$')) {
  dest <- if (name %in% c('file_audit.csv', 'excluded_replicates.csv')) paste0('lead_pip_', name) else name
  file.copy(file.path('simulation results/lead_pip_figures', name), file.path(tabdir, dest), overwrite = TRUE)
}

start <- function() par(mfrow = c(2,2), mar = c(3.3,3.9,2.1,.7), oma = c(2.2,.1,.4,.1),
  mgp = c(2,.65,0), tcl = -.25, family = 'sans', cex = .86)
panel <- function(pve, xlab, ylab, xlim = c(0,1), at = c(0,.5,1)) {
  plot(NA, xlim = xlim, ylim = c(0,1), axes = FALSE, xlab = xlab, ylab = '',
    main = paste0('PVE = ',pve*100,'%'), xaxs = 'i', yaxs = 'i')
  axis(1, at = at); axis(2, at = c(0,.25,.5,.75,1), las = 1); box(bty = 'l')
  mtext(ylab, 2, line = 2.9)
  abline(h = c(0,.25,.5,.75,1), col = '#EEEEEE')
}
save_plot <- function(stem, draw) {
  pdf(file.path(figdir,paste0(stem,'.pdf')), width = 7, height = 6.5, useDingbats = FALSE)
  tryCatch(draw(), finally = dev.off())
  png(file.path(stage,paste0(stem,'.png')), width = 1400, height = 1300, res = 200)
  tryCatch(draw(), finally = dev.off())
}
save_plot('lead_pip_accuracy_pooled_K', function() {
  start()
  for (pve in c(.1,.2,.3,.4)) {
    panel(pve,'Mean lead PIP in bin','Exact-match fraction')
    abline(0,1,lty=2,col='#777777')
    for (cl in cm_classes) {
      d <- cal[cal$pve == pve & cal$called_coding == cl, ]; d <- d[order(d$bin), ]
      lines(d$mean_pip,d$accuracy,col=lp_coding_colors[cl],lwd=1.2)
      points(d$mean_pip,d$accuracy,pch=16,col=lp_coding_colors[cl],cex=.55+.4*pmin(1,log10(d$n_leads)/3))
    }
  }
  lp_legend(lp_coding_labels,lp_coding_colors)
})
save_plot('additive_only_lead_coding_compact', function() {
  start()
  for (pve in c(.1,.2,.3,.4)) {
    panel(pve,'Number of causal SNPs','Share of unique leads',c(.7,5.3),1:5)
    d <- a$additive[a$additive$pve == pve, ]; d <- d[order(d$K), ]
    for(cl in cm_classes) lines(d$K,d[[paste0('share_unique_',cl)]],type='b',pch=16,col=lp_coding_colors[cl])
  }
  lp_legend(lp_coding_labels,lp_coding_colors)
})
save_plot('additive_only_pip_allocation_compact', function() {
  start()
  for (pve in c(.1,.2,.3,.4)) {
    panel(pve,'Number of causal SNPs','Share of total PIP mass',c(.7,5.3),1:5)
    for(j in seq_along(cm_classes)) {
      d <- mix[mix$pve == pve & mix$coding == cm_classes[j], ]; d <- d[order(d$K), ]
      x <- d$K + (j-2)*.06
      segments(x,d$lower,x,d$upper,col=lp_coding_colors[j])
      lines(x,d$estimated_share,type='b',pch=16,col=lp_coding_colors[j])
    }
  }
  lp_legend(lp_coding_labels,lp_coding_colors)
})
cat('Prepared three compact vector figures and checked numerical evidence tables.\n')
