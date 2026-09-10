args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
plot_dir <- file.path(root, 'plot', 'one_cs_dominant')
d <- read.csv(file.path(plot_dir, 'lead_snp_distance_ranking.csv'))
png(file.path(plot_dir, 'lead_snp_distance_overview.png'),
    width = 2700, height = 1350, res = 180)
par(mfrow = c(1, 2), oma = c(5.0, 0, 5, 0), family = 'sans',
    col.axis = '#444444', col.lab = '#333333', cex = 0.94)
par(mar = c(4, 11.5, 3, 1.5))
top <- d[10:1, ]
y <- barplot(top$distance_kb, horiz = TRUE, col = '#346b94', border = NA,
             xlim = c(0, 550), axes = FALSE,
             names.arg = paste(top$gene, top$tissue, sep = ' / '), las = 1,
             cex.names = .9, xlab = 'Lead-SNP distance (kb, GRCh38)')
axis(1, at = seq(0, 500, by = 100))
text(top$distance_kb + 8, y, sprintf('%.1f', top$distance_kb), adj = 0, cex = .85)
title(main = 'Ten largest physical shifts', adj = 0, cex.main = 1.15)
par(mar = c(4, 5, 3, 2))
plot(NA, xlim = c(-3.2, 3.2), ylim = c(-.02, 1.12), axes = FALSE,
     xlab = 'Lead-SNP distance (kb, logarithmic scale)',
     ylab = 'Mixed lead PIP (dominant coding)')
axis(1, at = -3:3, labels = c('0.001', '0.01', '0.1', '1', '10', '100', '1,000'))
axis(2, at = seq(0, 1, .25), las = 1)
abline(h = .5, v = 2, lty = 2, col = '#bbbbbb', lwd = .8)
g <- d$additive_n_cs == 1
points(log10(d$distance_kb[g]), d$lead_pip[g], pch = 16, cex = .58,
       col = adjustcolor('#8798a4', alpha.f = .6))
points(log10(d$distance_kb[!g]), d$lead_pip[!g], pch = 4, cex = .9,
       col = '#b36a35')
title(main = 'Distance and support among changed leads', adj = 0, cex.main = 1.15)
labels <- data.frame(
  gene = c('EIF3C', 'TSLP', 'CMTM6', 'TPSD1', 'CPXM1', 'MAN2C1'),
  tissue = c('Blood Vessel', 'Esophagus', 'Thyroid', 'Small Intestine', 'Spleen', 'Uterus'),
  x = c(2.32, 2.63, 1.43, -.12, .75, 2.50),
  y = c(.23, .55, .64, 1.07, .92, .74))
for (i in seq_len(nrow(labels))) {
  r <- d[d$gene == labels$gene[i] & d$tissue == labels$tissue[i], ]
  points(log10(r$distance_kb), r$lead_pip, pch = 21, bg = '#13786f',
         col = 'white', cex = 1.05)
  segments(log10(r$distance_kb), r$lead_pip, labels$x[i], labels$y[i],
           col = '#13786f', lwd = .7)
  text(labels$x[i], labels$y[i], labels$gene[i], col = '#07534d', cex = .82,
       pos = if (labels$gene[i] == 'EIF3C') 1 else 3, offset = .18)
}
legend('topleft', legend = c('One CS in both models (450)', 'Two additive CSs (16)'),
       pch = c(16, 4), col = c('#8798a4', '#b36a35'), bty = 'n', cex = .72)
mtext('Dominant one-CS examples: lead shifts and posterior support',
      outer = TRUE, side = 3, line = 2.4, cex = 1.45, font = 2)
mtext('1,340 examples   |   466 changed leads   |   46 shifts > 100 kb',
      outer = TRUE, side = 3, line = .6, cex = 1.1, col = '#444444')
mtext('Distance uses base-pair coordinates. PIP is coding-specific.',
      outer = TRUE, side = 1, line = 1.5, cex = .87, col = '#555555')
mtext('A large lead shift alone does not establish different credible sets or low LD.',
      outer = TRUE, side = 1, line = 3, cex = .87, col = '#555555')
invisible(dev.off())
