# Reproduce the dominant overview layout for the saved recessive cohort.
# Run with the project root as the first argument. Uses only base R.
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
plot_dir <- file.path(root, 'plot', 'one_cs_recessive')
all_cases <- read.csv(file.path(plot_dir, 'fit_mix_one_cs_recessive_plot_summary.csv'))
stopifnot(!anyDuplicated(all_cases[c('gene', 'tissue')]),
          all(all_cases$status == 'plotted'),
          all(all_cases$n_cs == 1),
          all(all_cases$lead_coding == 'recessive'))
parse_position <- function(ids) {
  stopifnot(all(grepl('^chr[^_]+_[0-9]+_[^_]+_[^_]+_b38(_|$)', ids)))
  pieces <- strsplit(ids, '_', fixed = TRUE)
  data.frame(chromosome = vapply(pieces, `[`, character(1), 1),
             position_bp = as.numeric(vapply(pieces, `[`, character(1), 2)))
}
add <- parse_position(all_cases$additive_lead_snp)
mix <- parse_position(all_cases$lead_snp)
stopifnot(all(add$chromosome == mix$chromosome))
all_cases$distance_bp <- abs(mix$position_bp - add$position_bp)
all_cases$distance_kb <- all_cases$distance_bp / 1000
all_cases$changed_lead <- all_cases$additive_lead_snp != all_cases$lead_snp
stopifnot(all(all_cases$changed_lead == !all_cases$same_lead_snp))
d <- all_cases[all_cases$changed_lead, ]
d <- d[order(-d$distance_bp, d$gene, d$tissue), ]
stopifnot(nrow(d) >= 10, all(d$distance_kb > 0),
          all(is.finite(d$lead_pip)), all(d$lead_pip >= 0 & d$lead_pip <= 1))

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
     ylab = 'Mixed lead PIP (recessive coding)')
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
  gene = c('ZDHHC2', 'ZNF354B', 'CFHR1', 'RASGRP3', 'NEK3', 'GUSB'),
  tissue = c('Blood', 'Skin', 'Spleen', 'Pancreas', 'Lung', 'Thyroid'),
  x = c(2.65, 2.23, 1.65, .28, 1.80, 2.76),
  y = c(.35, .82, .49, .83, .18, .12))
for (i in seq_len(nrow(labels))) {
  r <- d[d$gene == labels$gene[i] & d$tissue == labels$tissue[i], ]
  stopifnot(nrow(r) == 1)
  points(log10(r$distance_kb), r$lead_pip, pch = 21, bg = '#13786f',
         col = 'white', cex = 1.05)
  segments(log10(r$distance_kb), r$lead_pip, labels$x[i], labels$y[i],
           col = '#13786f', lwd = .7)
  text(labels$x[i], labels$y[i], labels$gene[i], col = '#07534d', cex = .82,
       pos = if (labels$gene[i] %in% c('NEK3', 'CFHR1')) 1 else 3, offset = .18)
}
legend_labels <- sprintf('One CS in both models (%d)', sum(g))
legend_pch <- 16
legend_col <- '#8798a4'
if (any(!g)) {
  legend_labels <- c(legend_labels, sprintf('Multiple additive CSs (%d)', sum(!g)))
  legend_pch <- c(legend_pch, 4)
  legend_col <- c(legend_col, '#b36a35')
}
legend('topleft', legend = legend_labels,
       pch = legend_pch, col = legend_col, bty = 'n', cex = .72)
mtext('Recessive one-CS examples: lead shifts and posterior support',
      outer = TRUE, side = 3, line = 2.4, cex = 1.45, font = 2)
mtext(sprintf('%s examples   |   %s changed leads   |   %s shifts > 100 kb',
              format(nrow(all_cases), big.mark = ','),
              format(nrow(d), big.mark = ','),
              format(sum(d$distance_bp > 100000), big.mark = ',')),
      outer = TRUE, side = 3, line = .6, cex = 1.1, col = '#444444')
mtext('Distance uses base-pair coordinates. PIP is coding-specific.',
      outer = TRUE, side = 1, line = 1.5, cex = .87, col = '#555555')
mtext('A large lead shift alone does not establish different credible sets or low LD.',
      outer = TRUE, side = 1, line = 3, cex = .87, col = '#555555')
invisible(dev.off())
cat(sprintf('Reviewed %d rows; %d changed leads; %d above 100 kb.\n',
            nrow(all_cases), nrow(d), sum(d$distance_bp > 100000)))
print(d[1:10, c('gene', 'tissue', 'distance_bp', 'lead_pip')], row.names = FALSE)
