# Overview of the completed four-panel plots. Only base R is needed, so the
# same helper can also redraw an overview from a saved plot-summary CSV.
plot_one_cs_lead_distance_overview <- function(
    plot_summary,
    expected_coding,
    output_file,
    top_n = 20L) {

  expected_coding <- match.arg(expected_coding, c("dominant", "recessive"))
  if (length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1L || top_n != as.integer(top_n)) {
    stop("top_n must be a positive integer.")
  }
  plot_summary <- as.data.frame(plot_summary)
  if (!"status" %in% names(plot_summary)) {
    stop("The plot summary has no status column.")
  }
  completed <- plot_summary[which(plot_summary$status == "plotted"), , drop = FALSE]
  n_completed <- nrow(completed)
  required <- c("gene", "tissue", "additive_lead_snp", "lead_snp",
                "additive_n_cs", "n_cs", "lead_coding", "lead_pip")
  if (n_completed > 0L) {
    missing <- setdiff(required, names(completed))
    if (length(missing) > 0L) {
      stop("The plot summary is missing: ", paste(missing, collapse = ", "))
    }
    if (anyDuplicated(completed[c("gene", "tissue")])) {
      stop("The plot summary contains duplicate plotted gene-tissue pairs.")
    }
    if (anyNA(completed[c("n_cs", "lead_coding")]) ||
        any(completed$n_cs != 1L | completed$lead_coding != expected_coding)) {
      stop("Expected one mixed CS with a ", expected_coding, " lead in every plotted case.")
    }
  } else {
    # The runner returns only audit columns when no cases succeed. Still
    # render an empty overview, rather than leaving a previous run's figure.
    completed <- data.frame(
      gene = character(), tissue = character(), additive_lead_snp = character(),
      lead_snp = character(), additive_n_cs = integer(), n_cs = integer(),
      lead_coding = character(), lead_pip = numeric()
    )
  }

  parse_snp <- function(ids) {
    pattern <- "^(chr[^_]+)_([0-9]+)_[^_]+_[^_]+_b38(?:_|$)"
    valid <- !is.na(ids) & grepl(pattern, ids, perl = TRUE)
    chromosome <- rep(NA_character_, length(ids))
    position <- rep(NA_real_, length(ids))
    chromosome[valid] <- sub("^(chr[^_]+)_.*$", "\\1", ids[valid])
    position[valid] <- as.numeric(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", ids[valid]))
    data.frame(chromosome = chromosome, position = position)
  }
  add <- parse_snp(completed$additive_lead_snp)
  mix <- parse_snp(completed$lead_snp)
  valid_coordinates <- !is.na(add$position) & !is.na(mix$position) &
    add$position > 0 & mix$position > 0 &
    !is.na(add$chromosome) & !is.na(mix$chromosome) &
    add$chromosome == mix$chromosome
  n_invalid <- sum(!valid_coordinates)
  if (n_invalid > 0L) {
    warning(n_invalid, " plotted lead pair(s) lack comparable GRCh38 coordinates; ",
            "excluded from the distance panels.")
  }
  completed$distance_bp <- abs(mix$position - add$position)
  completed$distance_bp[!valid_coordinates] <- NA_real_
  completed$distance_kb <- completed$distance_bp / 1000
  changed <- !is.na(completed$additive_lead_snp) & !is.na(completed$lead_snp) &
    completed$additive_lead_snp != completed$lead_snp
  d <- completed[which(changed & valid_coordinates), , drop = FALSE]
  d <- d[order(-d$distance_bp, d$gene, d$tissue), , drop = FALSE]
  top <- head(d, top_n)
  # Different alleles may occupy the same base-pair position. Keep them in
  # the ranking but do not place distance zero on the logarithmic x-axis.
  scatter_ok <- d$distance_bp > 0 & is.finite(d$lead_pip) &
    d$lead_pip >= 0 & d$lead_pip <= 1
  scatter <- d[which(scatter_ok), , drop = FALSE]

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(
    output_file, width = 2700, height = max(1350, 900 + 45 * min(top_n, 40L)), res = 180,
    type = if (capabilities("cairo")) "cairo" else getOption("bitmapType")
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(1, 2), oma = c(5, 0, 5, 0), family = "sans",
      col.axis = "#444444", col.lab = "#333333", cex = .94)
  par(mar = c(4, 11.5, 3, 1.5))
  if (nrow(top) > 0L) {
    bars <- top[rev(seq_len(nrow(top))), , drop = FALSE]
    x_max <- max(550, max(bars$distance_kb) * 1.15)
    y <- barplot(bars$distance_kb, horiz = TRUE, col = "#346b94", border = NA,
                 xlim = c(0, x_max), axes = FALSE,
                 names.arg = paste(bars$gene, bars$tissue, sep = " / "), las = 1,
                 cex.names = .9, xlab = "Lead-SNP distance (kb, GRCh38)")
    axis(1)
    text(bars$distance_kb + x_max * .015, y,
         sprintf("%.1f", bars$distance_kb), adj = 0, cex = .85)
    title(main = sprintf("%d largest physical shifts", nrow(top)), adj = 0, cex.main = 1.15)
  } else {
    plot.new()
    title(main = "Largest physical shifts", adj = 0, cex.main = 1.15)
    text(.5, .5, if (n_completed == 0L) "No completed four-panel plots" else
      "No changed leads with comparable coordinates", cex = .85)
  }

  par(mar = c(4, 5, 3, 2))
  log_limits <- c(-3.2, max(3.2, if (nrow(scatter) > 0L)
    max(log10(scatter$distance_kb)) + .2 else 3.2))
  plot(NA, xlim = log_limits, ylim = c(-.02, 1.12), axes = FALSE,
       xlab = "Lead-SNP distance (kb, logarithmic scale)",
       ylab = paste0("Mixed lead PIP (", expected_coding, " coding)"))
  ticks <- seq(-3, floor(log_limits[2]))
  axis(1, at = ticks, labels = format(10^ticks, scientific = FALSE, trim = TRUE,
                                     big.mark = ",", drop0trailing = TRUE))
  axis(2, at = seq(0, 1, .25), las = 1)
  abline(h = .5, v = 2, lty = 2, col = "#bbbbbb", lwd = .8)
  title(main = "Distance and support among changed leads", adj = 0, cex.main = 1.15)
  # Preserve the distinction between a one-CS comparison and a selected
  # additive lead from multiple CSs (or the no-CS fallback).
  groups <- list(
    list(keep = scatter$additive_n_cs == 1L, label = "One CS in both models", pch = 16, col = "#8798a4"),
    list(keep = scatter$additive_n_cs > 1L, label = "Multiple additive CSs", pch = 4, col = "#b36a35"),
    list(keep = scatter$additive_n_cs == 0L, label = "No additive CS", pch = 2, col = "#78579a"),
    list(keep = is.na(scatter$additive_n_cs), label = "Additive CS count unavailable", pch = 3, col = "#666666")
  )
  legend_labels <- character()
  legend_pch <- integer()
  legend_col <- character()
  for (group in groups) {
    rows <- which(group$keep)
    if (length(rows) == 0L) next
    points(log10(scatter$distance_kb[rows]), scatter$lead_pip[rows],
           pch = group$pch, cex = if (group$pch == 16) .58 else .9,
           col = adjustcolor(group$col, alpha.f = .65))
    legend_labels <- c(legend_labels, sprintf("%s (%d)", group$label, length(rows)))
    legend_pch <- c(legend_pch, group$pch)
    legend_col <- c(legend_col, group$col)
  }

  # Keep the examples annotated in the reviewed figures. Missing examples
  # are skipped, so a smaller or partially successful run is also supported.
  annotations <- if (expected_coding == "dominant") {
    data.frame(
      gene = c("EIF3C", "TSLP", "CMTM6", "TPSD1", "CPXM1", "MAN2C1"),
      tissue = c("Blood Vessel", "Esophagus", "Thyroid", "Small Intestine", "Spleen", "Uterus"),
      x = c(2.32, 2.63, 1.43, -.12, .75, 2.50),
      y = c(.23, .55, .64, 1.07, .92, .74), pos = c(1, 3, 3, 3, 3, 3)
    )
  } else {
    data.frame(
      gene = c("ZDHHC2", "ZNF354B", "CFHR1", "RASGRP3", "NEK3", "GUSB"),
      tissue = c("Blood", "Skin", "Spleen", "Pancreas", "Lung", "Thyroid"),
      x = c(2.65, 2.23, 1.65, .28, 1.80, 2.76),
      y = c(.35, .82, .49, .83, .18, .12), pos = c(3, 3, 1, 3, 1, 3)
    )
  }
  for (i in seq_len(nrow(annotations))) {
    r <- scatter[scatter$gene == annotations$gene[i] &
                   scatter$tissue == annotations$tissue[i], , drop = FALSE]
    if (nrow(r) != 1L) next
    points(log10(r$distance_kb), r$lead_pip, pch = 21, bg = "#13786f",
           col = "white", cex = 1.05)
    segments(log10(r$distance_kb), r$lead_pip, annotations$x[i], annotations$y[i],
             col = "#13786f", lwd = .7)
    text(annotations$x[i], annotations$y[i], annotations$gene[i],
         col = "#07534d", cex = .82, pos = annotations$pos[i], offset = .18)
  }
  if (length(legend_labels) > 0L) {
    legend("topleft", legend = legend_labels, pch = legend_pch, col = legend_col,
           bty = "n", cex = .72)
  } else {
    text(mean(log_limits), .5, "No changed leads with positive distance and valid PIP")
  }
  coding_title <- paste0(toupper(substr(expected_coding, 1, 1)), substring(expected_coding, 2))
  mtext(paste0(coding_title, " one-CS examples: lead shifts and posterior support"),
        outer = TRUE, side = 3, line = 2.4, cex = 1.45, font = 2)
  count_text <- function(x) format(x, big.mark = ",", trim = TRUE)
  mtext(sprintf("%s examples   |   %s changed leads   |   %s shifts > 100 kb",
                count_text(n_completed), count_text(sum(changed)),
                count_text(sum(d$distance_bp > 100000))),
        outer = TRUE, side = 3, line = .6, cex = 1.1, col = "#444444")
  omissions <- character()
  if (n_invalid > 0L) omissions <- c(omissions, paste(n_invalid, "pairs with unavailable distance"))
  if (any(!scatter_ok)) omissions <- c(omissions, paste(sum(!scatter_ok), "pairs omitted from log/PIP panel"))
  if (length(omissions) > 0L) {
    mtext(paste(omissions, collapse = "; "), outer = TRUE, side = 1, line = .1,
          cex = .8, col = "#555555")
  }
  mtext("Distance uses base-pair coordinates. PIP is coding-specific.",
        outer = TRUE, side = 1, line = 1.5, cex = .87, col = "#555555")
  mtext("A large lead shift alone does not establish different credible sets or low LD.",
        outer = TRUE, side = 1, line = 3, cex = .87, col = "#555555")
  message("Lead-SNP distance overview: ", output_file)
  invisible(list(distance_data = d, top_cases = top, n_completed = n_completed,
                 n_changed = sum(changed), n_invalid = n_invalid,
                 n_scatter = nrow(scatter), output_file = output_file))
}
