# Base-R plotting and auditable summaries for the tentative four-panel figure.
# Source em_cs_agreement_utils.R first: reuse its cohort validation and exact-first,
# one-to-one distance matching without changing the existing figure or its settings.
em_tentative_codings <- c("additive", "recessive", "dominant", "unknown")
em_tentative_colors <- c(additive = "#477A9B", recessive = "#CD7139",
                         dominant = "#197F73", unknown = "#888888")

em_tentative_coding <- function(x) {
  ifelse(x %in% em_tentative_codings, as.character(x), "unknown")
}

em_tentative_featured <- function(tissues) {
  eligible <- tissues[which(tissues$ranking_eligible), , drop = FALSE]
  low <- head(eligible[order(eligible$agreement_rate, -eligible$n_regions, eligible$tissue), ], 2L)
  high <- eligible[!eligible$tissue %in% low$tissue, , drop = FALSE]
  high <- head(high[order(-high$agreement_rate, -high$n_regions, high$tissue), ], 2L)
  low$selection <- rep("Most disagreement", nrow(low))
  high$selection <- rep("Most agreement", nrow(high))
  rbind(low, high)
}

summarize_em_cs_tentative <- function(res, cs, em_priors, association_threshold = 1e-8,
                                     minimum_mean_reads = 100, minimum_tissue_regions = 100L) {
  cs <- as.data.frame(cs, stringsAsFactors = FALSE)
  em_agreement_require(cs, c("cs_snps", "cs_size_snps"), "EM CS summary for panels C and D")
  base <- summarize_em_cs_agreement(res, cs, association_threshold, minimum_mean_reads,
                                    minimum_tissue_regions, em_priors = em_priors)
  regions <- base$region_audit
  tissues <- base$tissue_ranking
  # Descriptive Wilson intervals; repeated genes across tissues are not independent replications.
  z <- qnorm(.975); n <- tissues$n_regions; p <- tissues$agreement_rate
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  tissues$agreement_low95 <- ifelse(n > 0, centre - half, NA_real_)
  tissues$agreement_high95 <- ifelse(n > 0, centre + half, NA_real_)
  featured <- em_tentative_featured(tissues)
  tables <- lapply(featured$tissue, function(tissue) {
    d <- regions[regions$count_eligible & regions$tissue == tissue, ]
    levels <- sort(unique(c(0L, d$n_add, d$n_em)))
    tab <- as.data.frame(table(additive_cs_count = factor(d$n_add, levels),
                              em_cs_count = factor(d$n_em, levels)))
    data.frame(tissue = tissue, additive_cs_count = as.integer(as.character(tab$additive_cs_count)),
               em_cs_count = as.integer(as.character(tab$em_cs_count)), n_regions = tab$Freq)
  })
  tables <- if (length(tables)) do.call(rbind, tables) else data.frame(tissue = character(),
    additive_cs_count = integer(), em_cs_count = integer(), n_regions = integer())

  cs$region_id <- em_agreement_key(cs$gene, cs$tissue)
  cs$cs_id <- paste(cs$region_id, cs$model_key, cs$cs_number, sep = "|")
  pairs <- base$matched_cs_pairs; extras <- base$unpaired_cs
  wanted <- unique(c(pairs$additive_cs_id, pairs$em_cs_id, extras$cs_id))
  cs <- cs[cs$cs_id %in% wanted, , drop = FALSE]
  if (anyDuplicated(cs$cs_id)) stop("Duplicate CS IDs in selected summaries.")
  members <- lapply(cs$cs_snps, function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(character())
    unique(trimws(strsplit(s, ";", fixed = TRUE)[[1]]))
  })
  sizes <- suppressWarnings(as.numeric(cs$cs_size_snps))
  cs$membership_status <- vapply(seq_len(nrow(cs)), function(i) {
    m <- members[[i]]
    if (!length(m) || any(!nzchar(m))) return("missing_members")
    if (!all(em_agreement_coordinates(m)$valid)) return("invalid_member_id")
    if (!is.finite(sizes[i]) || length(m) != sizes[i]) return("size_mismatch")
    if (!cs$lead_snp[i] %in% m) return("lead_absent_from_members")
    "ok"
  }, character(1))
  cs$size_snps <- sizes
  cs$size_snps[cs$membership_status != "ok"] <- NA_real_
  a <- match(pairs$additive_cs_id, cs$cs_id); m <- match(pairs$em_cs_id, cs$cs_id)
  pairs$additive_cs_number <- cs$cs_number[a]; pairs$em_cs_number <- cs$cs_number[m]
  pairs$additive_size_snps <- cs$size_snps[a]; pairs$em_size_snps <- cs$size_snps[m]
  pairs$additive_lead_pip <- cs$lead_pip[a]
  pairs$em_lead_coding <- em_tentative_coding(pairs$em_lead_coding)
  pairs$membership_status <- ifelse(cs$membership_status[a] == "ok" & cs$membership_status[m] == "ok",
                                      "ok", "unavailable_membership")
  pairs$intersection_snps <- pairs$union_snps <- pairs$jaccard <- rep(NA_real_, nrow(pairs))
  for (i in which(pairs$membership_status == "ok")) {
    pairs$intersection_snps[i] <- length(intersect(members[[a[i]]], members[[m[i]]]))
    pairs$union_snps[i] <- length(union(members[[a[i]]], members[[m[i]]]))
    pairs$jaccard[i] <- pairs$intersection_snps[i] / pairs$union_snps[i]
  }
  pairs$size_ratio_em_to_additive <- pairs$em_size_snps / pairs$additive_size_snps
  pairs$log2_size_ratio <- log2(pairs$size_ratio_em_to_additive)
  e <- match(extras$cs_id, cs$cs_id)
  extras$size_snps <- cs$size_snps[e]
  extras$membership_status <- cs$membership_status[e]
  extras$lead_coding <- em_tentative_coding(extras$lead_coding)
  resolution <- do.call(rbind, lapply(em_tentative_codings, function(coding) {
    d <- pairs[pairs$em_lead_coding == coding & is.finite(pairs$log2_size_ratio), ]
    data.frame(coding = coding, n_pairs = nrow(d),
      median_size_ratio = if (nrow(d)) median(d$size_ratio_em_to_additive) else NA_real_,
      n_smaller_em = sum(d$size_ratio_em_to_additive < 1),
      n_equal = sum(d$size_ratio_em_to_additive == 1),
      n_larger_em = sum(d$size_ratio_em_to_additive > 1))
  }))
  extra_summary <- do.call(rbind, lapply(em_tentative_codings, function(coding) {
    d <- extras[extras$model == "EM SuSiE-mix" & extras$lead_coding == coding, ]
    valid <- is.finite(d$lead_pip) & d$lead_pip >= 0 & d$lead_pip <= 1
    data.frame(coding = coding, n_cs = nrow(d), n_valid_pip = sum(valid),
      n_pip_ge_half = sum(d$lead_pip[valid] >= .5), n_pip_ge_point9 = sum(d$lead_pip[valid] >= .9),
      n_no_susie_cs = sum(!d$other_model_has_cs),
      n_plotted = sum(valid & is.finite(d$size_snps) & d$size_snps > 0))
  }))
  selection <- as.data.frame(table(tissue = regions$tissue, status = regions$selection_status))
  names(selection)[3] <- "n_regions"
  selection <- selection[selection$n_regions > 0, ]
  list(schema = "em_cs_tentative_v1", iteration = base$iteration, settings = base$settings,
       n_regions = sum(regions$count_eligible), n_tissues = sum(tissues$n_regions > 0),
       tissue_ranking = tissues, selected_tissues = featured, tissue_cs_count_tables = tables,
       matched_cs_pairs = pairs, unpaired_cs = extras, resolution_summary = resolution,
       extra_em_support_summary = extra_summary, lead_distance_summary = base$lead_distance_summary,
       selection_summary = selection,
       region_audit = regions[c("gene", "tissue", "em_iteration", "n_add", "n_em", "selection_status",
                                "count_eligible", "matching_status", "n_pairs", "n_extra_add", "n_extra_em")],
       cs_membership_audit = cs[c("gene", "tissue", "model_key", "cs_number", "cs_id", "membership_status")])
}

em_tentative_distance_axis <- function(xmax) {
  ticks <- c(0, 10^(0:ceiling(xmax)))
  labels <- c("0", "1 bp", "10 bp", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb",
              "10 Mb", "100 Mb", "1 Gb")
  keep <- log10(1 + ticks) <= xmax
  axis(1, at = log10(1 + ticks[keep]), labels = labels[seq_along(ticks)][keep], cex.axis = .76)
}

plot_em_cs_tentative <- function(result, output_file = NULL, preview_label = NULL) {
  if (is.null(preview_label)) preview_label <- result$preview_label
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    if (tolower(tools::file_ext(output_file)) == "pdf") {
      grDevices::pdf(output_file, width = 18, height = 13, useDingbats = FALSE)
    } else grDevices::png(output_file, width = 3600, height = 2600, res = 200,
                          type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE, after = FALSE)
  layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 5), nrow = 2, byrow = TRUE), heights = c(1, 1))
  par(oma = c(4.3, .6, 4, .6), family = "sans", cex = 1,
      col.axis = "#444444", col.lab = "#333333")
  pairs <- result$matched_cs_pairs; extras <- result$unpaired_cs
  fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

  # A: four exact count tables, two agreement extremes in each direction.
  par(mar = c(1, 1.5, 3, 1))
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  title("A   Where do CS counts agree?", adj = 0, cex.main = 1.18)
  text(.5, .99, "Rows: SuSiE CSs  |  Columns: EM CSs  |  Cells: gene-tissue regions", cex = .8)
  if (!nrow(result$selected_tissues)) text(.5, .5, "No tissues meet the ranking threshold.")
  for (i in seq_len(nrow(result$selected_tissues))) {
    f <- result$selected_tissues[i, ]; tab <- result$tissue_cs_count_tables
    tab <- tab[tab$tissue == f$tissue, ]
    levels <- sort(unique(tab$additive_cs_count)); k <- length(levels)
    left <- if (i %% 2L) .055 else .56
    top <- if (i <= 2L) .78 else .31
    w <- .39 / k; h <- .28 / k
    color <- if (f$selection == "Most agreement") "#197F73" else "#B86338"
    text(left + .195, top + .14, f$tissue, font = 2, cex = 1, col = color)
    text(left + .195, top + .102, sprintf("%s | n = %s", f$selection, fmt(f$n_regions)), cex = .76)
    text(left + .195, top + .066, sprintf("Same: %s (%.1f%%) | Different: %s",
      fmt(f$same_cs_count), 100 * f$agreement_rate, fmt(f$different_cs_count)), cex = .75)
    for (a in seq_len(k)) for (m in seq_len(k)) {
      n <- tab$n_regions[tab$additive_cs_count == levels[a] & tab$em_cs_count == levels[m]]
      x <- left + (m - 1L) * w; y <- top - a * h
      fill <- grDevices::adjustcolor(if (a == m) "#197F73" else "#CD7139",
                                     alpha.f = if (n) .09 + .55 * sqrt(n / max(tab$n_regions)) else .025)
      rect(x, y, x + w, y + h, col = fill, border = "#e0e4e7", lwd = .5)
      text(x + w / 2, y + h / 2, if (levels[a] == 0 && levels[m] == 0) "-" else n,
           cex = min(.85, 6.8 / k))
    }
    text(left - .012, top - (seq_len(k) - .5) * h, levels, adj = 1, cex = min(.85, 6.8 / k))
    text(left + (seq_len(k) - .5) * w, top + .018, levels, cex = min(.85, 6.8 / k))
  }

  # B: exact identities are counted separately, retaining zero distances for different alleles.
  par(mar = c(5.3, 4.8, 3, 1))
  distances <- list("Paired, changed lead" = pairs$distance_bp[!pairs$same_lead_snp],
    "Extra EM to nearest SuSiE" = extras$nearest_other_distance_bp[extras$model == "EM SuSiE-mix"],
    "Extra SuSiE to nearest EM" = extras$nearest_other_distance_bp[extras$model == "SuSiE"])
  vals <- unlist(distances); vals <- vals[is.finite(vals)]
  xmax <- max(6.05, if (length(vals)) log10(1 + max(vals)) else 0)
  plot(NA, xlim = c(0, xmax), ylim = c(0, 1.04), axes = FALSE,
       xlab = "", ylab = "Cumulative fraction of CS leads")
  em_tentative_distance_axis(xmax)
  axis(2, at = seq(0, 1, .25), labels = paste0(seq(0, 100, 25), "%"), las = 1)
  abline(h = seq(0, 1, .25), col = "#e8e8e8")
  abline(v = log10(1 + c(1e4, 1e5)), col = "#cccccc", lty = 3)
  cols <- c("#477A9B", "#CD7139", "#8A6CA4"); labels <- character()
  for (j in seq_along(distances)) {
    d <- sort(distances[[j]][is.finite(distances[[j]])])
    if (length(d)) {
      x <- sort(unique(d)); y <- findInterval(x, d) / length(d)
      lines(c(0, log10(1 + x), xmax), c(0, y, 1), type = "s", col = cols[j], lty = j, lwd = 2)
    }
    labels <- c(labels, sprintf("%s (n = %s)", names(distances)[j], fmt(length(d))))
  }
  legend("bottomright", labels, col = cols, lty = 1:3, lwd = 2, cex = .77,
         bty = "o", box.col = NA, bg = "white")
  title("B   Lead retention and displacement: all tissues", adj = 0, cex.main = 1.18)
  mtext("Lead-SNP distance (logarithmic spacing; zero retained)", 1, 2.5, cex = .9)
  mtext(sprintf("Not in curves: %s identical pairs; %s extras with no other-model CS",
    fmt(sum(pairs$same_lead_snp)), fmt(sum(!extras$other_model_has_cs))), 1, 4.2, cex = .75)
  if (!length(vals)) text(xmax / 2, .5, "No usable lead distances")

  # C: changed paired leads only, so a large mass of shared leads cannot hide differences.
  par(mar = c(5.3, 4.8, 3.4, 1))
  changed <- pairs[!pairs$same_lead_snp, ]
  d <- changed[is.finite(changed$jaccard), ]
  plot(NA, xlim = c(0, xmax), ylim = c(-.05, 1.08), axes = FALSE,
       xlab = "", ylab = "CS SNP overlap (Jaccard index)")
  em_tentative_distance_axis(xmax); axis(2, at = seq(0, 1, .25), las = 1)
  abline(h = c(0, .5, 1), col = "#e8e8e8")
  abline(v = log10(1 + c(1e4, 1e5)), col = "#cccccc", lty = 3)
  if (nrow(d)) points(log10(1 + d$distance_bp), d$jaccard, pch = 16, cex = .55,
                      col = grDevices::adjustcolor(em_tentative_colors[d$em_lead_coding], alpha.f = .4))
  else text(xmax / 2, .5, "No changed pairs with valid CS membership")
  title("C   Does a shifted lead change the credible set?", adj = 0, cex.main = 1.18)
  mtext(sprintf("Changed pairs: %s | identical sets: %s | disjoint sets: %s",
    fmt(nrow(d)), fmt(sum(d$jaccard == 1)), fmt(sum(d$jaccard == 0))), 3, .35, cex = .79)
  present <- em_tentative_codings[em_tentative_codings %in% d$em_lead_coding]
  if (length(present)) legend("bottomleft", paste0(present, " EM lead"), pch = 16,
    col = em_tentative_colors[present], bty = "o", box.col = NA, bg = "white", cex = .77)
  mtext("Distance between paired leads (physical distance does not measure LD)", 1, 2.5, cex = .86)
  mtext(sprintf("Jaccard: 0 = disjoint; 1 = identical SNP sets. Missing membership: %s pairs.",
    fmt(nrow(changed) - nrow(d))), 1, 4.2, cex = .75)

  # D has two complementary views: matched-set resolution and unmatched-EM support.
  par(mar = c(5.3, 4.8, 3.4, .6))
  d <- pairs[is.finite(pairs$log2_size_ratio), ]
  groups <- em_tentative_codings[em_tentative_codings %in% d$em_lead_coding]
  yrange <- if (nrow(d)) range(c(-1, 1, d$log2_size_ratio)) else c(-1, 1)
  if (length(groups)) {
    values <- lapply(groups, function(g) d$log2_size_ratio[d$em_lead_coding == g])
    boxplot(values, at = seq_along(groups), axes = FALSE, outline = FALSE,
            ylim = yrange, col = grDevices::adjustcolor(em_tentative_colors[groups], .35),
            border = em_tentative_colors[groups], ylab = "CS size ratio: EM / SuSiE")
    axis(1, at = seq_along(groups), labels = paste0(tools::toTitleCase(groups), "\nn = ", lengths(values)),
         cex.axis = .73, tick = FALSE)
    # Show the full range so omitted boxplot outlier symbols do not conceal extremes.
    for (j in seq_along(values)) {
      segments(j, min(values[[j]]), j, max(values[[j]]), col = em_tentative_colors[groups[j]], lty = 3)
      points(rep(j, 2), range(values[[j]]), pch = 3, cex = .55, col = em_tentative_colors[groups[j]])
    }
  } else {
    plot(NA, xlim = c(0, 1), ylim = yrange, axes = FALSE, xlab = "", ylab = "CS size ratio: EM / SuSiE")
    text(.5, 0, "No valid paired CS sizes", cex = .85)
  }
  ticks <- pretty(yrange, n = 5); ticks <- ticks[ticks >= yrange[1] & ticks <= yrange[2]]
  axis(2, at = ticks, labels = format(signif(2^ticks, 2), trim = TRUE), las = 1, cex.axis = .75)
  abline(h = 0, col = "#555555", lty = 2)
  title("D   Resolution and support", adj = 0, cex.main = 1.12)
  mtext("Matched CSs, grouped by EM lead coding", 3, .35, cex = .77)
  mtext("Below 1: fewer SNPs in the EM CS", 1, 3.4, cex = .78)
  mtext(sprintf("Missing size: %s pairs | whiskers: 1.5 IQR; +: range",
                fmt(nrow(pairs) - nrow(d))), 1, 4.4, cex = .61)

  par(mar = c(5.3, 4.1, 3.4, 1))
  all_extra <- extras[extras$model == "EM SuSiE-mix", ]
  d <- all_extra[is.finite(all_extra$size_snps) & all_extra$size_snps > 0 &
                  is.finite(all_extra$lead_pip) & all_extra$lead_pip >= 0 & all_extra$lead_pip <= 1, ]
  size_max <- max(2, if (nrow(d)) log10(max(d$size_snps)) else 0)
  plot(NA, xlim = c(0, size_max), ylim = c(0, 1.05), axes = FALSE,
       xlab = "", ylab = "EM lead PIP (coding-specific)")
  ticks <- 0:floor(size_max)
  axis(1, at = ticks, labels = format(10^ticks, big.mark = ",", scientific = FALSE, trim = TRUE), cex.axis = .8)
  axis(2, at = seq(0, 1, .25), las = 1)
  abline(h = c(.5, .9), lty = 3, col = "#cccccc")
  if (nrow(d)) points(log10(d$size_snps), d$lead_pip, pch = 16, cex = .65,
                      col = grDevices::adjustcolor(em_tentative_colors[d$lead_coding], .5))
  else text(size_max / 2, .5, "No extra EM CSs with size and PIP", cex = .8)
  title("Additional EM CSs", adj = 0, cex.main = 1.12)
  mtext(sprintf("n = %s | PIP >= 0.9: %s", fmt(nrow(d)), fmt(sum(d$lead_pip >= .9))), 3, .35, cex = .8)
  present <- em_tentative_codings[em_tentative_codings %in% d$lead_coding]
  if (length(present)) legend("topright", tools::toTitleCase(present), pch = 16,
    col = em_tentative_colors[present], cex = .7, bty = "o", box.col = NA, bg = "white")
  mtext("Unique SNPs per CS (logarithmic scale)", 1, 2.5, cex = .79)
  mtext(sprintf("Unavailable size/PIP: %s; no SuSiE CS: %s",
    fmt(nrow(all_extra) - nrow(d)), fmt(sum(!all_extra$other_model_has_cs))), 1, 4.2, cex = .65)

  heading <- sprintf("SuSiE versus EM SuSiE-mix | iteration %d", result$iteration)
  if (!is.null(preview_label)) heading <- paste(preview_label, heading, sep = " - ")
  mtext(heading, outer = TRUE, side = 3, line = 2.05, cex = 1.55, font = 2)
  mtext(sprintf("%s eligible gene-tissue regions across %s tissues | tentative four-panel comparison",
    fmt(result$n_regions), fmt(result$n_tissues)), outer = TRUE, side = 3, line = .55, cex = 1)
  mtext(sprintf("Cohort: P < %s; mean reads >= %s; at least one CS in either model. Tissue ranking requires n >= %s.",
    format(result$settings$association_threshold, scientific = TRUE), result$settings$minimum_mean_reads,
    result$settings$minimum_tissue_regions), outer = TRUE, side = 1, line = 2.55, cex = .8)
  mtext("Pairing reserves identical lead IDs, then minimizes total distance without reusing CSs. Additional = unmatched under this rule.",
        outer = TRUE, side = 1, line = 1.35, cex = .8)
  mtext("SNP sets collapse genetic codings. Pairing does not establish biological equivalence; smaller sets do not establish better coverage.",
        outer = TRUE, side = 1, line = .15, cex = .8)
  invisible(result)
}

write_em_cs_tentative <- function(result, output_dir, show_plot = FALSE, preview_label = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- c("tissue_ranking", "selected_tissues", "tissue_cs_count_tables", "matched_cs_pairs",
              "unpaired_cs", "resolution_summary", "extra_em_support_summary", "lead_distance_summary",
              "selection_summary", "region_audit", "cs_membership_audit")
  for (name in tables) if (!is.null(result[[name]]))
    write.csv(result[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
  write.csv(data.frame(em_iteration = result$iteration, as.data.frame(result$settings)),
            file.path(output_dir, "analysis_settings.csv"), row.names = FALSE)
  # The small replot bundle omits the large region audit, which remains in CSV.
  bundle <- result; bundle$region_audit <- NULL
  bundle$preview_label <- if (!is.null(preview_label)) preview_label else result$preview_label
  saveRDS(bundle, file.path(output_dir, "tentative_plot_data.rds"))
  if (is.null(preview_label)) preview_label <- result$preview_label
  for (extension in c("png", "pdf")) plot_em_cs_tentative(result,
    file.path(output_dir, paste0("em_cs_comparison_tentative_4panel.", extension)), preview_label)
  if (isTRUE(show_plot)) plot_em_cs_tentative(result, preview_label = preview_label)
  message("Tentative figure and audit tables: ", normalizePath(output_dir, winslash = "/"))
  invisible(result)
}
