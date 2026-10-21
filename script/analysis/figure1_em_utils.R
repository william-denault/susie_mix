# Figure 1 v2: four count tables, two ECDFs and posterior coding support.
# Source em_analysis_utils.R and em_cs_agreement_utils.R first. Base R only.
fig1_fmt <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
fig1_pct <- function(x) if (is.finite(x)) sprintf("%.1f%%", 100 * x) else "NA"
fig1_members <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  unique(trimws(strsplit(x, ";", fixed = TRUE)[[1L]]))
}

summarize_figure1_em <- function(res, cs, priors, settings) {
  z <- summarize_em_cs_agreement(res, cs, settings$association_threshold,
    settings$minimum_mean_reads, settings$minimum_tissue_regions, em_priors = priors)
  r <- z$region_audit[z$region_audit$count_eligible, , drop = FALSE]
  if (!nrow(r)) stop("No regions pass the shared Figure 1 cohort filters.")
  z$schema <- "figure1_em_tables_ecdf_v2"
  z$settings <- settings
  z$cohort <- r
  p <- z$matched_cs_pairs
  z$metrics <- data.frame(
    n_regions = nrow(r), n_genes = length(unique(r$gene)), n_tissues = length(unique(r$tissue)),
    n_add_cs = sum(r$n_add), n_em_cs = sum(r$n_em),
    n_same_count = sum(r$n_add == r$n_em), n_different_count = sum(r$n_add != r$n_em),
    count_agreement_rate = mean(r$n_add == r$n_em),
    n_pairs = nrow(p), n_same_lead = sum(p$same_lead_snp), n_changed_lead = sum(!p$same_lead_snp),
    lead_retention_rate = if (nrow(p)) mean(p$same_lead_snp) else NA_real_,
    n_matching_unavailable_regions = sum(r$matching_status != "matched"),
    n_model_only_extra_cs = sum(!z$unpaired_cs$other_model_has_cs))
  # The selected four tissues affect A only; the headline, B, C and D use all tissues.
  z$tissue_count_table <- as.data.frame(table(tissue = r$tissue, n_add = r$n_add, n_em = r$n_em))
  names(z$tissue_count_table)[4L] <- "n_regions"
  z$tissue_count_table <- z$tissue_count_table[z$tissue_count_table$n_regions > 0, ]
  # Do not retain a prior-based panel under the new posterior panel's name.
  z$tissue_coding_priors <- NULL
  z$priors_used <- priors
  z
}

# Read each selected gene once. Sum the FULL component posterior over each
# coding, not only predictors inside the 95% CS, and not SNP PIPs.
fig1_collect_coding_support <- function(context, cohort, cs) {
  source(file.path(context$project_dir, "script/analysis/generate_summary_results.R"), local = TRUE)
  cs <- as.data.frame(cs, stringsAsFactors = FALSE)
  em_agreement_require(cs, c("component", "cs_snps", "cs_size_snps"), "Posterior coding CS summary")
  cs$region_id <- em_agreement_key(cs$gene, cs$tissue)
  cs <- cs[cs$model_key == "weighted_fit_mix" & cs$region_id %in% cohort$region_id, , drop = FALSE]
  if (anyDuplicated(cs[c("region_id", "cs_number")])) stop("Duplicate EM CS summary rows.")
  cs_parts <- split(cs, cs$region_id)
  by_gene <- split(cohort, cohort$gene)
  classes <- c("additive", "recessive", "dominant")
  rows <- list()
  audit <- cohort[c("region_id", "gene", "tissue", "n_em")]
  audit$n_components_used <- 0L
  audit$status <- ifelse(audit$n_em > 0, "pending", "no_em_cs")
  message("Reading posterior coding support for ", length(by_gene), " cohort genes...")
  for (g in seq_along(by_gene)) {
    r <- by_gene[[g]]; gene <- names(by_gene)[g]
    if (!gene %in% context$manifest$gene) stop("Cohort gene absent from the EM manifest: ", gene)
    file <- file.path(context$iteration_dir, "results", paste0(gene, ".rds"))
    fits <- readRDS(file)
    if (!is.list(fits) || !is.null(fits[["error"]]) ||
        !isTRUE(attr(fits, "em_iteration", exact = TRUE) == context$iteration))
      stop("Invalid selected-iteration EM gene result: ", gene)
    if (is.null(names(fits)) || anyNA(names(fits)) || anyDuplicated(names(fits)))
      stop("Invalid EM tissue names for ", gene)
    for (i in seq_len(nrow(r))) {
      label <- paste(gene, r$tissue[i], sep = " / ")
      x <- fits[[r$tissue[i]]]
      if (!is.list(x)) stop("Missing EM tissue fit: ", label)
      fit <- x$weighted_fit_mix
      if (!is.matrix(fit$alpha) || !isTRUE(fit$converged)) stop("Missing or unconverged EM fit: ", label)
      map <- as.data.frame(x$mix_predictor_map)
      validate_predictor_map(fit, map, label)
      if (anyNA(map$coding) || any(!map$coding %in% classes) || anyDuplicated(map$predictor_name))
        stop("Invalid EM predictor codings/names: ", label)
      expected_prior <- unlist(context$priors[match(r$tissue[i], context$priors$tissue),
        c("pi_add", "pi_rec", "pi_dom")], use.names = FALSE)
      observed_prior <- x$mix_coding_prior[classes]
      if (length(observed_prior) != 3L || anyNA(expected_prior) ||
          !isTRUE(all.equal(unname(observed_prior), expected_prior, tolerance = 1e-10)))
        stop("Saved EM prior does not match selected iteration: ", label)
      sets <- get_cs(fit)
      part <- cs_parts[[r$region_id[i]]]
      if (is.null(part)) part <- cs[FALSE, , drop = FALSE]
      if (length(sets) != r$n_em[i] || nrow(part) != r$n_em[i] ||
          !setequal(part$cs_number, seq_along(sets)))
        stop("Saved EM CS count and summary disagree: ", label, ". Regenerate EM summaries.")
      if (!length(sets)) next
      components <- fit$sets$cs_index
      if (length(components) != length(sets) || anyNA(components) ||
          any(components < 1 | components > nrow(fit$alpha) | components != floor(components)) ||
          anyDuplicated(components)) stop("Invalid CS-to-component mapping: ", label)
      alpha <- fit$alpha[components, , drop = FALSE]
      if (any(!is.finite(alpha)) || any(alpha < 0 | alpha > 1) ||
          any(abs(rowSums(alpha) - 1) > 1e-6))
        stop("Reported-CS component posteriors must be finite probabilities summing to one: ", label)
      for (j in seq_along(sets)) {
        s <- part[part$cs_number == j, , drop = FALSE]
        indices <- sets[[j]]
        if (!is.numeric(indices) || !length(indices) || anyNA(indices) ||
            any(indices < 1 | indices > nrow(map) | indices != floor(indices)) || anyDuplicated(indices))
          stop("Invalid CS predictor indices: ", label)
        lead <- get_cs_lead_predictor(fit, sets, j)
        ids <- unique(as.character(map$snp[indices]))
        if (!isTRUE(s$component == components[j]) ||
            !isTRUE(s$lead_snp == map$snp[lead]) || !isTRUE(s$lead_coding == map$coding[lead]) ||
            !isTRUE(s$cs_size_snps == length(ids)) || !setequal(fig1_members(s$cs_snps), ids) ||
            length(fit$pip) != nrow(map) || !is.finite(fit$pip[lead]) ||
            !isTRUE(all.equal(as.numeric(s$lead_pip), as.numeric(fit$pip[lead]), tolerance = 1e-8)))
          stop("Saved EM CS and summary disagree: ", label, " / CS ", j, ". Regenerate EM summaries.")
        support <- vapply(classes, function(coding) sum(alpha[j, map$coding == coding]), numeric(1))
        rows[[length(rows) + 1L]] <- data.frame(
          region_id = r$region_id[i], gene = gene, tissue = r$tissue[i], em_iteration = context$iteration,
          cs_number = j, component = components[j], lead_snp = s$lead_snp, lead_coding = s$lead_coding,
          lead_pip = s$lead_pip, cs_size_snps = length(ids), alpha_sum = sum(alpha[j, ]),
          posterior_additive = support[1], posterior_recessive = support[2], posterior_dominant = support[3],
          posterior_mass_inside_cs = sum(alpha[j, indices]), row.names = NULL)
      }
      k <- match(r$region_id[i], audit$region_id)
      audit$n_components_used[k] <- length(sets); audit$status[k] <- "included"
    }
    if (g %% 100L == 0L || g == length(by_gene)) message("  Read ", g, " / ", length(by_gene), " genes.")
  }
  support <- if (length(rows)) do.call(rbind, rows) else data.frame(
    region_id = character(), gene = character(), tissue = character(), em_iteration = integer(),
    cs_number = integer(), component = integer(), lead_snp = character(), lead_coding = character(),
    lead_pip = numeric(), cs_size_snps = integer(), alpha_sum = numeric(),
    posterior_additive = numeric(), posterior_recessive = numeric(), posterior_dominant = numeric(),
    posterior_mass_inside_cs = numeric())
  stopifnot(nrow(support) == sum(cohort$n_em), all(audit$n_components_used == audit$n_em))
  list(cs = support, regions = audit)
}

fig1_add_coding_summary <- function(z, support) {
  z$coding_cs_support <- support$cs
  z$coding_region_audit <- support$regions
  columns <- paste0("posterior_", c("additive", "recessive", "dominant"))
  aggregate_one <- function(tissue, r, d) {
    means <- if (nrow(d)) colMeans(d[columns]) else setNames(rep(NA_real_, 3L), columns)
    data.frame(tissue = tissue, n_cohort_regions = nrow(r), n_regions_with_em_cs = sum(r$n_em > 0),
      n_em_cs = nrow(d), mean_additive = unname(means[1]), mean_recessive = unname(means[2]),
      mean_dominant = unname(means[3]), mean_nonadditive = unname(means[2] + means[3]))
  }
  z$coding_tissue_support <- do.call(rbind, lapply(sort(unique(z$cohort$tissue)), function(tissue)
    aggregate_one(tissue, z$cohort[z$cohort$tissue == tissue, ],
      support$cs[support$cs$tissue == tissue, , drop = FALSE])))
  z$coding_overall_support <- aggregate_one("All tissues", z$cohort, support$cs)
  z
}

fig1_distance_axis <- function(xmax) {
  ticks <- c(0, 10^(0:8))
  labels <- c("0", "1 bp", "10 bp", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  keep <- log10(1 + ticks) <= xmax
  axis(1, at = log10(1 + ticks[keep]), labels = labels[keep], cex.axis = .75)
}

# Coordinates of a right-continuous ECDF, including its jump at zero.
fig1_ecdf_steps <- function(values, xmax) {
  d <- sort(values[is.finite(values)])
  if (!length(d)) return(data.frame(x = numeric(), y = numeric()))
  x <- sort(unique(d))
  data.frame(x = c(0, log10(1 + x), xmax), y = c(0, findInterval(x, d) / length(d), 1))
}

plot_figure1_em <- function(z, file = NULL) {
  if (!identical(z$schema, "figure1_em_tables_ecdf_v2") || is.null(z$coding_tissue_support))
    stop("This plot requires the v2 tables/ECDF/posterior-support bundle; rerun on RCC.")
  if (!is.null(file)) {
    if (tolower(tools::file_ext(file)) == "pdf") pdf(file, width = 16, height = 12.5, useDingbats = FALSE)
    else png(file, width = 4000, height = 3125, res = 250,
      type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    on.exit(dev.off(), add = TRUE)
  }
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE, after = FALSE)
  scale <- if (is.null(file)) min(1, dev.size("in") / c(16, 12.5)) else 1
  par(family = "sans", cex = scale, mex = scale, mgp = c(2.1, .7, 0),
      col.axis = "#444444", col.lab = "#333333", las = 1, xaxs = "i", yaxs = "i")
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4))
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(.5, .975, paste0("SuSiE versus EM SuSiE-mix | iteration ", z$iteration), font = 2, cex = 1.5)
  m <- z$metrics
  text(.5, .949, sprintf("P < %s | mean reads >= %s | >=1 CS in either model | %s regions across %s tissues",
    format(z$settings$association_threshold, scientific = TRUE), z$settings$minimum_mean_reads,
    fig1_fmt(m$n_regions), m$n_tissues), cex = .94, col = "#555555")
  titles <- c("A   CS-count agreement: four high-disagreement tissues",
              "B   How far do leads differ? All tissues",
              "C   Distance to TSS for the extra leads in B",
              "D   Posterior coding support in the same cohort")
  for (i in 1:4) text(if (i %% 2L) .018 else .533, if (i <= 2L) .911 else .466,
    titles[i], adj = 0, font = 2, cex = 1.03)
  text(.022, .879, sprintf("%s same CS count overall  |  %s / %s regions",
    fig1_pct(m$count_agreement_rate), fig1_fmt(m$n_same_count), fig1_fmt(m$n_regions)),
    adj = 0, font = 2, cex = 1.03, col = "#287D72")
  text(.553, .879, sprintf("%s of matched leads retained  |  %s / %s pairs",
    fig1_pct(m$lead_retention_rate), fig1_fmt(m$n_same_lead), fig1_fmt(m$n_pairs)),
    adj = 0, font = 2, cex = 1.03, col = "#346B94")
  text(.5, .039, sprintf("Same parent cohort in A-D. Table ranking: CS-count disagreement rate, minimum n = %s; B-D use all cohort tissues.",
    z$settings$minimum_tissue_regions), cex = .76, col = "#555555")
  text(.5, .024, sprintf("Distance matching is descriptive, not evidence of shared signals or LD. Matching unavailable for %s selected regions.",
    fig1_fmt(m$n_matching_unavailable_regions)), cex = .76, col = "#555555")
  if (!is.null(z$preview_label)) text(.5, .008, z$preview_label, font = 2, col = "#B32929", cex = .77)

  # A: literal count tables from the established comparison, with a pooled headline.
  par(fig = c(.012, .495, .492, .856), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(.5, .981, "Rows: SuSiE CSs  |  Columns: EM CSs  |  Cells: region counts", cex = .75)
  if (!nrow(z$selected_tissues)) text(.5, .5, "No tissue meets the table-ranking threshold.")
  for (i in seq_len(nrow(z$selected_tissues))) {
    f <- z$selected_tissues[i, ]; tab <- z$tissue_cs_count_tables[z$tissue_cs_count_tables$tissue == f$tissue, ]
    levels <- sort(unique(tab$additive_cs_count)); k <- length(levels)
    left <- if (i %% 2L) .055 else .56; top <- if (i <= 2L) .79 else .33
    w <- .39 / k; h <- .285 / k
    text(left + .195, top + .13, f$tissue, font = 2, cex = .9, col = "#B75C32")
    text(left + .195, top + .09, sprintf("Disagreement rank %s | n = %s", f$disagreement_rank, fig1_fmt(f$n_regions)), cex = .67)
    text(left + .195, top + .052, sprintf("Same: %s (%s) | Different: %s", f$same_cs_count,
      fig1_pct(f$agreement_rate), f$different_cs_count), cex = .67)
    for (a in seq_len(k)) for (b in seq_len(k)) {
      n <- tab$n_regions[tab$additive_cs_count == levels[a] & tab$em_cs_count == levels[b]]
      x <- left + (b - 1L) * w; y <- top - a * h
      fill <- adjustcolor(if (a == b) "#287D72" else "#C16A40",
        alpha.f = if (n) .10 + .58 * sqrt(n / max(tab$n_regions)) else .025)
      rect(x, y, x + w, y + h, col = fill, border = "#E0E4E7", lwd = .5)
      text(x + w / 2, y + h / 2, if (levels[a] == 0 && levels[b] == 0) "-" else n, cex = min(.8, 6.6 / k))
    }
    text(left - .015, top - (seq_len(k) - .5) * h, levels, adj = 1, cex = min(.8, 6.6 / k))
    text(left + (seq_len(k) - .5) * w, top + .015, levels, cex = min(.8, 6.6 / k))
  }

  # B and C retain the populations and ECDF definitions in the preferred figure.
  colors <- c("#346B94", "#B75C32", "#7B68A5")
  ecdf_panel <- function(values, fig, ylab, xlabel, legends, cols, ltys) {
    par(fig = fig, mar = c(3.3, 4.5, .6, 1), new = TRUE)
    finite <- unlist(values); finite <- finite[is.finite(finite)]
    xmax <- max(log10(1 + 1e6), if (length(finite)) log10(1 + max(finite)) else 0)
    plot(NA, xlim = c(0, xmax), ylim = c(0, 1.035), axes = FALSE,
      xlab = xlabel, ylab = "", cex.lab = .85)
    mtext(ylab, side = 2, line = 3.2, las = 0, cex = .85 * scale)
    fig1_distance_axis(xmax)
    axis(2, at = seq(0, 1, .25), labels = paste0(seq(0, 100, 25), "%"), cex.axis = .82)
    abline(h = seq(0, 1, .25), col = "#EEEEEE")
    abline(v = log10(1 + c(1e4, 1e5)), lty = 3, col = "#CCCCCC")
    for (j in seq_along(values)) {
      steps <- fig1_ecdf_steps(values[[j]], xmax)
      if (nrow(steps)) lines(steps$x, steps$y, type = "s", col = cols[j], lty = ltys[j], lwd = 2)
    }
    legend("bottomright", legends, col = cols, lty = ltys, lwd = 2, bty = "o", box.col = NA,
      bg = "white", cex = .72)
    if (!length(finite)) text(xmax / 2, .55, "No eligible distances")
  }
  labels <- sprintf("%s (n = %s)", names(z$distances),
    vapply(z$distances, function(d) fig1_fmt(sum(is.finite(d))), character(1)))
  ecdf_panel(z$distances, c(.51, .985, .545, .852), "Cumulative fraction of CS leads",
    "Lead-SNP distance (logarithmic spacing; zero retained)", labels, colors, 1:3)
  par(fig = c(.55, .985, .492, .537), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0, .77, sprintf("Curves exclude %s identical pairs and %s extra CSs in model-only regions.",
    fig1_fmt(m$n_same_lead), fig1_fmt(m$n_model_only_extra_cs)), adj = 0, cex = .72)
  text(0, .25, "Changed pairs use one-to-one matches; extra CSs use the nearest other-model lead.", adj = 0, cex = .7, col = "#555555")

  labels <- sprintf("%s (TSS n = %s / %s in B)", names(z$tss_distances),
    vapply(z$tss_distances, function(d) fig1_fmt(sum(is.finite(d))), character(1)),
    vapply(z$tss_distances, function(d) fig1_fmt(length(d)), character(1)))
  ecdf_panel(z$tss_distances, c(.01, .495, .105, .436), "Cumulative fraction of extra CS leads",
    "Absolute lead-SNP distance to gene TSS (logarithmic spacing)", labels, colors[2:3], 2:3)
  par(fig = c(.045, .495, .058, .097), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0, .78, sprintf("Same extra CSs as B; %s missing TSS distances excluded.",
    fig1_fmt(sum(z$tss_distance_summary$n_missing_tss))), adj = 0, cex = .73)
  text(0, .24, "Model-only regions excluded from both ECDFs. TSS proximity does not establish accuracy.", adj = 0, cex = .69, col = "#555555")

  # D: equal weight per reported-CS component, pooled over the shared cohort.
  ov <- z$coding_overall_support
  par(fig = c(.555, .985, .407, .444), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0, .8, sprintf("All %s EM CSs: additive %s | recessive %s | dominant %s",
    fig1_fmt(ov$n_em_cs), fig1_pct(ov$mean_additive), fig1_pct(ov$mean_recessive), fig1_pct(ov$mean_dominant)),
    adj = 0, cex = .78)
  text(0, .22, "Mean component posterior probabilities; additive mass retained.", adj = 0, cex = .74, col = "#555555")
  par(fig = c(.515, .985, .105, .405), mar = c(3.3, 8.6, .4, 2.9), new = TRUE)
  d <- z$coding_tissue_support
  d <- d[order(-d$mean_nonadditive, d$tissue, na.last = TRUE), ]
  n <- nrow(d); yy <- rev(seq_len(n))
  xmax <- max(c(.01, d$mean_recessive, d$mean_dominant), na.rm = TRUE) * 1.08
  plot(NA, xlim = c(0, xmax), ylim = c(.2, n + 2.8), axes = FALSE,
    xlab = "Mean posterior coding support", ylab = "", cex.lab = .85)
  ticks <- pretty(c(0, xmax), n = 5); ticks <- ticks[ticks >= 0 & ticks <= xmax]
  axis(1, at = ticks, labels = paste0(format(100 * ticks, trim = TRUE), "%"), cex.axis = .8)
  # axis() silently skips crowded labels; explicitly draw every cohort tissue.
  text(0, yy, d$tissue, pos = 2, xpd = NA, cex = min(.74, 17.5 / max(1, n)))
  abline(v = ticks, col = "#EEEEEE"); abline(h = yy, col = "#F3F3F3")
  cc <- c("#9A6BB0", "#C18B35")
  for (j in 1:2) {
    v <- d[[c("mean_recessive", "mean_dominant")[j]]]
    segments(0, yy + c(.13, -.13)[j], v, yy + c(.13, -.13)[j], col = adjustcolor(cc[j], .45), lwd = 1.2)
    points(v, yy + c(.13, -.13)[j], pch = c(16, 15)[j], col = cc[j], cex = .67)
  }
  text(xmax, yy, fig1_fmt(d$n_em_cs), pos = 4, xpd = NA, cex = min(.7, 17 / max(1, n)))
  text(xmax, n + 1, "CSs", pos = 4, xpd = NA, cex = .72)
  if (any(d$n_em_cs == 0)) text(0, yy[d$n_em_cs == 0], "No EM CS", pos = 4, cex = .65, col = "#777777")
  legend("top", c("Recessive", "Dominant"), col = cc, pch = c(16, 15), ncol = 2, bty = "n", cex = .8)
  par(fig = c(.555, .985, .058, .097), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0, .78, "Each reported EM CS contributes its full component posterior; counts shown at right.", adj = 0, cex = .71)
  text(0, .24, "This is posterior support among reported CSs, not the fraction of nonadditive regions.", adj = 0, cex = .7, col = "#555555")
  invisible(z)
}

fig1_caption <- function(z) {
  m <- z$metrics; d <- z$coding_overall_support
  c("Figure 1. SuSiE versus EM SuSiE-mix.",
    sprintf("One parent cohort: %s gene-tissue regions across %s tissues, iteration %s. Minimum marginal additive P < %s, mean reads >= %s, >=1 CS in either model. These are gene-tissue regions, not independent genes.",
      m$n_regions, m$n_tissues, z$iteration, z$settings$association_threshold, z$settings$minimum_mean_reads),
    sprintf("A. Overall CS-count agreement: %s / %s (%s). Four tissues with highest disagreement rate among tissues with >=%s cohort regions; ties use larger region count, then tissue name. The ranking cutoff affects only table selection, not other panels or the headline. Rows are SuSiE counts; columns are EM counts; cells are exact region counts. The 0/0 cell is outside the cohort.",
      m$n_same_count, m$n_regions, fig1_pct(m$count_agreement_rate), z$settings$minimum_tissue_regions),
    sprintf("B. Exact biological lead retention: %s / %s one-to-one matched CS pairs (%s), ignoring coding. Reserve identical variant IDs first, then minimize total physical distance for remaining CSs. ECDFs display changed pairs and extra CSs with a counterpart-model CS. Extra-CS distances are to the nearest other-model lead and are not one-to-one matches. %s model-only extra CSs have no defined nearest lead and are omitted. %s selected regions have unavailable matching and are audited.",
      m$n_same_lead, m$n_pairs, fig1_pct(m$lead_retention_rate), m$n_model_only_extra_cs, m$n_matching_unavailable_regions),
    "B/C. ECDFs use log10(1 + distance in bp), retaining zero. Distinct alleles at the same coordinate are different leads with zero physical separation. Physical matching does not establish shared causal signals or LD.",
    sprintf("C. Absolute lead-to-gene-TSS distance for precisely the two extra-CS groups in B, excluding %s missing distances. Model-only regions remain excluded. The two ECDFs are descriptive, unmatched groups; TSS proximity is not a ground-truth accuracy measure.",
      sum(z$tss_distance_summary$n_missing_tss)),
    sprintf("D. Posterior coding support from the selected iteration's saved weighted_fit_mix fits in the same cohort. For each component l with a reported CS, q(l,c) = sum_j alpha(l,j) over ALL predictors j with coding c (additive, recessive, dominant), including predictors outside the reported CS. Average q equally across reported-CS components within each tissue. No PIP threshold, hard lead-coding assignment, CS truncation/renormalization, or refitting is used. %s EM CSs from %s / %s cohort regions contribute. Regions with zero EM CSs have no component to average and remain in the region audit; tissues with zero EM CSs have undefined support, not zero support.",
      d$n_em_cs, d$n_regions_with_em_cs, d$n_cohort_regions),
    sprintf("D pooled support, weighted by CS count: additive %s; recessive %s; dominant %s. Additive probability remains in the denominator although only recessive and dominant probabilities are shown by tissue. Tissues are ordered by mean recessive+dominant posterior support. These are model posterior probabilities among reported-CS components, not the proportion of loci that truly have nonadditive effects.",
      fig1_pct(d$mean_additive), fig1_pct(d$mean_recessive), fig1_pct(d$mean_dominant)),
    "The selected iteration is resolved once from completion markers, manifest and priors. Latest completed does not itself establish convergence across EM iterations. Coding priors are used only for provenance checks and are not panel D. Saved fit counts, component indices, memberships, lead IDs/codings and PIPs must match the EM summary.",
    if (!is.null(z$preview_label)) z$preview_label)
}

write_figure1_em <- function(z, output_dir, show_plot = FALSE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- c("metrics", "cohort", "region_audit", "tissue_ranking", "selected_tissues", "tissue_cs_count_tables",
    "tissue_count_table", "matched_cs_pairs", "unpaired_cs", "lead_retention_summary", "lead_distance_summary",
    "tss_distance_summary", "orphan_cs_rows", "coding_cs_support", "coding_region_audit",
    "coding_tissue_support", "coding_overall_support", "priors_used")
  for (nm in tables) write.csv(z[[nm]], file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
  writeLines(fig1_caption(z), file.path(output_dir, "figure1_caption.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  saveRDS(z, file.path(output_dir, "figure1_plot_data.rds"))
  plot_figure1_em(z, file.path(output_dir, "figure1_em.pdf"))
  plot_figure1_em(z, file.path(output_dir, "figure1_em.png"))
  if (show_plot) tryCatch(plot_figure1_em(z), error = function(e)
    message("Files saved; on-screen preview failed: ", conditionMessage(e), ". Open figure1_em.pdf."))
  message("Figure 1 saved to: ", normalizePath(output_dir, winslash = "/"))
  invisible(output_dir)
}

run_figure1_em <- function(settings) {
  if (is.null(settings$plot_data_file)) {
    context <- em_analysis_context(settings$project_dir, settings$iteration)
    summaries <- em_load_summary(context)
    message("Building the shared cohort, four tissue tables and ECDFs...")
    result <- summarize_figure1_em(summaries$res, summaries$cs, context$priors, settings)
    support <- fig1_collect_coding_support(context, result$cohort, summaries$cs)
    result <- fig1_add_coding_summary(result, support)
    result$preview_label <- settings$preview_label
    result$provenance <- list(iteration = context$iteration, iteration_dir = context$iteration_dir,
      summary_metadata = summaries$metadata, generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      coding_method = "Full alpha row summed by coding; equal weight per reported-CS component")
  } else {
    result <- readRDS(settings$plot_data_file)
    if (!identical(result$schema, "figure1_em_tables_ecdf_v2"))
      stop("Older/incompatible Figure 1 bundle: rerun on RCC with plot_data_file = NULL to calculate posterior support.")
    if (as.character(settings$iteration) != "latest" && as.character(settings$iteration) != as.character(result$iteration))
      stop("Requested iteration differs from the exported plot data.")
    message("Replotting frozen iteration ", result$iteration, "; all cohort settings come from the export.")
  }
  output <- settings$output_dir
  if (is.null(output)) output <- file.path(settings$project_dir, "results_em", sprintf("iteration_%03d", result$iteration),
    "descriptive_results", "figure1", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  result$output_dir <- output
  write_figure1_em(result, output, isTRUE(settings$show_plot))
  message("Shared cohort: ", fig1_fmt(result$metrics$n_regions), " regions; count agreement ",
    fig1_pct(result$metrics$count_agreement_rate), "; matched-lead retention ",
    fig1_pct(result$metrics$lead_retention_rate), ".")
  invisible(result)
}

