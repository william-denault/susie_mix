# Figure 1 calculations and graphics. Base R except lazy GTEx reconstruction.
# Source em_analysis_utils.R and em_cs_agreement_utils.R first.
fig1_coding_colors <- c(additive = "#397DA2", recessive = "#C65A30",
                        dominant = "#168578", unknown = "#858585")
fig1_outcome_colors <- c(same_leads = "#BCD5D2", changed_leads = "#4A87A6",
  fewer_em = "#E2B581", susie_only = "#BA7544", more_em = "#8470A8",
  em_only = "#554472", unknown_leads = "#B9B9B9")
fig1_outcome_labels <- c(same_leads = "Same count + same leads", changed_leads = "Same count + changed leads",
  fewer_em = "Fewer EM CSs (both >0)", susie_only = "SuSiE only",
  more_em = "More EM CSs (both >0)", em_only = "EM only", unknown_leads = "Same count; leads unavailable")
fig1_fmt <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
fig1_members <- function(s) {
  if (is.na(s) || !nzchar(trimws(s))) return(character())
  unique(trimws(strsplit(s, ";", fixed = TRUE)[[1L]]))
}

summarize_figure1_em <- function(res, cs, priors, settings) {
  em_agreement_require(cs, c("cs_snps", "cs_size_snps"), "Figure 1 CS summary")
  base <- summarize_em_cs_agreement(res, cs, settings$association_threshold,
    settings$minimum_mean_reads, minimum_tissue_regions = 1L, em_priors = priors)
  audit <- base$region_audit
  r <- audit[audit$count_eligible, , drop = FALSE]
  if (!nrow(r)) stop("No regions pass the shared Figure 1 cohort filters.")
  r$outcome <- "unknown_leads"
  same <- r$n_add == r$n_em
  r$outcome[which(same & r$lead_status == "ok" & r$n_shared_leads == r$n_add)] <- "same_leads"
  r$outcome[which(same & r$lead_status == "ok" & r$n_shared_leads < r$n_add)] <- "changed_leads"
  r$outcome[r$n_em < r$n_add & r$n_em > 0] <- "fewer_em"
  r$outcome[r$n_em == 0] <- "susie_only"
  r$outcome[r$n_em > r$n_add & r$n_add > 0] <- "more_em"
  r$outcome[r$n_add == 0] <- "em_only"
  categories <- names(fig1_outcome_colors)
  overall <- data.frame(outcome = categories, n_regions = as.integer(table(factor(r$outcome, categories))))
  overall$fraction <- overall$n_regions / nrow(r)
  tissue <- as.data.frame(table(tissue = r$tissue, outcome = factor(r$outcome, categories)))
  names(tissue)[3L] <- "n_regions"
  tissue$denominator <- ave(tissue$n_regions, tissue$tissue, FUN = sum)
  tissue$fraction <- tissue$n_regions / tissue$denominator
  counts <- as.data.frame(table(tissue = r$tissue, n_add = r$n_add, n_em = r$n_em))
  names(counts)[4L] <- "n_regions"
  counts <- counts[counts$n_regions > 0, ]

  # Only the CSs in this parent cohort are carried into panels B/C and the audit.
  cs <- as.data.frame(cs, stringsAsFactors = FALSE)
  cs$region_id <- em_agreement_key(cs$gene, cs$tissue)
  cs <- cs[cs$region_id %in% r$region_id & cs$model_key %in% c("susie_add", "weighted_fit_mix"), ]
  cs$cs_id <- paste(cs$region_id, cs$model_key, cs$cs_number, sep = "|")
  if (anyDuplicated(cs$cs_id)) stop("Duplicate CS IDs in the Figure 1 cohort.")
  members <- lapply(cs$cs_snps, fig1_members)
  cs$membership_status <- vapply(seq_len(nrow(cs)), function(i) {
    ids <- members[[i]]
    if (!length(ids) || !all(em_agreement_coordinates(ids)$valid)) return("missing_or_invalid_members")
    if (!is.finite(cs$cs_size_snps[i]) || length(ids) != cs$cs_size_snps[i]) return("size_mismatch")
    if (!cs$lead_snp[i] %in% ids) return("lead_absent_from_members")
    "ok"
  }, character(1))
  p <- base$matched_cs_pairs
  a <- match(p$additive_cs_id, cs$cs_id); m <- match(p$em_cs_id, cs$cs_id)
  if (anyNA(c(a, m))) stop("Matched CS absent from cohort summary.")
  p$additive_lead_pip <- cs$lead_pip[a]
  p$additive_size_snps <- cs$cs_size_snps[a]; p$em_size_snps <- cs$cs_size_snps[m]
  p$em_lead_coding[!p$em_lead_coding %in% names(fig1_coding_colors)] <- "unknown"
  p$valid_pip <- is.finite(p$em_lead_pip) & p$em_lead_pip >= 0 & p$em_lead_pip <= 1
  p$jaccard <- rep(NA_real_, nrow(p))
  for (i in which(cs$membership_status[a] == "ok" & cs$membership_status[m] == "ok"))
    p$jaccard[i] <- length(intersect(members[[a[i]]], members[[m[i]]])) /
      length(union(members[[a[i]]], members[[m[i]]]))
  p$overlap <- ifelse(is.na(p$jaccard), "Unknown", ifelse(p$jaccard == 1, "Identical",
    ifelse(p$jaccard == 0, "Disjoint", "Partial")))
  p$add_abs_tss_bp <- abs(cs$distance_to_tss_kb[a]) * 1000
  p$em_abs_tss_bp <- abs(cs$distance_to_tss_kb[m]) * 1000
  p$tss_status <- ifelse(is.finite(p$add_abs_tss_bp) & is.finite(p$em_abs_tss_bp),
                          "ok", "missing_distance")
  # Reject stale workhorse distances attached to a different component lead.
  for (idx in list(a, m)) {
    if (all(c("distance_source", "workhorse_lead_snp") %in% names(cs))) {
      bad <- grepl("workhorse", cs$distance_source[idx], ignore.case = TRUE) &
        (is.na(cs$workhorse_lead_snp[idx]) | cs$workhorse_lead_snp[idx] != cs$lead_snp[idx])
      p$tss_status[which(bad)] <- "workhorse_lead_mismatch"
    }
    if ("chromosome_matches_tss" %in% names(cs))
      p$tss_status[which(cs$chromosome_matches_tss[idx] %in% FALSE)] <- "tss_chromosome_mismatch"
    if ("lead_position" %in% names(cs)) {
      pos <- em_agreement_coordinates(cs$lead_snp[idx])$pos
      p$tss_status[which(is.finite(cs$lead_position[idx]) & cs$lead_position[idx] != pos)] <- "lead_coordinate_mismatch"
    }
    if ("tss_position" %in% names(cs)) {
      expected <- abs(em_agreement_coordinates(cs$lead_snp[idx])$pos - cs$tss_position[idx])
      observed <- abs(cs$distance_to_tss_kb[idx]) * 1000
      p$tss_status[which(is.finite(expected) & is.finite(observed) & abs(expected - observed) > 1)] <- "tss_distance_mismatch"
    }
  }
  if ("tss_position" %in% names(cs)) {
    bad <- is.finite(cs$tss_position[a]) & is.finite(cs$tss_position[m]) & cs$tss_position[a] != cs$tss_position[m]
    p$tss_status[which(bad)] <- "different_tss_annotations"
  }
  p$tss_delta_bp <- ifelse(p$tss_status == "ok", p$em_abs_tss_bp - p$add_abs_tss_bp, NA_real_)
  p$panel_b <- !p$same_lead_snp & p$valid_pip
  p$panel_c <- !p$same_lead_snp & p$tss_status == "ok"
  changed <- p[!p$same_lead_snp, , drop = FALSE]
  tc <- p[p$panel_c, , drop = FALSE]
  metrics <- data.frame(n_regions = nrow(r), n_genes = length(unique(r$gene)), n_tissues = length(unique(r$tissue)),
    n_add_cs = sum(r$n_add), n_em_cs = sum(r$n_em), n_same_count = sum(same),
    n_matching_unavailable_regions = sum(r$matching_status != "matched"), n_pairs = nrow(p),
    n_same_lead = sum(p$same_lead_snp), n_changed_lead = nrow(changed),
    n_changed_zero_bp = sum(changed$distance_bp == 0), n_panel_b = sum(p$panel_b),
    n_changed_missing_pip = sum(!changed$valid_pip),
    n_changed_over10kb = sum(changed$distance_bp > 1e4),
    n_changed_over100kb = sum(changed$distance_bp > 1e5), n_panel_c = nrow(tc),
    n_changed_missing_tss = sum(changed$tss_status != "ok"),
    n_tss_closer = sum(tc$tss_delta_bp < 0), n_tss_equal = sum(tc$tss_delta_bp == 0),
    n_tss_farther = sum(tc$tss_delta_bp > 0),
    median_tss_delta_kb = if (nrow(tc)) median(tc$tss_delta_bp) / 1000 else NA_real_)
  list(schema = "figure1_em_v1", iteration = base$iteration, settings = settings, metrics = metrics,
    overall = overall, tissue_outcomes = tissue, tissue_count_table = counts,
    cohort = r, region_audit = audit, matched_cs_pairs = p, unpaired_cs = base$unpaired_cs,
    membership_audit = cs[c("gene", "tissue", "model_key", "cs_number", "membership_status")],
    priors = priors)
}

fig1_reconstruct_example <- function(context, tissue_result, settings) {
  needed <- c("data.table", "matrixStats")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("GTEx reconstruction requires installed R packages: ", paste(absent, collapse = ", "))
  source(file.path(context$project_dir, "script/analysis/fit_mix_expression_plot_utils.R"), local = TRUE)
  d <- settings$datadir
  files <- c(subject_pheno_file = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz",
    sample_attr_file = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz",
    expr_file = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
    gtf_file = "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
  paths <- setNames(file.path(d, files), names(files))
  if (any(!file.exists(paths))) stop("Missing GTEx inputs: ", paste(paths[!file.exists(paths)], collapse = ", "))
  plink <- file.path(d, "plink2")
  if (!file.exists(plink)) stop("PLINK2 executable not found: ", plink)
  shared <- do.call(load_fit_mix_expression_inputs, c(as.list(paths), list(gene_annotation_function =
    file.path(context$project_dir, "script/scan_tissue_attempt/get_gene_annotations.R"))))
  gene <- prepare_fit_mix_gene_data(settings$example_gene, shared, context$baseline_dir,
    temp_dir = tempdir(), plink_exec = plink,
    geno_file = file.path(d, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"),
    cis_window = 5e5, min_maf_plink = 0, min_maf = .05, hwe_thresh = 1e-8)
  # The shared reconstruction helper calls this slot susie_mix. Explicitly use EM.
  tissue_result$susie_mix <- tissue_result$weighted_fit_mix
  prepare_fit_mix_tissue_data(gene, settings$example_gene, settings$example_tissue,
    tissue_result, min_samples = 0L, min_n_rec = 5L)
}

fig1_example <- function(context, summary_cs, result, settings, reconstruct = fig1_reconstruct_example) {
  id <- em_agreement_key(settings$example_gene, settings$example_tissue)
  row <- result$cohort[result$cohort$region_id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop(settings$example_gene, " / ", settings$example_tissue,
    " is outside the shared cohort. Filters will not be relaxed for panel D.")
  source(file.path(context$project_dir, "script/analysis/generate_summary_results.R"), local = TRUE)
  join <- em_join_gene; environment(join) <- environment()
  x <- join(file.path(context$iteration_dir, "results", paste0(settings$example_gene, ".rds")), context)[[settings$example_tissue]]
  if (is.null(x)) stop("Example is absent from the selected EM result.")
  leads <- points <- list()
  for (model in c("susie_add", "weighted_fit_mix")) {
    fit <- x[[model]]
    map <- as.data.frame(x[[if (model == "susie_add") "add_predictor_map" else "mix_predictor_map"]])
    map <- map[order(map$predictor_index), , drop = FALSE]
    validate_predictor_map(fit, map, paste("Example", model))
    sets <- get_cs(fit)
    if (length(sets) != if (model == "susie_add") row$n_add else row$n_em)
      stop("Example CS count differs between summary and saved fit; regenerate EM summaries.")
    if (length(fit$pip) != nrow(map) || any(!is.finite(fit$pip) | fit$pip < 0 | fit$pip > 1))
      stop("Invalid saved example predictor PIPs.")
    coords <- em_agreement_coordinates(map$snp)
    if (any(!coords$valid) || length(unique(coords$chr)) != 1L) stop("Invalid example predictor coordinates.")
    points[[model]] <- data.frame(model = model, snp = map$snp, coding = map$coding,
      position = coords$pos, chromosome = coords$chr, pip = fit$pip)
    for (i in seq_along(sets)) {
      k <- get_cs_lead_predictor(fit, sets, i)
      ids <- unique(map$snp[sets[[i]]])
      s <- summary_cs[summary_cs$gene == settings$example_gene & summary_cs$tissue == settings$example_tissue &
                        summary_cs$model_key == model & summary_cs$cs_number == i, , drop = FALSE]
      if (nrow(s) != 1L || !isTRUE(s$lead_snp == map$snp[k]) || !isTRUE(s$lead_coding == map$coding[k]) ||
          !isTRUE(all.equal(as.numeric(s$lead_pip), as.numeric(fit$pip[k]), tolerance = 1e-8)) ||
          !setequal(fig1_members(s$cs_snps), ids))
        stop("Example CS summary disagrees with saved ", model, " CS ", i, "; regenerate summaries.")
      leads[[length(leads) + 1L]] <- data.frame(model = model, cs_number = i, predictor_index = k,
        predictor_name = map$predictor_name[k], snp = map$snp[k], coding = map$coding[k], pip = fit$pip[k],
        size_snps = length(ids), position = coords$pos[k],
        left = min(em_agreement_coordinates(ids)$pos), right = max(em_agreement_coordinates(ids)$pos))
    }
  }
  leads <- do.call(rbind, leads); points <- do.call(rbind, points)
  em <- leads[leads$model == "weighted_fit_mix", , drop = FALSE]
  choices <- em[em$coding == settings$example_coding, , drop = FALSE]
  if (!nrow(choices)) stop("Example has no ", settings$example_coding, " EM lead in iteration ", context$iteration, ".")
  shared <- choices$snp %in% leads$snp[leads$model == "susie_add"]
  focal <- choices[order(!shared, -choices$pip, choices$cs_number)[1L], , drop = FALSE]
  message("Reconstructing ", settings$example_gene, " / ", settings$example_tissue, " from GTEx inputs...")
  d <- reconstruct(context, x, settings)
  expected <- as.data.frame(x$mix_predictor_map)
  expected <- expected[order(expected$predictor_index), ]
  if (!identical(colnames(d$geno_mix), as.character(expected$predictor_name)) ||
      !isTRUE(length(d$y) == x$n_ind) || nrow(d$geno_mix) != length(d$y) ||
      nrow(d$geno_for_counts) != length(d$y) || any(!is.finite(d$y)) ||
      any(!is.finite(d$geno_mix))) stop("Reconstructed example sample/predictor mismatch.")
  if (!focal$snp %in% colnames(d$geno_for_counts)) stop("Focal SNP absent from reconstructed genotypes.")
  g <- as.numeric(d$geno_for_counts[, focal$snp])
  if (anyNA(g) || any(!g %in% 0:2) || !all(0:2 %in% g)) stop("Example requires observed genotype groups 0, 1 and 2.")
  focal_code <- switch(focal$coding, additive = g, recessive = as.numeric(g == 2), dominant = as.numeric(g >= 1))
  if (!isTRUE(all.equal(unname(d$geno_mix[, focal$predictor_index]), focal_code)))
    stop("Focal coding does not agree with reconstructed minor-allele dosage.")
  # Same conditional-expression convention as existing example scripts: fit all
  # EM CS leads jointly, then subtract the other leads' fitted contributions.
  other <- setdiff(unique(em$predictor_index), focal$predictor_index)
  nuisance <- d$geno_mix[, other, drop = FALSE]
  design <- cbind(intercept = 1, focal = focal_code, nuisance)
  joint <- lm.fit(design, d$y)
  if (joint$rank != ncol(design) || any(!is.finite(joint$coefficients)))
    stop("Conditional example expression has a rank-deficient lead design.")
  y <- as.numeric(d$y - if (length(other)) nuisance %*% joint$coefficients[-c(1L, 2L)] else 0)
  fit_a <- lm.fit(cbind(1, g), y)$coefficients
  fit_r <- lm.fit(cbind(1, as.numeric(g == 2)), y)$coefficients
  groups <- do.call(rbind, lapply(0:2, function(i) data.frame(dosage = i, n = sum(g == i),
    mean = mean(y[g == i]), sd = sd(y[g == i]))))
  groups$additive_ols <- fit_a[1L] + fit_a[2L] * groups$dosage
  groups$recessive_ols <- fit_r[1L] + fit_r[2L] * as.numeric(groups$dosage == 2)
  list(gene = settings$example_gene, tissue = settings$example_tissue, n = length(y),
    focal = focal, shared_lead = focal$snp %in% leads$snp[leads$model == "susie_add"],
    leads = leads, regional = points, groups = groups,
    expression = data.frame(dosage = g, transformed_expression = as.numeric(d$y), conditional_expression = y),
    adjusted_for = expected$predictor_name[other], conditioning_coefficients = setNames(joint$coefficients, colnames(design)),
    region = row, chromosome = unique(points$chromosome))
}

fig1_distance_axis <- function(side, max_bp) {
  ticks <- c(0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
  labels <- c("0", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  keep <- ticks <= max_bp
  axis(side, at = log10(1 + ticks[keep]), labels = labels[keep], las = 1, cex.axis = .8)
}

plot_figure1_em <- function(z, file = NULL) {
  if (!identical(z$schema, "figure1_em_v1") || is.null(z$example)) stop("Incomplete Figure 1 bundle.")
  if (!is.null(file)) {
    if (tolower(tools::file_ext(file)) == "pdf") pdf(file, width = 15, height = 12, useDingbats = FALSE)
    else png(file, width = 3750, height = 3000, res = 250,
      type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    on.exit(dev.off(), add = TRUE)
  }
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE, after = FALSE)
  drawing_scale <- if (is.null(file)) min(1, dev.size("in") / c(15, 12)) else 1
  par(family = "sans", col.axis = "#454545", col.lab = "#333333", cex = drawing_scale, mex = drawing_scale,
      xaxs = "i", yaxs = "i", las = 1, mgp = c(2.0, .65, 0))
  # Fixed page coordinates keep the four panel titles aligned across exports.
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4))
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(.5, .976, "Fine-mapping under additive and mixed genetic models", font = 2, cex = 1.55)
  m <- z$metrics; ex <- z$example; p <- z$matched_cs_pairs
  subtitle <- sprintf("EM iteration %s  |  %s gene-tissue regions  |  %s genes  |  %s tissues",
    z$iteration, fig1_fmt(m$n_regions), fig1_fmt(m$n_genes), m$n_tissues)
  text(.5, .949, subtitle, cex = 1.02, col = "#555555")
  text(.018, .914, "A", adj = 0, font = 2, cex = 1.4)
  text(.048, .914, "Signal architecture across tissues", adj = 0, font = 2, cex = 1.15)
  text(.52, .914, "B", adj = 0, font = 2, cex = 1.4)
  text(.55, .914, "Lead changes and posterior support", adj = 0, font = 2, cex = 1.15)
  text(.018, .464, "C", adj = 0, font = 2, cex = 1.4)
  text(.048, .464, "Do lead changes move toward the TSS?", adj = 0, font = 2, cex = 1.15)
  text(.52, .464, "D", adj = 0, font = 2, cex = 1.4)
  text(.55, .464, paste0(ex$gene, ": ", if (ex$shared_lead) "shared variant, different coding" else "genotype and expression"),
       adj = 0, font = 2, cex = 1.15)
  text(.5, .026, sprintf("Shared cohort: additive P < %s; mean reads >= %s; >=1 CS in either model. B/C: matched CS subsets; D: one included region.",
    format(z$settings$association_threshold, scientific = TRUE), z$settings$minimum_mean_reads), cex = .78, col = "#555555")
  if (!is.null(z$preview_label)) text(.5, .009, z$preview_label, font = 2, col = "#B32929", cex = .9)

  # A: one outcome per region, with every tissue represented and no ranking cutoff.
  par(fig = c(.005, .495, .578, .893), mar = c(1.8, 7.9, .15, 3.1), new = TRUE)
  t <- z$tissue_outcomes
  tids <- unique(as.character(t$tissue))
  agree <- vapply(tids, function(s) sum(t$fraction[t$tissue == s & t$outcome %in% c("same_leads", "changed_leads", "unknown_leads")]), numeric(1))
  tids <- tids[order(agree, tids)]
  n <- length(tids)
  plot.new(); plot.window(c(0, 1), c(0, n + 3))
  abline(v = c(.25, .5, .75, 1), col = "#EEEEEE")
  bar <- function(d, y, h) {
    right <- cumsum(d$fraction); left <- c(0, head(right, -1))
    rect(left, y - h, right, y + h, border = NA, col = fig1_outcome_colors[as.character(d$outcome)])
  }
  bar(z$overall, n + 2.1, .48)
  text(0, n + 2.1, "All tissues", pos = 2, xpd = NA, font = 2, cex = .83)
  text(1, n + 2.1, fig1_fmt(m$n_regions), pos = 4, xpd = NA, font = 2, cex = .78)
  for (i in seq_along(tids)) {
    d <- t[as.character(t$tissue) == tids[i], ]
    y <- n - i + 1
    bar(d, y, .34)
    text(0, y, tids[i], pos = 2, xpd = NA, cex = min(.83, 20 / n))
    text(1, y, fig1_fmt(unique(d$denominator)), pos = 4, xpd = NA, cex = min(.77, 19 / n))
  }
  axis(1, at = seq(0, 1, .25), labels = paste0(seq(0, 100, 25), "%"), cex.axis = .8)
  par(fig = c(.02, .49, .493, .574), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  active <- z$overall$outcome[z$overall$n_regions > 0]
  labels <- sprintf("%s (%.1f%%)", fig1_outcome_labels[active], 100 * z$overall$fraction[match(active, z$overall$outcome)])
  legend("topleft", legend = labels, fill = fig1_outcome_colors[active], border = NA,
         bty = "n", ncol = 2, cex = .78, x.intersp = .6, y.intersp = 1.1)
  text(.01, .15, "Rows ordered by CS-count agreement; numbers give region counts.", adj = 0, cex = .75, col = "#555555")

  # B: all paired leads in the retention strip; changed leads in the scatter.
  par(fig = c(.56, .975, .82, .897), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  retained <- if (m$n_pairs) m$n_same_lead / m$n_pairs else 0
  rect(0, .45, retained, .67, col = fig1_outcome_colors["same_leads"], border = NA)
  rect(retained, .45, 1, .67, col = fig1_outcome_colors["changed_leads"], border = NA)
  text(0, .93, sprintf("%s / %s paired leads retained (%.1f%%)", fig1_fmt(m$n_same_lead),
    fig1_fmt(m$n_pairs), 100 * retained), adj = 0, font = 2, cex = .9)
  text(0, .13, sprintf("%s changed leads; %s shown with a valid PIP", fig1_fmt(m$n_changed_lead), fig1_fmt(m$n_panel_b)), adj = 0, cex = .8)
  par(fig = c(.515, .985, .57, .819), mar = c(3, 4.3, .4, .8), new = TRUE)
  b <- p[p$panel_b, , drop = FALSE]
  max_bp <- max(c(1e5, b$distance_bp)) * 1.08
  plot(NA, xlim = c(0, log10(1 + max_bp)), ylim = c(-.02, 1.09), axes = FALSE,
       xlab = "Lead-SNP separation (log1p scale)", ylab = "EM lead PIP (coding-specific)", cex.lab = .9)
  fig1_distance_axis(1, max_bp); axis(2, at = c(0, .5, 1), cex.axis = .8)
  abline(h = c(0, .5, 1), col = "#E8E8E8", lty = c(1, 2, 1))
  abline(v = log10(1 + c(1e4, 1e5)), col = "#E8E8E8", lty = 2)
  shapes <- c(Identical = 16, Partial = 17, Disjoint = 4, Unknown = 1)
  if (nrow(b)) points(log10(1 + b$distance_bp), b$em_lead_pip, pch = shapes[b$overlap],
                     col = adjustcolor(fig1_coding_colors[b$em_lead_coding], .7), cex = .8)
  # At most three distinct genes; explicit simple placement avoids label crowds.
  label <- b[b$em_lead_pip >= .9 & b$distance_bp > 1e4, , drop = FALSE]
  label <- label[order(-label$em_lead_pip, -label$distance_bp, label$gene), ]
  label <- head(label[!duplicated(label$gene), ], 3L)
  if (nrow(label)) {
    x <- log10(1 + label$distance_bp); y <- label$em_lead_pip
    ly <- seq(1.055, .86, length.out = nrow(label))
    lx <- pmax(.3, pmin(log10(1 + max_bp) - .65, x - .5))
    segments(x, y, lx, ly, col = "#888888", lwd = .6)
    text(lx, ly, label$gene, cex = .7)
  }
  if (!nrow(b)) text(log10(1 + max_bp) / 2, .6, "No changed leads with a valid PIP")
  par(fig = c(.55, .98, .49, .565), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  cl <- if (any(p$em_lead_coding == "unknown")) 1:4 else 1:3
  legend("topleft", c("Additive", "Recessive", "Dominant", "Unknown coding")[cl], col = fig1_coding_colors[cl], pch = 16, ncol = length(cl),
    bty = "n", cex = .8, x.intersp = .5)
  legend(.0, .57, paste(c("Identical", "Partial", "Disjoint", "Unknown"), "CSs"),
    pch = shapes, ncol = 4, bty = "n", cex = .73, x.intersp = .5)
  text(0, .05, sprintf("%s unpaired CSs; %s regions with unavailable matching (audited).",
    fig1_fmt(nrow(z$unpaired_cs)), fig1_fmt(m$n_matching_unavailable_regions)), adj = 0, cex = .73, col = "#555555")

  # C: paired distances for the SAME changed-lead population, not two unmatched distributions.
  par(fig = c(.015, .495, .108, .436), mar = c(3.4, 5.3, 2.3, .8), new = TRUE)
  cdat <- p[p$panel_c, , drop = FALSE]
  max_tss <- max(c(1e5, cdat$add_abs_tss_bp, cdat$em_abs_tss_bp)) * 1.12
  lim <- c(0, log10(1 + max_tss))
  plot(NA, xlim = lim, ylim = lim, axes = FALSE, xlab = "SuSiE |lead - TSS| (log1p scale)",
       ylab = "", cex.lab = .9)
  mtext("EM |lead - TSS| (log1p scale)", side = 2, line = 3.6, las = 0, cex = .9 * drawing_scale)
  fig1_distance_axis(1, max_tss); fig1_distance_axis(2, max_tss)
  abline(0, 1, col = "#666666", lty = 2)
  if (nrow(cdat)) points(log10(1 + cdat$add_abs_tss_bp), log10(1 + cdat$em_abs_tss_bp),
    col = adjustcolor(fig1_coding_colors[cdat$em_lead_coding], .65), pch = shapes[cdat$overlap], cex = .85)
  text(lim[2] * .04, lim[2] * .91, "EM farther from TSS", adj = 0, col = "#666666", cex = .8)
  text(lim[2] * .96, lim[2] * .10, "EM closer to TSS", adj = 1, col = "#666666", cex = .8)
  title(sub = NULL, main = sprintf("%s changed pairs with both TSS distances (%s unavailable)",
    fig1_fmt(m$n_panel_c), fig1_fmt(m$n_changed_missing_tss)), adj = 0, font.main = 1, cex.main = .85)
  par(fig = c(.06, .49, .055, .10), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0, .8, sprintf("Closer: %s   |   Equal: %s   |   Farther: %s   |   Median change: %s kb",
    m$n_tss_closer, m$n_tss_equal, m$n_tss_farther, format(round(m$median_tss_delta_kb, 1), trim = TRUE)), adj = 0, cex = .8)
  text(0, .28, "Colors and shapes as in B. Proximity to TSS does not establish causal accuracy.", adj = 0, cex = .74, col = "#555555")

  # D: saved predictor PIPs and conditional expression from exactly this cohort region.
  par(fig = c(.52, .985, .287, .434), mar = c(2.9, 4.1, 1.9, .8), new = TRUE)
  span <- range(c(ex$leads$left, ex$leads$right))
  mid <- mean(span); half <- max(diff(span) / 2 + 1e4, 2.5e4)
  span <- c(mid - half, mid + half)
  pp <- ex$regional[ex$regional$position >= span[1] & ex$regional$position <= span[2], ]
  plot(NA, xlim = span / 1e6, ylim = c(0, 1.1), xlab = paste0("Chromosome ", ex$chromosome, " position (Mb, GRCh38)"),
    ylab = "Predictor PIP", axes = FALSE, cex.lab = .83)
  axis(1, cex.axis = .74); axis(2, at = c(0, .5, 1), cex.axis = .74)
  abline(h = .5, col = "#E6E6E6", lty = 2)
  abline(v = ex$focal$position / 1e6, col = "#D0D0D0", lty = 2)
  for (model in c("susie_add", "weighted_fit_mix")) {
    d <- pp[pp$model == model, ]
    cols <- if (model == "susie_add") rep("#777777", nrow(d)) else unname(fig1_coding_colors[d$coding])
    points(d$position / 1e6, d$pip, col = adjustcolor(cols, .65),
      pch = if (model == "susie_add") 1 else 16, cex = if (model == "susie_add") .75 else .65)
  }
  title(main = sprintf("%s | %s samples | %s SuSiE CS / %s EM CSs", ex$tissue,
    fig1_fmt(ex$n), ex$region$n_add, ex$region$n_em), adj = 0, font.main = 1, cex.main = .85)
  legend("topright", c("SuSiE", "EM A", "EM R", "EM D"), pch = c(1, 16, 16, 16),
    col = c("#777777", fig1_coding_colors[1:3]), ncol = 4, bty = "n", cex = .63, x.intersp = .5)
  par(fig = c(.52, .985, .096, .28), mar = c(2.8, 4.1, 2.1, .8), new = TRUE)
  d <- ex$expression; g <- ex$groups
  yr <- range(c(d$conditional_expression, g$additive_ols, g$recessive_ols)); pad <- max(.15, diff(yr) * .08)
  plot(NA, xlim = c(-.4, 2.4), ylim = yr + c(-pad, pad), axes = FALSE,
    xlab = "Minor-allele dosage", ylab = "Conditional expression", cex.lab = .83)
  axis(1, at = 0:2, labels = paste0(0:2, "  (n=", g$n, ")"), cex.axis = .77)
  axis(2, cex.axis = .74)
  # Deterministic anonymous jitter: does not alter the user's random-number state.
  jitter <- ((seq_len(nrow(d)) * .61803398875) %% 1 - .5) * .24
  points(d$dosage + jitter, d$conditional_expression, pch = 16, col = adjustcolor("#777777", .3), cex = .55)
  lines(0:2, g$additive_ols, col = fig1_coding_colors["additive"], lwd = 2, lty = 2)
  lines(0:2, g$recessive_ols, col = fig1_coding_colors["recessive"], lwd = 2)
  points(0:2, g$mean, pch = 21, bg = "white", col = "#222222", cex = 1.05, lwd = 1.3)
  title(main = sprintf("%s | EM %s PIP %.3f", sub("_b38.*$", "", ex$focal$snp),
    ex$focal$coding, ex$focal$pip), adj = 0, font.main = 1, cex.main = .84)
  legend("topleft", c("Additive OLS", "Recessive OLS", "Group mean"),
    col = c(fig1_coding_colors[c("additive", "recessive")], "#222222"), lty = c(2, 1, NA),
    pch = c(NA, NA, 21), pt.bg = "white", ncol = 3, bty = "n", cex = .67, x.intersp = .5)
  par(fig = c(.55, .985, .047, .09), mar = rep(0, 4), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  extra <- ex$leads[ex$leads$model == "weighted_fit_mix" & ex$leads$cs_number != ex$focal$cs_number, ]
  line <- if (nrow(extra)) paste(sprintf("CS %s: %s SNPs, lead PIP %.3f", extra$cs_number,
    extra$size_snps, extra$pip), collapse = "; ") else "No other EM CS in this region."
  text(0, .83, paste0("Other EM CSs: ", line), adj = 0, cex = .73)
  text(0, .31, "Expression adjusted for other EM CS leads; OLS curves are descriptive.", adj = 0, cex = .73, col = "#555555")
  invisible(z)
}

fig1_caption <- function(z) {
  m <- z$metrics
  c("Figure 1. Fine-mapping under additive and mixed genetic models.",
    sprintf("One parent cohort: %s gene-tissue regions, %s genes, %s tissues; EM iteration %s. Minimum marginal additive P < %s, mean read count >= %s, and at least one CS in either model.",
      m$n_regions, m$n_genes, m$n_tissues, z$iteration, z$settings$association_threshold, z$settings$minimum_mean_reads),
    "A. Mutually exclusive region outcomes, pooled and for all eligible tissues. Same leads means the same multiset of biological variant IDs, ignoring coding; it does not imply identical CS membership. Tissue order uses CS-count agreement, not expression sample size.",
    sprintf("B. %s one-to-one CS pairs: %s identical biological leads and %s changed leads; %s changed pairs have valid coding-specific EM PIPs. Reserve identical lead IDs first, then minimize total physical distance among remaining CSs within each region. Extra CSs are not forced into pairs: %s are exported separately. Matching is descriptive and is not evidence of shared causal signals or LD.",
      m$n_pairs, m$n_same_lead, m$n_changed_lead, m$n_panel_b, nrow(z$unpaired_cs)),
    "B/C. Colors denote the EM lead's coding. Shapes denote biological-SNP CS membership: identical (Jaccard=1), partial (0<Jaccard<1), disjoint (Jaccard=0), or unknown. Coding-specific PIP is not the posterior probability that a CS exists. Log1p axes retain zero distance, including different alleles at the same position.",
    sprintf("C. The changed-lead pairs in B's parent population, paired by gene and tissue; %s have two usable TSS distances and %s do not. EM minus SuSiE absolute TSS distance: %s closer, %s equal, %s farther; median %.3f kb. This panel has no PIP threshold. TSS proximity alone is not evidence of improved causal localization.",
      m$n_panel_c, m$n_changed_missing_tss, m$n_tss_closer, m$n_tss_equal, m$n_tss_farther, m$median_tss_delta_kb),
    sprintf("D. %s / %s, an included cohort region with %s samples. Regional points are saved predictor PIPs, not summed across codings. The window includes all CS members plus padding. The focal %s lead is %s, PIP %.6f. Genotypes are minor-allele dosage from the original reconstruction pipeline.",
      z$example$gene, z$example$tissue, z$example$n, z$example$focal$coding, z$example$focal$snp, z$example$focal$pip),
    "D. Reconstructed phenotype uses the existing library-size adjustment, log1p, rank-normalization, sex residualization and second rank-normalization pipeline. For conditional expression, regress this phenotype jointly on the EM CS lead predictors; subtract contributions of all other EM leads. Lines are additive and recessive ordinary least-squares fits to this same adjusted phenotype; circles are observed group means. This selected example illustrates coding, not independent validation or a formal model-comparison test.",
    "EM coding priors were learned on the full EM fitting data, not re-estimated in this display cohort. Latest means latest completed iteration verified by manifests, markers and priors; convergence across EM iterations must be established separately. The selected iteration and summary provenance are saved with the figure.",
    if (!is.null(z$preview_label)) z$preview_label)
}

write_figure1_em <- function(z, output_dir, show_plot = FALSE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- c("metrics", "overall", "tissue_outcomes", "tissue_count_table", "cohort", "region_audit",
              "matched_cs_pairs", "unpaired_cs", "membership_audit", "priors")
  for (nm in tables) write.csv(z[[nm]], file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
  for (nm in c("leads", "regional", "groups", "expression"))
    write.csv(z$example[[nm]], file.path(output_dir, paste0("example_", nm, ".csv")), row.names = FALSE)
  writeLines(fig1_caption(z), file.path(output_dir, "figure1_caption.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  saveRDS(z, file.path(output_dir, "figure1_plot_data.rds"))
  plot_figure1_em(z, file.path(output_dir, "figure1_em.pdf"))
  plot_figure1_em(z, file.path(output_dir, "figure1_em.png"))
  if (show_plot) tryCatch(plot_figure1_em(z), error = function(e)
    message("Files saved, but the on-screen preview failed: ", conditionMessage(e),
      ". Open figure1_em.pdf or enlarge the Plots pane and call plot_figure1_em(figure1_result)."))
  message("Figure 1 saved to: ", normalizePath(output_dir, winslash = "/"))
  invisible(output_dir)
}

run_figure1_em <- function(settings, reconstruct = fig1_reconstruct_example) {
  if (is.null(settings$plot_data_file)) {
    context <- em_analysis_context(settings$project_dir, settings$iteration)
    summaries <- em_load_summary(context)
    message("Building the shared Figure 1 cohort and matching CSs...")
    result <- summarize_figure1_em(summaries$res, summaries$cs, context$priors, settings)
    result$preview_label <- settings$preview_label
    result$example <- fig1_example(context, summaries$cs, result, settings, reconstruct)
    result$provenance <- list(iteration = context$iteration, iteration_dir = context$iteration_dir,
      summary_metadata = summaries$metadata, generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      example_baseline_file = file.path(context$baseline_dir, paste0(settings$example_gene, ".rds")),
      example_em_file = file.path(context$iteration_dir, "results", paste0(settings$example_gene, ".rds")),
      reconstruction = "fit_mix_expression_plot_utils.R; MAF > .05; HWE >= 1e-8; min_n_rec=5; cis_window=500kb")
  } else {
    result <- readRDS(settings$plot_data_file)
    if (!identical(result$schema, "figure1_em_v1")) stop("Unsupported Figure 1 plot-data file.")
    if (as.character(settings$iteration) != "latest" && as.character(settings$iteration) != as.character(result$iteration))
      stop("Requested iteration differs from the exported plot data.")
    message("Replotting frozen EM iteration ", result$iteration, "; cohort and example come from the export.")
  }
  output <- settings$output_dir
  if (is.null(output)) output <- file.path(settings$project_dir, "results_em", sprintf("iteration_%03d", result$iteration),
    "descriptive_results", "figure1", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  result$output_dir <- output
  write_figure1_em(result, output, settings$show_plot)
  message("Shared cohort: ", fig1_fmt(result$metrics$n_regions), " regions; ", result$metrics$n_tissues,
    " tissues. Example: ", result$example$gene, " / ", result$example$tissue, ".")
  invisible(result)
}
