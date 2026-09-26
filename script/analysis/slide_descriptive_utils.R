# Source slide_analysis_utils.R before this file.
slide_tss_distribution <- function(cs, tissues = NULL, limit_kb = 200, bandwidth_kb = 10) {
  rows <- list()
  for (tissue in c("All tissues", tissues)) for (model in c("SuSiE", "SuSiE-slide")) {
    z <- cs[cs$model == model & (tissue == "All tissues" | cs$tissue == tissue), , drop = FALSE]
    rows[[length(rows) + 1L]] <- summarize_tss_kde(
      z, model, tissue, plot_limit_kb = limit_kb, bandwidth_kb = bandwidth_kb)
  }
  do.call(rbind, rows)
}

run_slide_descriptive_results <- function(project_dir,
    summary_dir = file.path(project_dir, "results_slide/summary"),
    output_dir = file.path(project_dir, "results_slide/descriptive_results"),
    association_threshold = 1e-8, minimum_mean_reads = 100,
    tss_bandwidth_kb = 10, tss_plot_limit_kb = 200) {
  summaries <- slide_load_summary(summary_dir)
  res <- summaries$res
  primary <- slide_primary(res, association_threshold, minimum_mean_reads)
  cs <- summaries$cs
  cs <- cs[slide_region_key(cs$gene, cs$tissue) %in% slide_region_key(primary$gene, primary$tissue), , drop = FALSE]
  # Incomplete CS tables must not turn an unavailable comparator into zero CSs.
  for (model in c("susie_add", "fit_slide")) {
    counts <- table(slide_region_key(cs$gene[cs$model_key == model], cs$tissue[cs$model_key == model]))
    observed <- as.integer(counts[slide_region_key(primary$gene, primary$tissue)])
    observed[is.na(observed)] <- 0L
    expected <- primary[[if (model == "susie_add") "ncs_susie" else "ncs_slide"]]
    if (any(observed != expected)) stop("CS counts disagree with tissue summary; regenerate slide summaries together.")
  }
  selection <- select_tss_disagreement(cs, mixed_model = "SuSiE-slide",
                                        plot_limit_kb = tss_plot_limit_kb)
  tissues <- sort(unique(primary$tissue))
  tss <- slide_tss_distribution(selection$selected, tissues, limit_kb = tss_plot_limit_kb,
                                bandwidth_kb = tss_bandwidth_kb)
  coding_levels <- c("additive", "partial_recessive", "recessive", "partial_dominant", "dominant")
  coding_colors <- c("#4584AD", "#EBA77F", "#CC603D", "#AC9CC8", "#795294")
  summaries_by_tissue <- codings <- agreements <- list()
  pct <- function(x) if (length(x) && any(!is.na(x))) 100 * mean(x, na.rm = TRUE) else NA_real_
  for (tissue in c("All tissues", tissues)) {
    z <- primary[tissue == "All tissues" | primary$tissue == tissue, , drop = FALSE]
    s <- cs[cs$model_key == "fit_slide" & (tissue == "All tissues" | cs$tissue == tissue), , drop = FALSE]
    summaries_by_tissue[[length(summaries_by_tissue) + 1L]] <- data.frame(
      tissue = tissue, n_pairs = nrow(z), n_genes = length(unique(z$gene)),
      additive_cs = sum(z$ncs_susie), slide_cs = sum(z$ncs_slide),
      additive_pairs_with_cs = sum(z$ncs_susie > 0), slide_pairs_with_cs = sum(z$ncs_slide > 0),
      slide_only_pairs = sum(z$ncs_slide > 0 & z$ncs_susie == 0),
      additive_only_pairs = sum(z$ncs_susie > 0 & z$ncs_slide == 0),
      pct_slide_more_cs = pct(z$ncs_slide > z$ncs_susie),
      pct_same_lead_snp_set = pct(z$same_lead_snp_set),
      n_lead_agreement_comparable = sum(!is.na(z$same_lead_snp_set)),
      median_elbo_difference = if (any(is.finite(z$dif_elbo_slide_vs_add)))
        median(z$dif_elbo_slide_vs_add, na.rm = TRUE) else NA_real_,
      n_slide_not_converged = sum(z$converged_slide %in% FALSE),
      n_additive_not_converged = sum(z$converged_susie %in% FALSE))
    for (coding in coding_levels) codings[[length(codings) + 1L]] <- data.frame(
      tissue = tissue, coding = coding, n_cs = sum(s$lead_coding == coding),
      n_slide_cs = nrow(s), percentage = safe_percentage(sum(s$lead_coding == coding), nrow(s)))
    for (scope in c("all", "one_cs_in_both")) {
      a <- if (scope == "all") z else z[z$ncs_susie == 1L & z$ncs_slide == 1L, , drop = FALSE]
      agreements[[length(agreements) + 1L]] <- data.frame(tissue = tissue, scope = scope,
        n_pairs = nrow(a), n_with_leads_in_both = sum(!is.na(a$same_lead_snp_set)),
        n_same_biological_lead_set = sum(a$same_lead_snp_set %in% TRUE),
        n_different_biological_lead_set = sum(a$same_lead_snp_set %in% FALSE),
        pct_same_lead_set = pct(a$same_lead_snp_set),
        n_same_cs_snp_sets = sum(a$same_cs_snp_sets %in% TRUE),
        n_with_any_cs_overlap = sum(a$overlap_snp > 0))
    }
  }
  tissue_summary <- do.call(rbind, summaries_by_tissue)
  coding_summary <- do.call(rbind, codings)
  agreement_summary <- do.call(rbind, agreements)
  permutation <- slide_bind_rows(lapply(list(all = res, primary = primary), function(z) {
    paired <- z[which(z$has_susie_add_perm & z$has_fit_slide_perm), , drop = FALSE]
    data.frame(n_pairs = nrow(z), n_pairs_with_both_permutations = nrow(paired),
      additive_permutation_cs = sum(paired$perm_cs_susie), slide_permutation_cs = sum(paired$perm_cs_slide),
      additive_permutation_pairs_with_cs = sum(paired$perm_cs_susie > 0),
      slide_permutation_pairs_with_cs = sum(paired$perm_cs_slide > 0))
  }))
  permutation$subset <- c("all", "primary")
  cs_counts <- as.data.frame(table(additive_cs = primary$ncs_susie, slide_cs = primary$ncs_slide))
  cs_counts <- cs_counts[cs_counts$Freq > 0L, , drop = FALSE]
  overview <- data.frame(n_pairs_available = nrow(res), n_pairs_primary = nrow(primary),
    association_threshold = association_threshold, minimum_mean_reads = minimum_mean_reads,
    delta_tolerance = unique(res$delta_tolerance)[1],
    tss_agreement_rule = "same biological CS lead SNP within gene/tissue; coding ignored")
  tables <- list(overall_summary = overview, tissue_summary = tissue_summary,
    coding_summary = coding_summary, lead_agreement_summary = agreement_summary,
    cs_count_comparison = cs_counts, permutation_summary = permutation,
    primary_gene_tissue_pairs = primary, slide_cs_delta = cs[cs$model_key == "fit_slide", , drop = FALSE],
    tss_cs_agreement_audit = selection$audit, tss_cs_selection_summary = selection$summary,
    tss_distance_distribution = tss[tss$tissue == "All tissues", ],
    tissue_tss_distance_distribution = tss[tss$tissue != "All tissues", ])
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tables)) write.csv(tables[[nm]], file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
  plot_tss <- function(tissue) {
    plot_tss_kde(tss[tss$tissue == tissue, ],
                 model_colors = c(SuSiE = "#4584AD", `SuSiE-slide` = "#CC603D"),
                 plot_limit_kb = tss_plot_limit_kb,
                 main_title = "TSS distance: disagreeing CS leads", compact = TRUE)
  }
  plot_coding <- function(tissue) {
    z <- coding_summary[coding_summary$tissue == tissue, ]
    barplot(z$n_cs, names.arg = c("Additive", "Partial\nrecessive", "Recessive", "Partial\ndominant", "Dominant"),
      col = coding_colors, border = NA, ylab = "SuSiE-slide CS count", cex.names = .75,
      main = "Coding of the CS lead", ylim = c(0, max(1, z$n_cs) * 1.15), cex.main = .95)
  }
  draw_overview <- function() {
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(1.5, 0, 3, 0))
    z <- tissue_summary[1, ]
    barplot(c(z$additive_cs, z$slide_cs), names.arg = c("SuSiE", "SuSiE-slide"),
      col = c("#4584AD", "#CC603D"), border = NA, ylab = "Credible sets", main = "Primary analysis set",
      ylim = c(0, max(1, z$additive_cs, z$slide_cs) * 1.15))
    mtext(sprintf("%d gene-tissue pairs", nrow(primary)), side = 3, line = .2, cex = .8)
    a <- agreement_summary[agreement_summary$tissue == "All tissues" & agreement_summary$scope == "all", ]
    barplot(c(a$n_same_biological_lead_set, a$n_different_biological_lead_set,
               a$n_pairs - a$n_with_leads_in_both), names.arg = c("Same leads", "Different leads", "No CS in\none/both"),
      col = c("#468578", "#CC603D", "#ADB5BC"), border = NA, ylab = "Gene-tissue pairs",
      main = "Biological lead-SNP agreement", cex.names = .8,
      ylim = c(0, max(1, a$n_pairs) * 1.1))
    plot_coding("All tissues")
    plot_tss("All tissues")
    mtext("SuSiE-slide versus additive SuSiE", outer = TRUE, side = 3, line = 1, cex = 1.3, font = 2)
    mtext(sprintf("P < %.1e; mean reads >= %g. Coding uses component-specific delta.",
                  association_threshold, minimum_mean_reads), outer = TRUE, side = 1, cex = .75)
  }
  pdf(file.path(output_dir, "coding_and_tss_summary_2x2.pdf"), width = 12, height = 9)
  tryCatch(draw_overview(), finally = dev.off())
  png(file.path(output_dir, "coding_and_tss_summary_2x2.png"), width = 1800, height = 1350, res = 150)
  tryCatch(draw_overview(), finally = dev.off())
  pdf(file.path(output_dir, "tissue_specific_coding_and_tss_4x2.pdf"), width = 12, height = 14)
  tryCatch({
    if (!length(tissues)) { plot.new(); text(.5, .5, "No gene-tissue pairs pass the primary filters") }
    for (i in seq_along(tissues)) {
      if ((i - 1L) %% 4L == 0L) par(mfrow = c(4, 2), mar = c(3.5, 4.5, 3.5, 1))
      plot_coding(tissues[i]); mtext(tissues[i], side = 3, line = 2, font = 2, cex = .9)
      plot_tss(tissues[i])
    }
  }, finally = dev.off())
  print(tissue_summary[1, ])
  cat("TSS plots exclude shared biological CS leads (coding ignored):\n")
  print(selection$summary[selection$summary$tissue == "All tissues", ])
  message("Slide descriptive results: ", output_dir)
  invisible(tables)
}
