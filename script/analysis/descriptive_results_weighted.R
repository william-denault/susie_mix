# Descriptive analysis dedicated to weighted_fit_mix results.
#
# Run generate_summary_results.R first, then run this script. Outputs are
# written to descriptive_results_weighted/ and do not replace the original
# unweighted descriptive analysis.

library(data.table)

summary_file <- "/project2/mstephens/wdenault/susie_mix/res_summary.RData"
cs_summary_file <- "/project2/mstephens/wdenault/susie_mix/res_cs_summary.RData"
output_dir <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "descriptive_results_weighted/"
)

load(summary_file)
load(cs_summary_file)

if (!exists("res_summary")) {
  stop("res_summary.RData does not contain res_summary.")
}
if (!exists("res_cs_summary")) {
  stop("res_cs_summary.RData does not contain res_cs_summary.")
}

res <- as.data.table(res_summary)
cs_res <- as.data.table(res_cs_summary)

association_threshold <- 1e-8
minimum_mean_reads <- 100
tss_plot_limit_kb <- 200
tss_bin_width_kb <- 10

required_columns <- c(
  "gene",
  "tissue",
  "min_pv",
  "mean_count",
  "has_weighted_fit_mix",
  "ncs_susie",
  "ncs_susie_mix",
  "ncs_weighted_fit_mix",
  "n_add",
  "n_rec",
  "n_dom",
  "n_add_weighted",
  "n_rec_weighted",
  "n_dom_weighted",
  "overlap_snp_add_weighted",
  "mix_weighted_overlap_snp",
  "mix_weighted_same_cs_count",
  "mix_weighted_same_cs_snp_sets",
  "mix_weighted_same_lead_snp_set",
  "dif_elbo_weighted_vs_add",
  "dif_elbo_weighted_vs_mix",
  "log_lik_add",
  "log_lik_weighted_mix",
  "converged_weighted_fit_mix"
)

required_cs_columns <- c(
  "gene",
  "tissue",
  "model",
  "model_key",
  "lead_snp",
  "lead_coding",
  "distance_to_tss_kb"
)

missing_columns <- setdiff(required_columns, names(res))
missing_cs_columns <- setdiff(required_cs_columns, names(cs_res))
if (length(missing_columns) > 0L) {
  stop(
    "res_summary is missing weighted columns: ",
    paste(missing_columns, collapse = ", "),
    ". Rerun generate_summary_results.R."
  )
}
if (length(missing_cs_columns) > 0L) {
  stop(
    "res_cs_summary is missing weighted columns: ",
    paste(missing_cs_columns, collapse = ", "),
    ". Rerun generate_summary_results.R."
  )
}

safe_ratio_percent <- function(numerator, denominator) {
  if (length(denominator) == 0L || is.na(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  100 * numerator / denominator
}

safe_percent <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  100 * mean(x, na.rm = TRUE)
}

make_metric <- function(metric, count, denominator) {
  data.table(
    metric = metric,
    count = count,
    denominator = denominator,
    percentage = safe_ratio_percent(count, denominator)
  )
}

coding_levels <- c("Additive", "Recessive", "Dominant")
pattern_levels <- c("none", "A", "R", "D", "AR", "AD", "RD", "ADR")


# ============================================================
# Weighted analysis sets: same thresholds as descriptive results.R
# ============================================================

weighted_available <- res[
  has_weighted_fit_mix %in% TRUE & !is.na(ncs_weighted_fit_mix)
]

if (nrow(weighted_available) == 0L) {
  stop(
    paste0(
      "No weighted_fit_mix results were found. The RDS files must be created ",
      "by the weighted workhorse before generate_summary_results.R is run."
    )
  )
}

weighted_strong <- weighted_available[
  is.finite(min_pv) & min_pv < association_threshold
]

weighted_primary <- weighted_available[
  is.finite(min_pv) &
    min_pv < association_threshold &
    is.finite(mean_count) &
    mean_count >= minimum_mean_reads
]

primary_keys <- paste(weighted_primary$gene, weighted_primary$tissue, sep = "\r")
primary_cs <- cs_res[
  model_key %in% c("susie_add", "weighted_fit_mix") &
    paste(gene, tissue, sep = "\r") %in% primary_keys
]
weighted_cs <- primary_cs[model_key == "weighted_fit_mix"]

# A weighted fit with no CS must contribute zero counts, not be discarded.
coding_count_sum <- rowSums(
  weighted_primary[, .(n_add_weighted, n_rec_weighted, n_dom_weighted)]
)
if (any(coding_count_sum != weighted_primary$ncs_weighted_fit_mix)) {
  stop(
    paste0(
      "n_add_weighted + n_rec_weighted + n_dom_weighted does not equal ",
      "ncs_weighted_fit_mix. Regenerate both summaries from the same RDS files."
    )
  )
}

weighted_cs_counts <- weighted_cs[, .N, by = .(gene, tissue)]
expected_cs_counts <- weighted_primary[, .(
  gene,
  tissue,
  expected = ncs_weighted_fit_mix
)]
expected_cs_counts[
  weighted_cs_counts,
  on = .(gene, tissue),
  observed := i.N
]
expected_cs_counts[is.na(observed), observed := 0L]
if (any(expected_cs_counts$expected != expected_cs_counts$observed)) {
  stop(
    paste0(
      "Weighted rows in res_cs_summary do not match ncs_weighted_fit_mix. ",
      "Regenerate both summaries together."
    )
  )
}

lead_summary <- function(model_key_value, count_name, lead_name) {
  out <- primary_cs[
    model_key == model_key_value,
    .(
      cs_rows = .N,
      valid_leads = sum(!is.na(lead_snp) & nzchar(trimws(lead_snp))),
      lead_set = paste(sort(lead_snp), collapse = ";")
    ),
    by = .(gene, tissue)
  ]
  setnames(out, c("cs_rows", "valid_leads", "lead_set"),
           c(count_name, paste0(count_name, "_valid"), lead_name))
  out
}

agreement_data <- copy(weighted_primary)
for (definition in list(
  c("susie_add", "additive_cs_rows", "additive_lead_snp_set"),
  c("weighted_fit_mix", "weighted_cs_rows", "weighted_lead_snp_set")
)) {
  info <- lead_summary(definition[1], definition[2], definition[3])
  row_match <- match(
    paste(agreement_data$gene, agreement_data$tissue, sep = "\r"),
    paste(info$gene, info$tissue, sep = "\r")
  )
  count_values <- info[[definition[2]]][row_match]
  valid_values <- info[[paste0(definition[2], "_valid")]][row_match]
  lead_values <- info[[definition[3]]][row_match]
  count_values[is.na(count_values)] <- 0L
  valid_values[is.na(valid_values)] <- 0L
  lead_values[is.na(lead_values)] <- ""
  set(agreement_data, j = definition[2], value = as.integer(count_values))
  set(agreement_data, j = paste0(definition[2], "_valid"),
      value = as.integer(valid_values))
  set(agreement_data, j = definition[3], value = lead_values)
}

if (any(
  agreement_data$additive_cs_rows != agreement_data$ncs_susie |
  agreement_data$weighted_cs_rows != agreement_data$ncs_weighted_fit_mix |
  agreement_data$additive_cs_rows_valid != agreement_data$additive_cs_rows |
  agreement_data$weighted_cs_rows_valid != agreement_data$weighted_cs_rows
)) {
  stop(
    paste0(
      "Additive or weighted CS rows do not match the gene-tissue summary, ",
      "or a biological lead SNP is missing."
    )
  )
}

agreement_data[
  , additive_weighted_same_lead_snp_set :=
    ncs_susie > 0 &
    ncs_susie == ncs_weighted_fit_mix &
    additive_lead_snp_set == weighted_lead_snp_set
]


# ============================================================
# Overall weighted credible-set results
# ============================================================

describe_weighted_subset <- function(x, label) {
  data.table(
    subset = label,
    gene_tissue_pairs = nrow(x),
    genes = uniqueN(x$gene),
    tissues = uniqueN(x$tissue),
    total_weighted_cs = sum(x$ncs_weighted_fit_mix, na.rm = TRUE),
    pairs_with_weighted_cs = sum(x$ncs_weighted_fit_mix > 0, na.rm = TRUE),
    pct_with_weighted_cs = safe_percent(x$ncs_weighted_fit_mix > 0),
    converged_weighted_fits = sum(x$converged_weighted_fit_mix %in% TRUE),
    pct_converged_weighted_fits = safe_percent(x$converged_weighted_fit_mix)
  )
}

overall_summary <- rbindlist(list(
  describe_weighted_subset(weighted_available, "All weighted gene-tissue pairs"),
  describe_weighted_subset(
    weighted_strong,
    sprintf("P < %.1e", association_threshold)
  ),
  describe_weighted_subset(
    weighted_primary,
    sprintf(
      "P < %.1e and mean reads >= %g",
      association_threshold,
      minimum_mean_reads
    )
  )
))


# ============================================================
# Additive, recessive and dominant findings in the weighted fit
# ============================================================

coding_summary <- data.table(
  coding = coding_levels,
  credible_sets = c(
    sum(weighted_primary$n_add_weighted, na.rm = TRUE),
    sum(weighted_primary$n_rec_weighted, na.rm = TRUE),
    sum(weighted_primary$n_dom_weighted, na.rm = TRUE)
  )
)
coding_summary[
  , percentage := safe_ratio_percent(credible_sets, sum(credible_sets))
]

unweighted_coding_summary <- data.table(
  coding = coding_levels,
  credible_sets = c(
    sum(weighted_primary$n_add, na.rm = TRUE),
    sum(weighted_primary$n_rec, na.rm = TRUE),
    sum(weighted_primary$n_dom, na.rm = TRUE)
  )
)
unweighted_coding_summary[, model := "Unweighted mixed"]
coding_summary_for_comparison <- copy(coding_summary)
coding_summary_for_comparison[, model := "Weighted mixed"]

coding_comparison <- rbindlist(list(
  unweighted_coding_summary,
  coding_summary_for_comparison[, .(coding, credible_sets, model)]
))
setcolorder(coding_comparison, c("model", "coding", "credible_sets"))
coding_comparison[
  , percentage_within_model := safe_ratio_percent(credible_sets, sum(credible_sets)),
  by = model
]

coding_count_changes <- dcast(
  coding_comparison,
  coding ~ model,
  value.var = "credible_sets"
)
coding_count_changes[
  , `:=`(
    difference_weighted_minus_unweighted =
      `Weighted mixed` - `Unweighted mixed`,
    ratio_weighted_to_unweighted = fifelse(
      `Unweighted mixed` > 0,
      `Weighted mixed` / `Unweighted mixed`,
      NA_real_
    )
  )
]

pattern_data <- copy(weighted_primary)
pattern_data[
  , coding_pattern := paste0(
    fifelse(n_add_weighted > 0, "A", ""),
    fifelse(n_rec_weighted > 0, "R", ""),
    fifelse(n_dom_weighted > 0, "D", "")
  )
]
pattern_data[coding_pattern == "", coding_pattern := "none"]
pattern_data[coding_pattern == "ARD", coding_pattern := "ADR"]
pattern_data[
  , number_coding_types :=
    (n_add_weighted > 0) + (n_rec_weighted > 0) + (n_dom_weighted > 0)
]

coding_pattern_summary <- merge(
  data.table(coding_pattern = pattern_levels),
  pattern_data[, .(count = .N), by = coding_pattern],
  by = "coding_pattern",
  all.x = TRUE,
  sort = FALSE
)
coding_pattern_summary[is.na(count), count := 0L]
coding_pattern_summary[
  , percentage := safe_ratio_percent(count, nrow(pattern_data))
]
setorder(coding_pattern_summary, -count)

tissues <- sort(unique(pattern_data$tissue))
tissue_coding_patterns <- merge(
  CJ(tissue = tissues, coding_pattern = pattern_levels, unique = TRUE),
  pattern_data[, .(count = .N), by = .(tissue, coding_pattern)],
  by = c("tissue", "coding_pattern"),
  all.x = TRUE,
  sort = FALSE
)
tissue_coding_patterns[is.na(count), count := 0L]
tissue_coding_patterns[
  , `:=`(
    total_gene_tissue_pairs = sum(count),
    percentage = safe_ratio_percent(count, sum(count)),
    pattern_order = match(coding_pattern, pattern_levels)
  ),
  by = tissue
]
setorder(tissue_coding_patterns, tissue, pattern_order)

additive_regions <- pattern_data[n_add_weighted > 0][order(gene, tissue)]
recessive_regions <- pattern_data[n_rec_weighted > 0][order(gene, tissue)]
dominant_regions <- pattern_data[n_dom_weighted > 0][order(gene, tissue)]
adr_regions <- pattern_data[coding_pattern == "ADR"][order(gene, tissue)]


# ============================================================
# Tissue-specific weighted coding and additive enrichment
# ============================================================

recognized_codings <- c("additive", "dominant", "recessive")
tissue_additivity <- weighted_cs[
  , .(
    total_weighted_cs = .N,
    classified_weighted_cs = sum(lead_coding %in% recognized_codings),
    additive_cs = sum(lead_coding == "additive", na.rm = TRUE),
    dominant_cs = sum(lead_coding == "dominant", na.rm = TRUE),
    recessive_cs = sum(lead_coding == "recessive", na.rm = TRUE),
    unclassified_cs = sum(is.na(lead_coding) |
                            !lead_coding %in% recognized_codings)
  ),
  by = tissue
]
tissue_additivity[, nonadditive_cs := dominant_cs + recessive_cs]

total_additive <- sum(tissue_additivity$additive_cs)
total_nonadditive <- sum(tissue_additivity$nonadditive_cs)
total_classified <- total_additive + total_nonadditive
overall_additive_proportion <- if (total_classified > 0L) {
  total_additive / total_classified
} else {
  NA_real_
}

tissue_additivity[
  , additive_proportion := additive_cs / classified_weighted_cs
]
z_value <- qnorm(0.975)
tissue_additivity[
  , `:=`(
    additive_ci_lower = (
      additive_proportion + z_value^2 / (2 * classified_weighted_cs) -
        z_value * sqrt(
          additive_proportion * (1 - additive_proportion) /
            classified_weighted_cs +
            z_value^2 / (4 * classified_weighted_cs^2)
        )
    ) / (1 + z_value^2 / classified_weighted_cs),
    additive_ci_upper = (
      additive_proportion + z_value^2 / (2 * classified_weighted_cs) +
        z_value * sqrt(
          additive_proportion * (1 - additive_proportion) /
            classified_weighted_cs +
            z_value^2 / (4 * classified_weighted_cs^2)
        )
    ) / (1 + z_value^2 / classified_weighted_cs)
  )
]
tissue_additivity[!is.finite(additive_ci_lower), additive_ci_lower := NA_real_]
tissue_additivity[!is.finite(additive_ci_upper), additive_ci_upper := NA_real_]

if (nrow(tissue_additivity) > 0L) {
  fisher_results <- rbindlist(lapply(seq_len(nrow(tissue_additivity)), function(i) {
    contingency <- matrix(c(
      tissue_additivity$additive_cs[i],
      tissue_additivity$nonadditive_cs[i],
      total_additive - tissue_additivity$additive_cs[i],
      total_nonadditive - tissue_additivity$nonadditive_cs[i]
    ), nrow = 2L, byrow = TRUE)
    test <- tryCatch(fisher.test(contingency), error = function(e) NULL)
    if (is.null(test)) {
      return(data.table(additive_odds_ratio = NA_real_, fisher_p_value = NA_real_))
    }
    estimate <- as.numeric(test$estimate)
    if (length(estimate) == 0L) estimate <- NA_real_
    data.table(additive_odds_ratio = estimate[1], fisher_p_value = test$p.value)
  }))
  tissue_additivity <- cbind(tissue_additivity, fisher_results)
  tissue_additivity[
    , `:=`(
      additive_percentage = 100 * additive_proportion,
      additive_ci_lower_percentage = 100 * additive_ci_lower,
      additive_ci_upper_percentage = 100 * additive_ci_upper,
      overall_additive_percentage = 100 * overall_additive_proportion,
      additive_difference_percentage_points =
        100 * (additive_proportion - overall_additive_proportion),
      fisher_fdr = p.adjust(fisher_p_value, method = "BH")
    )
  ]
  setorder(tissue_additivity, -additive_difference_percentage_points)
}


# ============================================================
# Additive versus weighted mixed-model agreement
# ============================================================

n_additive <- agreement_data$ncs_susie
n_weighted <- weighted_primary$ncs_weighted_fit_mix
shared_snps <- agreement_data$overlap_snp_add_weighted
same_lead <- agreement_data$additive_weighted_same_lead_snp_set

agreement_levels <- c(
  "Same CS count; same lead SNP(s)",
  "Same CS count; different lead SNP(s); SNP overlap",
  "Different CS count",
  "Same CS count; no SNP overlap",
  "Neither model reported a CS",
  "Additive model only",
  "Weighted mixed model only"
)

agreement_category <- fcase(
  n_additive == 0 & n_weighted == 0, agreement_levels[5],
  n_additive > 0 & n_weighted == 0, agreement_levels[6],
  n_additive == 0 & n_weighted > 0, agreement_levels[7],
  n_additive == n_weighted & shared_snps > 0 & same_lead,
  agreement_levels[1],
  n_additive == n_weighted & shared_snps > 0 & !same_lead,
  agreement_levels[2],
  n_additive == n_weighted & shared_snps == 0, agreement_levels[4],
  default = agreement_levels[3]
)

agreement_summary <- data.table(
  agreement_category = agreement_levels,
  count = vapply(
    agreement_levels,
    function(level) sum(agreement_category == level),
    integer(1L)
  )
)
agreement_summary[
  , percentage := safe_ratio_percent(count, sum(count))
]

both_models <- n_additive > 0 & n_weighted > 0
agreement_statistics <- rbindlist(list(
  make_metric(
    "Same number of credible sets",
    sum(n_additive == n_weighted),
    nrow(weighted_primary)
  ),
  make_metric(
    "Same CS count and same biological lead SNP set",
    sum(n_additive == n_weighted & n_weighted > 0 & same_lead),
    nrow(weighted_primary)
  ),
  make_metric(
    "Different number of credible sets",
    sum(n_additive != n_weighted),
    nrow(weighted_primary)
  ),
  make_metric(
    "Both models reported at least one CS",
    sum(both_models),
    nrow(weighted_primary)
  ),
  make_metric(
    "Any shared SNP when both models reported CSs",
    sum(both_models & shared_snps > 0),
    sum(both_models)
  )
))

cs_count_table <- as.data.table(with(
  weighted_primary,
  table(
    additive_cs = ncs_susie,
    weighted_mixed_cs = ncs_weighted_fit_mix
  )
))
setnames(cs_count_table, "N", "count")
cs_count_table[
  , percentage := safe_ratio_percent(count, sum(count))
]

fit_comparison_summary <- data.table(
  metric = c(
    "Median weighted-minus-additive ELBO",
    "Mean weighted-minus-additive ELBO",
    "Median weighted-minus-unweighted ELBO",
    "Mean weighted-minus-unweighted ELBO",
    "Median additive log-likelihood metric",
    "Median weighted log-likelihood metric"
  ),
  value = c(
    median(weighted_primary$dif_elbo_weighted_vs_add, na.rm = TRUE),
    mean(weighted_primary$dif_elbo_weighted_vs_add, na.rm = TRUE),
    median(weighted_primary$dif_elbo_weighted_vs_mix, na.rm = TRUE),
    mean(weighted_primary$dif_elbo_weighted_vs_mix, na.rm = TRUE),
    median(weighted_primary$log_lik_add, na.rm = TRUE),
    median(weighted_primary$log_lik_weighted_mix, na.rm = TRUE)
  )
)

# This smaller table addresses the prior-sensitivity question directly,
# without replacing the primary additive-versus-weighted replication.
prior_sensitivity_statistics <- rbindlist(list(
  make_metric(
    "Same number of credible sets",
    sum(weighted_primary$ncs_susie_mix == n_weighted),
    nrow(weighted_primary)
  ),
  make_metric(
    "Same biological-SNP credible sets",
    sum(weighted_primary$mix_weighted_same_cs_snp_sets %in% TRUE),
    nrow(weighted_primary)
  ),
  make_metric(
    "Same biological lead-SNP set",
    sum(weighted_primary$mix_weighted_same_lead_snp_set %in% TRUE),
    nrow(weighted_primary)
  ),
  make_metric(
    "Any shared SNP when both models reported CSs",
    sum(
      weighted_primary$ncs_susie_mix > 0 &
        n_weighted > 0 &
        weighted_primary$mix_weighted_overlap_snp > 0
    ),
    sum(weighted_primary$ncs_susie_mix > 0 & n_weighted > 0)
  )
))


# ============================================================
# Weighted lead-SNP distance to the TSS
# ============================================================

bin_centers <- seq(
  -tss_plot_limit_kb,
  tss_plot_limit_kb,
  by = tss_bin_width_kb
)
bin_breaks <- c(
  bin_centers - tss_bin_width_kb / 2,
  tail(bin_centers, 1) + tss_bin_width_kb / 2
)

make_tss_summary <- function(
    x,
    model_name,
    tissue_name = "All tissues") {
  distance <- x[
    is.finite(distance_to_tss_kb),
    distance_to_tss_kb
  ]
  plotted <- distance[
    distance >= min(bin_breaks) & distance <= max(bin_breaks)
  ]
  counts <- hist(
    plotted,
    breaks = bin_breaks,
    plot = FALSE,
    include.lowest = TRUE,
    right = FALSE
  )$counts
  data.table(
    tissue = tissue_name,
    model = model_name,
    distance_to_tss_kb = bin_centers,
    bin_lower_kb = head(bin_breaks, -1),
    bin_upper_kb = tail(bin_breaks, -1),
    count = counts,
    proportion = if (length(plotted) > 0L) counts / length(plotted) else NA_real_,
    percentage = if (length(plotted) > 0L) 100 * counts / length(plotted) else NA_real_,
    n_cs_total = length(distance),
    n_cs_in_window = length(plotted)
  )
}

tss_distance_summary <- rbindlist(list(
  make_tss_summary(
    primary_cs[model_key == "susie_add"],
    "SuSiE"
  ),
  make_tss_summary(
    weighted_cs,
    "Weighted SuSiE-mix"
  )
))
tissue_tss_distance_summary <- rbindlist(lapply(tissues, function(tissue_name) {
  rbindlist(list(
    make_tss_summary(
      primary_cs[tissue == tissue_name & model_key == "susie_add"],
      "SuSiE",
      tissue_name
    ),
    make_tss_summary(
      weighted_cs[tissue == tissue_name],
      "Weighted SuSiE-mix",
      tissue_name
    )
  ))
}))


# ============================================================
# Save weighted-only tables and figures
# ============================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tables <- list(
  weighted_overall_summary = overall_summary,
  weighted_mixed_coding_summary = coding_summary,
  weighted_vs_unweighted_coding_summary = coding_comparison,
  weighted_vs_unweighted_coding_count_changes = coding_count_changes,
  weighted_coding_patterns = coding_pattern_summary,
  weighted_tissue_coding_patterns = tissue_coding_patterns[
    , setdiff(names(tissue_coding_patterns), "pattern_order"), with = FALSE
  ],
  weighted_tissue_additivity_summary = tissue_additivity,
  weighted_vs_additive_agreement_categories = agreement_summary,
  weighted_vs_additive_agreement_statistics = agreement_statistics,
  weighted_vs_additive_cs_count_table = cs_count_table,
  weighted_vs_unweighted_prior_sensitivity = prior_sensitivity_statistics,
  weighted_fit_comparison_summary = fit_comparison_summary,
  weighted_tss_distance_distribution = tss_distance_summary,
  weighted_tissue_tss_distance_distribution = tissue_tss_distance_summary,
  weighted_additive_gene_tissue_pairs = additive_regions,
  weighted_recessive_gene_tissue_pairs = recessive_regions,
  weighted_dominant_gene_tissue_pairs = dominant_regions,
  weighted_adr_gene_tissue_pairs = adr_regions,
  weighted_primary_analysis_set = pattern_data
)

for (table_name in names(tables)) {
  fwrite(
    tables[[table_name]],
    file.path(output_dir, paste0(table_name, ".csv"))
  )
}

pdf(
  file.path(output_dir, "weighted_coding_patterns.pdf"),
  width = 8,
  height = 5
)
plot_order <- order(-coding_pattern_summary$count)
plot_data <- coding_pattern_summary[plot_order]
heights <- fifelse(is.finite(plot_data$percentage), plot_data$percentage, 0)
positions <- barplot(
  heights,
  names.arg = plot_data$coding_pattern,
  col = "#009E73",
  border = NA,
  xlab = "Coding pattern",
  ylab = "Percentage of gene-tissue pairs",
  main = "Weighted mixed-model coding patterns"
)
text(positions, heights, labels = plot_data$count, pos = 3, cex = 0.8)
invisible(dev.off())

pdf(
  file.path(output_dir, "weighted_lead_snp_tss_distribution.pdf"),
  width = 8,
  height = 5
)
valid_tss <- is.finite(tss_distance_summary$proportion)
y_max <- if (any(valid_tss)) {
  max(0.01, 1.1 * max(tss_distance_summary$proportion[valid_tss]))
} else {
  1
}
plot(
  NA_real_,
  NA_real_,
  xlim = c(-tss_plot_limit_kb, tss_plot_limit_kb),
  ylim = c(0, y_max),
  xlab = "Distance to TSS (kb)",
  ylab = "Proportion of credible sets",
  main = "Additive versus weighted mixed CS leads"
)
abline(v = 0, col = "gray65", lty = 3)
model_colors <- c(
  SuSiE = "#E69F00",
  `Weighted SuSiE-mix` = "#009E73"
)
for (model_name in names(model_colors)) {
  temp <- tss_distance_summary[model == model_name]
  lines(
    temp$distance_to_tss_kb,
    temp$proportion,
    lwd = 2.5,
    col = model_colors[[model_name]]
  )
}
legend(
  "topright",
  legend = names(model_colors),
  col = unname(model_colors),
  lwd = 2.5,
  bty = "n"
)
invisible(dev.off())


# ============================================================
# Console summary
# ============================================================

cat("\nWeighted descriptive analysis complete.\n")
cat("Primary gene-tissue analyses:", nrow(weighted_primary), "\n")
cat("Weighted credible sets:", sum(coding_summary$credible_sets), "\n")
cat("Additive:", coding_summary[coding == "Additive", credible_sets], "\n")
cat("Recessive:", coding_summary[coding == "Recessive", credible_sets], "\n")
cat("Dominant:", coding_summary[coding == "Dominant", credible_sets], "\n")
cat(
  "Same CS count as additive model:",
  sprintf("%.1f%%", safe_percent(n_additive == n_weighted)),
  "\n"
)
cat("Outputs saved in:", output_dir, "\n")
