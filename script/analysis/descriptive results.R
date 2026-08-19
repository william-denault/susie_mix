library(data.table)

load(
  "/project2/mstephens/wdenault/susie_mix/res_summary.RData"
)

res <- as.data.table(res_summary)

# ============================================================
# Analysis thresholds
# ============================================================

association_threshold <- 1e-6
minimum_mean_reads <- 100

strong_label <- sprintf(
  "P < %.1e",
  association_threshold
)

primary_label <- sprintf(
  "P < %.1e and mean reads >= %g",
  association_threshold,
  minimum_mean_reads
)

required_columns <- c(
  "gene",
  "tissue",
  "min_pv",
  "mean_count",
  "ncs_susie",
  "ncs_susie_mix",
  "overlap_snp",
  "n_add_cs_snps",
  "n_mix_cs_snps",
  "overlap_coding",
  "n_add_cs_snps_not_retained",
  "n_add_cs_snps_other_coding_only",
  "n_shared_snps_preferred_additive",
  "n_shared_snps_preferred_dominant",
  "n_shared_snps_preferred_recessive",
  "n_shared_snps_preferred_nonadditive",
  "n_shared_snps_preference_ambiguous",
  "perm_cs_susie",
  "perm_cs_susie_mix",
  "dif_elbo",
  "dif_elbo_perm",
  "log_lik_add",
  "log_lik_mix",
  "n_add",
  "n_rec",
  "n_dom",
  "n_add_perm",
  "n_rec_perm",
  "n_dom_perm"
)

missing_columns <- setdiff(
  required_columns,
  names(res)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", "),
    ". Rerun summarize_susie_results.R with the new overlap calculations."
  )
}


# ============================================================
# Helper functions
# ============================================================

safe_percent <- function(x) {

  if (
    length(x) == 0L ||
    all(is.na(x))
  ) {
    return(NA_real_)
  }

  100 * mean(
    x,
    na.rm = TRUE
  )
}


safe_ratio_percent <- function(
    numerator,
    denominator) {

  percentage <- 100 * numerator / denominator

  invalid_denominator <- (
    is.na(denominator) |
      denominator <= 0
  )

  percentage[invalid_denominator] <- NA_real_

  percentage
}


# ============================================================
# Define analysis sets
# ============================================================

res_strong <- res[
  is.finite(min_pv) &
    min_pv < association_threshold
]

res_idx <- res[
  is.finite(min_pv) &
    min_pv < association_threshold &
    is.finite(mean_count) &
    mean_count >= minimum_mean_reads &
    !is.na(ncs_susie) &
    !is.na(ncs_susie_mix) &
    !is.na(overlap_snp)
]

cat("All gene-tissue pairs:", nrow(res), "\n")
cat("Strong association pairs:", nrow(res_strong), "\n")
cat(
  "Strong, sufficiently expressed pairs:",
  nrow(res_idx),
  "\n"
)
cat(
  "Distinct genes in primary set:",
  uniqueN(res_idx$gene),
  "\n"
)
cat(
  "Distinct tissues in primary set:",
  uniqueN(res_idx$tissue),
  "\n\n"
)


# ============================================================
# General credible-set summaries
# ============================================================

describe_subset <- function(x, label) {

  data.table(
    subset = label,
    gene_tissue_pairs = nrow(x),
    genes = uniqueN(x$gene),
    tissues = uniqueN(x$tissue),

    total_additive_cs = sum(
      x$ncs_susie,
      na.rm = TRUE
    ),

    total_mixed_cs = sum(
      x$ncs_susie_mix,
      na.rm = TRUE
    ),

    pct_with_additive_cs = safe_percent(
      x$ncs_susie > 0
    ),

    pct_with_mixed_cs = safe_percent(
      x$ncs_susie_mix > 0
    ),

    total_permuted_additive_cs = sum(
      x$perm_cs_susie,
      na.rm = TRUE
    ),

    total_permuted_mixed_cs = sum(
      x$perm_cs_susie_mix,
      na.rm = TRUE
    ),

    pct_with_permuted_additive_cs = safe_percent(
      x$perm_cs_susie > 0
    ),

    pct_with_permuted_mixed_cs = safe_percent(
      x$perm_cs_susie_mix > 0
    )
  )
}


overall_summary <- rbindlist(
  list(
    describe_subset(
      res,
      "All gene-tissue pairs"
    ),
    describe_subset(
      res_strong,
      strong_label
    ),
    describe_subset(
      res_idx,
      primary_label
    )
  )
)

print(overall_summary)


# ============================================================
# SNP overlap and coding-aware overlap
# ============================================================

summarize_overlap_and_coding <- function(
    x,
    label) {

  total_add_cs_snps <- sum(
    x$n_add_cs_snps,
    na.rm = TRUE
  )

  total_mix_cs_snps <- sum(
    x$n_mix_cs_snps,
    na.rm = TRUE
  )

  total_shared_snps <- sum(
    x$overlap_snp,
    na.rm = TRUE
  )

  total_not_retained <- sum(
    x$n_add_cs_snps_not_retained,
    na.rm = TRUE
  )

  total_additive_coding_overlap <- sum(
    x$overlap_coding,
    na.rm = TRUE
  )

  total_other_coding_only <- sum(
    x$n_add_cs_snps_other_coding_only,
    na.rm = TRUE
  )

  preferred_additive <- sum(
    x$n_shared_snps_preferred_additive,
    na.rm = TRUE
  )

  preferred_dominant <- sum(
    x$n_shared_snps_preferred_dominant,
    na.rm = TRUE
  )

  preferred_recessive <- sum(
    x$n_shared_snps_preferred_recessive,
    na.rm = TRUE
  )

  preferred_nonadditive <- sum(
    x$n_shared_snps_preferred_nonadditive,
    na.rm = TRUE
  )

  preference_ambiguous <- sum(
    x$n_shared_snps_preference_ambiguous,
    na.rm = TRUE
  )

  data.table(
    subset = label,
    gene_tissue_pairs = nrow(x),

    additive_cs_snp_occurrences = total_add_cs_snps,
    mixed_cs_snp_occurrences = total_mix_cs_snps,

    shared_snp_occurrences = total_shared_snps,

    pct_additive_cs_snps_retained_any_coding = (
      safe_ratio_percent(
        total_shared_snps,
        total_add_cs_snps
      )
    ),

    additive_cs_snps_not_retained = total_not_retained,

    pct_additive_cs_snps_not_retained = (
      safe_ratio_percent(
        total_not_retained,
        total_add_cs_snps
      )
    ),

    additive_coding_overlap = (
      total_additive_coding_overlap
    ),

    pct_additive_cs_snps_retained_as_additive = (
      safe_ratio_percent(
        total_additive_coding_overlap,
        total_add_cs_snps
      )
    ),

    additive_cs_snps_retained_other_coding_only = (
      total_other_coding_only
    ),

    pct_additive_cs_snps_retained_other_coding_only = (
      safe_ratio_percent(
        total_other_coding_only,
        total_add_cs_snps
      )
    ),

    pct_shared_snps_other_coding_only = (
      safe_ratio_percent(
        total_other_coding_only,
        total_shared_snps
      )
    ),

    shared_snps_preferred_additive = preferred_additive,
    shared_snps_preferred_dominant = preferred_dominant,
    shared_snps_preferred_recessive = preferred_recessive,
    shared_snps_preferred_nonadditive = preferred_nonadditive,
    shared_snps_preference_ambiguous = preference_ambiguous,

    pct_shared_snps_preferred_nonadditive = (
      safe_ratio_percent(
        preferred_nonadditive,
        total_shared_snps
      )
    )
  )
}


overlap_coding_summary <- rbindlist(
  list(
    summarize_overlap_and_coding(
      res,
      "All gene-tissue pairs"
    ),
    summarize_overlap_and_coding(
      res_strong,
      strong_label
    ),
    summarize_overlap_and_coding(
      res_idx,
      primary_label
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

print(overlap_coding_summary)


coding_preference_summary <- data.table(
  preferred_coding = c(
    "Additive",
    "Dominant",
    "Recessive",
    "Ambiguous or unresolved"
  ),
  shared_snp_occurrences = c(
    sum(
      res_idx$n_shared_snps_preferred_additive,
      na.rm = TRUE
    ),
    sum(
      res_idx$n_shared_snps_preferred_dominant,
      na.rm = TRUE
    ),
    sum(
      res_idx$n_shared_snps_preferred_recessive,
      na.rm = TRUE
    ),
    sum(
      res_idx$n_shared_snps_preference_ambiguous,
      na.rm = TRUE
    )
  )
)

coding_preference_summary[
  ,
  percentage := safe_ratio_percent(
    shared_snp_occurrences,
    sum(shared_snp_occurrences)
  )
]

print(coding_preference_summary)


# ============================================================
# Define agreement categories
# ============================================================

agreement_levels <- c(
  "Same CS count; SNP overlap",
  "Different CS count",
  "Same CS count; no SNP overlap",
  "Neither model reported a CS",
  "Additive model only",
  "Mixed model only"
)

res_idx[
  ,
  agreement_category := fcase(
    ncs_susie == 0 &
      ncs_susie_mix == 0,
    "Neither model reported a CS",

    ncs_susie > 0 &
      ncs_susie_mix == 0,
    "Additive model only",

    ncs_susie == 0 &
      ncs_susie_mix > 0,
    "Mixed model only",

    ncs_susie == ncs_susie_mix &
      overlap_snp > 0,
    "Same CS count; SNP overlap",

    ncs_susie == ncs_susie_mix &
      overlap_snp == 0,
    "Same CS count; no SNP overlap",

    default = "Different CS count"
  )
]

res_idx[
  ,
  agreement_category := factor(
    agreement_category,
    levels = agreement_levels
  )
]

agreement_summary <- data.table(
  agreement_category = agreement_levels,
  count = vapply(
    agreement_levels,
    function(category) {
      sum(
        as.character(res_idx$agreement_category) == category,
        na.rm = TRUE
      )
    },
    integer(1L)
  )
)

agreement_summary[
  ,
  percentage := safe_ratio_percent(
    count,
    sum(count)
  )
]

print(agreement_summary)


# ============================================================
# Main agreement statistics
# ============================================================

both_models_cs <- res_idx[
  ncs_susie > 0 &
    ncs_susie_mix > 0
]

one_cs_each <- res_idx[
  ncs_susie == 1 &
    ncs_susie_mix == 1
]


make_metric <- function(
    metric,
    count,
    denominator) {

  data.table(
    metric = metric,
    count = count,
    denominator = denominator,
    percentage = safe_ratio_percent(
      count,
      denominator
    )
  )
}


agreement_statistics <- rbindlist(
  list(
    make_metric(
      "Same number of credible sets",
      sum(
        res_idx$ncs_susie ==
          res_idx$ncs_susie_mix
      ),
      nrow(res_idx)
    ),

    make_metric(
      "Different number of credible sets",
      sum(
        res_idx$ncs_susie !=
          res_idx$ncs_susie_mix
      ),
      nrow(res_idx)
    ),

    make_metric(
      "Both models reported at least one CS",
      nrow(both_models_cs),
      nrow(res_idx)
    ),

    make_metric(
      "Any shared SNP when both models reported CSs",
      sum(
        both_models_cs$overlap_snp > 0
      ),
      nrow(both_models_cs)
    ),

    make_metric(
      "No shared SNP when both models reported CSs",
      sum(
        both_models_cs$overlap_snp == 0
      ),
      nrow(both_models_cs)
    ),

    make_metric(
      "Exactly one CS from each model",
      nrow(one_cs_each),
      nrow(res_idx)
    ),

    make_metric(
      "One CS each with at least one shared SNP",
      sum(
        one_cs_each$overlap_snp > 0
      ),
      nrow(one_cs_each)
    ),

    make_metric(
      "One CS each with no shared SNP",
      sum(
        one_cs_each$overlap_snp == 0
      ),
      nrow(one_cs_each)
    )
  )
)

print(agreement_statistics)


# ============================================================
# Cross-tabulation of credible-set counts
# ============================================================

cs_count_table <- as.data.table(
  with(
    res_idx,
    table(
      additive_cs = ncs_susie,
      mixed_cs = ncs_susie_mix
    )
  )
)

setnames(
  cs_count_table,
  "N",
  "count"
)

cs_count_table[
  ,
  percentage := safe_ratio_percent(
    count,
    sum(count)
  )
]

print(cs_count_table)


# ============================================================
# Coding composition of mixed-model credible sets
# ============================================================

coding_summary <- data.table(
  coding = c(
    "Additive",
    "Recessive",
    "Dominant"
  ),
  credible_sets = c(
    sum(res_idx$n_add, na.rm = TRUE),
    sum(res_idx$n_rec, na.rm = TRUE),
    sum(res_idx$n_dom, na.rm = TRUE)
  )
)

coding_summary[
  ,
  percentage := safe_ratio_percent(
    credible_sets,
    sum(credible_sets)
  )
]

print(coding_summary)


res_idx[
  ,
  coding_pattern := paste0(
    ifelse(!is.na(n_add) & n_add > 0, "A", ""),
    ifelse(!is.na(n_rec) & n_rec > 0, "R", ""),
    ifelse(!is.na(n_dom) & n_dom > 0, "D", "")
  )
]

res_idx[
  coding_pattern == "",
  coding_pattern := "none"
]

# Use ADR as the canonical label for analyses containing all three
# coding types, while retaining the existing two-type labels.
res_idx[
  coding_pattern == "ARD",
  coding_pattern := "ADR"
]

res_idx[
  ,
  number_coding_types := (
    (!is.na(n_add) & n_add > 0) +
      (!is.na(n_rec) & n_rec > 0) +
      (!is.na(n_dom) & n_dom > 0)
  )
]

coding_pattern_levels <- c(
  "none",
  "A",
  "R",
  "D",
  "AR",
  "AD",
  "RD",
  "ADR"
)

coding_pattern_summary <- merge(
  data.table(
    coding_pattern = coding_pattern_levels
  ),
  res_idx[
    ,
    .(count = .N),
    by = coding_pattern
  ],
  by = "coding_pattern",
  all.x = TRUE,
  sort = FALSE
)[
  is.na(count),
  count := 0L
][
  ,
  percentage := safe_ratio_percent(
    count,
    nrow(res_idx)
  )
][
  order(-count)
]

print(coding_pattern_summary)

adr_regions <- res_idx[
  coding_pattern == "ADR",
  .(
    gene,
    tissue,
    min_pv,
    mean_count,
    ncs_susie_mix,
    n_add,
    n_dom,
    n_rec
  )
][
  order(gene, tissue)
]

cat(
  "Gene-tissue pairs with additive, dominant, and recessive coding (ADR):",
  nrow(adr_regions),
  "\n"
)

cat(
  "Percentage with ADR coding:",
  safe_ratio_percent(
    nrow(adr_regions),
    nrow(res_idx)
  ),
  "\n"
)

cat(
  "Gene-tissue pairs with at least two coding types:",
  sum(
    res_idx$number_coding_types >= 2
  ),
  "\n"
)

cat(
  "Percentage with at least two coding types:",
  safe_percent(
    res_idx$number_coding_types >= 2
  ),
  "\n"
)


# ============================================================
# Permutation control
# ============================================================

permutation_summary <- data.table(
  model = c(
    "Additive",
    "Mixed"
  ),
  total_credible_sets = c(
    sum(
      res_idx$perm_cs_susie,
      na.rm = TRUE
    ),
    sum(
      res_idx$perm_cs_susie_mix,
      na.rm = TRUE
    )
  ),
  analyses_with_at_least_one_cs = c(
    sum(
      res_idx$perm_cs_susie > 0,
      na.rm = TRUE
    ),
    sum(
      res_idx$perm_cs_susie_mix > 0,
      na.rm = TRUE
    )
  )
)

permutation_summary[
  ,
  percentage_with_at_least_one_cs := (
    safe_ratio_percent(
      analyses_with_at_least_one_cs,
      nrow(res_idx)
    )
  )
]

print(permutation_summary)


# ============================================================
# Discordant one-CS cases
# ============================================================

res_1cs <- copy(
  res_idx[
    ncs_susie == 1 &
      ncs_susie_mix == 1 &
      overlap_snp == 0
  ]
)

res_1cs[
  ,
  likelihood_statistic := (
    -2 *
      (
        log_lik_add -
          log_lik_mix
      )
  )
]

setorderv(
  res_1cs,
  cols = c(
    "dif_elbo",
    "min_pv"
  ),
  order = c(
    -1L,
    1L
  ),
  na.last = TRUE
)

candidate_columns <- intersect(
  c(
    "gene",
    "tissue",
    "min_pv",
    "mean_count",
    "n_ind",
    "n_SNP",
    "n_mix_predictor",
    "dif_elbo",
    "likelihood_statistic",
    "n_add",
    "n_rec",
    "n_dom",
    "n_add_cs_snps",
    "n_mix_cs_snps",
    "overlap_snp",
    "overlap_coding",
    "n_add_cs_snps_other_coding_only",
    "n_shared_snps_preferred_nonadditive"
  ),
  names(res_1cs)
)

top_discordant_candidates <- head(
  res_1cs[
    ,
    ..candidate_columns
  ],
  50
)

print(top_discordant_candidates)

dominant_candidates <- res_1cs[
  n_dom > 0
]

recessive_candidates <- res_1cs[
  n_rec > 0
]

additive_candidates <- res_1cs[
  n_add > 0
]


# ============================================================
# Tissue-specific disagreement and coding preference
# ============================================================

res_idx[
  ,
  different_cs_number := (
    ncs_susie != ncs_susie_mix
  )
]

res_idx[
  ,
  no_snp_overlap_when_both_detect := (
    ncs_susie > 0 &
      ncs_susie_mix > 0 &
      overlap_snp == 0
  )
]

tissue_summary <- res_idx[
  ,
  .(
    total = .N,

    different_cs_number = sum(
      different_cs_number
    ),

    no_snp_overlap = sum(
      no_snp_overlap_when_both_detect
    ),

    shared_snp_occurrences = sum(
      overlap_snp,
      na.rm = TRUE
    ),

    shared_snps_preferred_nonadditive = sum(
      n_shared_snps_preferred_nonadditive,
      na.rm = TRUE
    )
  ),
  by = tissue
]

tissue_summary[
  ,
  pct_different_cs_number := (
    safe_ratio_percent(
      different_cs_number,
      total
    )
  )
]

tissue_summary[
  ,
  pct_no_snp_overlap := (
    safe_ratio_percent(
      no_snp_overlap,
      total
    )
  )
]

tissue_summary[
  ,
  pct_shared_snps_preferred_nonadditive := (
    safe_ratio_percent(
      shared_snps_preferred_nonadditive,
      shared_snp_occurrences
    )
  )
]

setorder(
  tissue_summary,
  -pct_no_snp_overlap
)

print(tissue_summary)


# ============================================================
# Descriptive plots
# ============================================================

plot_coding_pattern_summary <- function(x) {

  plot_data <- copy(x)
  setorder(plot_data, -count)

  valid_percentage <- is.finite(
    plot_data$percentage
  )

  if (!any(valid_percentage)) {
    plot.new()
    title("Mixed-model coding patterns")
    text(
      0.5,
      0.5,
      "No gene-tissue pairs in the analysis set"
    )
    return(invisible(NULL))
  }

  bar_heights <- ifelse(
    valid_percentage,
    plot_data$percentage,
    0
  )

  y_max <- max(
    1,
    1.15 * max(bar_heights)
  )

  pattern_colors <- c(
    A = "#0072B2",
    R = "#009E73",
    D = "#D55E00",
    AR = "#56B4E9",
    AD = "#CC79A7",
    RD = "#E69F00",
    ADR = "#6A3D9A",
    none = "#999999"
  )

  bar_colors <- unname(
    pattern_colors[
      plot_data$coding_pattern
    ]
  )

  par(mar = c(5, 5, 4, 1))

  bar_positions <- barplot(
    bar_heights,
    names.arg = plot_data$coding_pattern,
    ylim = c(0, y_max),
    xlab = "Coding pattern",
    ylab = "Percentage of gene-tissue pairs",
    main = "Mixed-model coding patterns",
    col = bar_colors,
    border = NA
  )

  text(
    x = bar_positions,
    y = bar_heights,
    labels = plot_data$count,
    pos = 3,
    cex = 0.8
  )

  mtext(
    "A = additive; R = recessive; D = dominant; counts above bars",
    side = 3,
    line = 0.25,
    adj = 1,
    cex = 0.7
  )

  invisible(bar_positions)
}


par(
  mfrow = c(2, 2),
  mar = c(5, 5, 4, 1)
)

likelihood_values <- -2 * (
  res_idx$log_lik_add -
    res_idx$log_lik_mix
)

hist(
  likelihood_values,
  xlim = c(-20, 50),
  nclass = 1000,
  main = "-2(log lik additive - log lik mixed)",
  xlab = "-2(log lik additive - log lik mixed)",
  border = "gray80"
)

abline(
  v = 0,
  col = "red",
  lty = 2,
  lwd = 2
)

permuted_elbo <- res_idx$dif_elbo_perm[
  is.finite(res_idx$dif_elbo_perm)
]

hist(
  permuted_elbo,
  xlim = c(-20, 50),
  nclass = 1000,
  main = "ELBO difference, permuted phenotype",
  xlab = "ELBO mixed - ELBO additive",
  border = "gray80"
)

abline(
  v = 0,
  col = "red",
  lty = 2,
  lwd = 2
)

plot(
  -log10(res_idx$min_pv),
  pch = 16,
  cex = 0.5,
  col = "steelblue",
  xlab = "Gene-tissue analysis",
  ylab = expression(-log[10](P)),
  main = "Marginal association strength"
)

abline(
  h = -log10(association_threshold),
  col = "red",
  lty = 2,
  lwd = 2
)

agreement_plot_labels <- c(
  "Same count;\nSNP overlap",
  "Different\nCS count",
  "Same count;\nno SNP overlap",
  "Neither model\nreported a CS",
  "Additive\nmodel only",
  "Mixed\nmodel only"
)

par(mar = c(8, 5, 4, 1))

barplot(
  agreement_summary$percentage,
  names.arg = agreement_plot_labels,
  cex.names = 0.65,
  ylab = "Percentage",
  main = "Agreement between fine-mapping models",
  col = "steelblue"
)

par(mfrow = c(1, 1))


par(
  mfrow = c(2, 2),
  mar = c(5, 5, 4, 1)
)




agreement_plot_labels <- c(
  "Same count;\nSNP overlap",
  "Different\nCS count",
  "Same count;\nno SNP overlap",
  "Neither model\nreported a CS",
  "Additive\nmodel only",
  "Mixed\nmodel only"
)

par(mar = c(8, 5, 4, 1))

barplot(
  agreement_summary$percentage,
  names.arg = agreement_plot_labels,
  cex.names = 0.65,
  ylab = "Percentage",
  main = "Agreement between fine-mapping models",
  col = "steelblue"
)
if (
  sum(
    coding_preference_summary$shared_snp_occurrences
  ) > 0
) {

  barplot(
    coding_preference_summary$percentage,
    names.arg = coding_preference_summary$preferred_coding,
    las = 2,
    cex.names = 0.8,
    ylab = "Percentage of shared SNP occurrences",
    main = "Mixed-model preferred coding",
    col = c(
      "steelblue",
      "darkorange",
      "firebrick",
      "gray70"
    )
  )
}

plot_coding_pattern_summary(
  coding_pattern_summary
)


par(mfrow = c(1, 1))

# ============================================================
# Save descriptive results
# ============================================================

output_dir <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "descriptive_results/"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pdf(
  file.path(
    output_dir,
    "coding_patterns_barplot.pdf"
  ),
  width = 8,
  height = 5
)

plot_coding_pattern_summary(
  coding_pattern_summary
)

invisible(dev.off())

fwrite(
  overall_summary,
  file.path(
    output_dir,
    "overall_summary.csv"
  )
)

fwrite(
  overlap_coding_summary,
  file.path(
    output_dir,
    "overlap_coding_summary.csv"
  )
)

fwrite(
  coding_preference_summary,
  file.path(
    output_dir,
    "shared_snp_coding_preference.csv"
  )
)

fwrite(
  agreement_summary,
  file.path(
    output_dir,
    "agreement_categories.csv"
  )
)

fwrite(
  agreement_statistics,
  file.path(
    output_dir,
    "agreement_statistics.csv"
  )
)

fwrite(
  cs_count_table,
  file.path(
    output_dir,
    "credible_set_count_table.csv"
  )
)

fwrite(
  coding_summary,
  file.path(
    output_dir,
    "mixed_coding_summary.csv"
  )
)

fwrite(
  coding_pattern_summary,
  file.path(
    output_dir,
    "coding_patterns.csv"
  )
)

fwrite(
  adr_regions,
  file.path(
    output_dir,
    "adr_gene_tissue_pairs.csv"
  )
)

fwrite(
  permutation_summary,
  file.path(
    output_dir,
    "permutation_summary.csv"
  )
)

fwrite(
  tissue_summary,
  file.path(
    output_dir,
    "tissue_summary.csv"
  )
)

fwrite(
  top_discordant_candidates,
  file.path(
    output_dir,
    "top_discordant_candidates.csv"
  )
)

fwrite(
  res_idx,
  file.path(
    output_dir,
    "primary_analysis_set.csv"
  )
)


# ============================================================
# Results-sentence template
# ============================================================

same_cs_percentage <- safe_percent(
  res_idx$ncs_susie ==
    res_idx$ncs_susie_mix
)

one_cs_overlap_percentage <- safe_percent(
  one_cs_each$overlap_snp > 0
)

one_cs_no_overlap_percentage <- safe_percent(
  one_cs_each$overlap_snp == 0
)

primary_overlap_summary <- overlap_coding_summary[
  subset == primary_label
]

cat("\nSuggested results sentence:\n\n")

cat(
  sprintf(
    paste0(
      "After restricting the analysis to gene-tissue pairs with ",
      "a marginal association P value below %.1e and a mean read ",
      "count of at least %d, %s gene-tissue pairs representing %s ",
      "genes were retained. The additive and mixed-coding SuSiE ",
      "models inferred the same number of credible sets in %.1f%% ",
      "of analyses. Among the %s gene-tissue pairs for which both ",
      "models inferred exactly one credible set, %.1f%% shared at ",
      "least one biological SNP, whereas %.1f%% showed no SNP-level ",
      "overlap. Across additive credible-set SNP occurrences, %.1f%% ",
      "were retained by the mixed model under any coding and %.1f%% ",
      "were retained specifically under additive coding. Among SNPs ",
      "shared by both analyses, %.1f%% received their largest ",
      "mixed-model predictor PIP under dominant or recessive coding."
    ),
    association_threshold,
    minimum_mean_reads,
    format(nrow(res_idx), big.mark = ","),
    format(uniqueN(res_idx$gene), big.mark = ","),
    same_cs_percentage,
    format(nrow(one_cs_each), big.mark = ","),
    one_cs_overlap_percentage,
    one_cs_no_overlap_percentage,
    primary_overlap_summary$
      pct_additive_cs_snps_retained_any_coding,
    primary_overlap_summary$
      pct_additive_cs_snps_retained_as_additive,
    primary_overlap_summary$
      pct_shared_snps_preferred_nonadditive
  ),
  "\n"
)
