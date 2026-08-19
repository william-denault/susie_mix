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

required_columns <- c(
  "gene",
  "tissue",
  "min_pv",
  "mean_count",
  "ncs_susie",
  "ncs_susie_mix",
  "overlap",
  "perm_cs_susie",
  "perm_cs_susie_mix",
  "dif_elbo",
  "dif_elbo_perm",
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
    paste(missing_columns, collapse = ", ")
  )
}


# ============================================================
# Define the primary analysis set
# ============================================================

res_strong <- res[
 # is.finite(min_pv) &
    min_pv < association_threshold
]

res_idx <- res[
  #is.finite(min_pv) &
    min_pv < association_threshold &
    is.finite(mean_count) &
    mean_count >= minimum_mean_reads &
    !is.na(ncs_susie) &
    !is.na(ncs_susie_mix) &
    !is.na(overlap)
]

cat("All gene-tissue pairs:", nrow(res), "\n")
cat(
  "Strong association pairs:",
  nrow(res_strong),
  "\n"
)
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
      "P < 5e-8"
    ),
    describe_subset(
      res_idx,
      "P < 5e-8 and mean reads >= 100"
    )
  )
)

print(overall_summary)


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
      overlap > 0,
    "Same CS count; SNP overlap",

    ncs_susie == ncs_susie_mix &
      overlap == 0,
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

# Retain categories with zero observations.
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
  percentage := if (sum(count) > 0L) {
    100 * count / sum(count)
  } else {
    NA_real_
  }
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
    percentage = if (denominator > 0L) {
      100 * count / denominator
    } else {
      NA_real_
    }
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
      "Any s
      hared SNP when both models reported CSs",
      sum(
        both_models_cs$overlap > 0
      ),
      nrow(both_models_cs)
    ),

    make_metric(
      "No shared SNP when both models reported CSs",
      sum(
        both_models_cs$overlap == 0
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
        one_cs_each$overlap > 0
      ),
      nrow(one_cs_each)
    ),

    make_metric(
      "One CS each with no shared SNP",
      sum(
        one_cs_each$overlap == 0
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
  percentage := 100 * count / sum(count)
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
  percentage := (
    100 * credible_sets /
      sum(credible_sets)
  )
]

print(coding_summary)


# Coding combinations found within each gene-tissue analysis.
res_idx[
  ,
  coding_pattern := paste0(
    ifelse(n_add > 0, "A", ""),
    ifelse(n_rec > 0, "R", ""),
    ifelse(n_dom > 0, "D", "")
  )
]

res_idx[
  coding_pattern == "",
  coding_pattern := "none"
]

res_idx[
  ,
  number_coding_types := (
    (n_add > 0) +
      (n_rec > 0) +
      (n_dom > 0)
  )
]

coding_pattern_summary <- res_idx[
  ,
  .(
    count = .N
  ),
  by = coding_pattern
][
  ,
  percentage := 100 * count / sum(count)
][
  order(-count)
]

print(coding_pattern_summary)

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
    100 *
      analyses_with_at_least_one_cs /
      nrow(res_idx)
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
      overlap == 0
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
    "overlap"
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


# Separate candidates by mixed-model coding.
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
# Tissue-specific disagreement
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
      overlap == 0
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
    )
  ),
  by = tissue
]

tissue_summary[
  ,
  pct_different_cs_number := (
    100 *
      different_cs_number /
      total
  )
]

tissue_summary[
  ,
  pct_no_snp_overlap := (
    100 *
      no_snp_overlap /
      total
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

par(
  mfrow = c(2, 2),
  mar = c(5, 5, 4, 1)
)

elbo_values <- -2* (res_idx$log_lik_add-res_idx$log_lik_mix)

hist(
  elbo_values,
  xlim=c(-20,50),
  nclass = 1000,
  main = paste0(
    "-2( log lik additive - log lik mixed)"
  ),
  xlab = "-2( log lik additive - log lik mixed)",

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
 # breaks = "FD",
 xlim=c(-20,50),
 nclass = 1000,
  main = paste0(
    "-2( log lik additive - log lik mixed)\n permuted"
  ),
  xlab = "-2( log lik additive - log lik mixed)",

  border = "gray80"
)

abline(
  v = 0,
  col = "red",
  lty = 2,
  lwd = 2
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
  "Same CS count;\nSNP overlap",
  "Different\nCS count",
  "Same CS count;\nno SNP overlap",
  "Neither model\nreported a CS",
  "Additive\nmodel only",
  "Mixed\nmodel only"
)

par(mar = c(9, 5, 4, 1))

barplot(
  agreement_summary$percentage,
  names.arg = agreement_plot_labels,
  #las = 2,
  cex.names = 0.70,
  ylab = "Percentage",
  main = "Agreement between fine-mapping models",
  col = "steelblue"
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

fwrite(
  overall_summary,
  file.path(
    output_dir,
    "overall_summary.csv"
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
  one_cs_each$overlap > 0
)

one_cs_no_overlap_percentage <- safe_percent(
  one_cs_each$overlap == 0
)

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
      "overlap."
    ),
    association_threshold,
    minimum_mean_reads,
    format(nrow(res_idx), big.mark = ","),
    format(uniqueN(res_idx$gene), big.mark = ","),
    same_cs_percentage,
    format(nrow(one_cs_each), big.mark = ","),
    one_cs_overlap_percentage,
    one_cs_no_overlap_percentage
  ),
  "\n"
)
