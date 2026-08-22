library(data.table)

load(
  "/project2/mstephens/wdenault/susie_mix/res_summary.RData"
)

load(
  "/project2/mstephens/wdenault/susie_mix/res_cs_summary.RData"
)

res <- as.data.table(res_summary)

if (exists("res_cs_summary")) {
  cs_res <- as.data.table(res_cs_summary)
} else if (exists("cs_summary")) {
  # Backward compatibility with an earlier object name.
  cs_res <- as.data.table(cs_summary)
} else {
  stop(
    "res_cs_summary.RData does not contain res_cs_summary."
  )
}

# ============================================================
# Analysis thresholds
# ============================================================

association_threshold <- 1e-8
minimum_mean_reads <- 100
tss_plot_limit_kb <- 200
tss_bin_width_kb <- 10

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

required_cs_columns <- c(
  "gene",
  "tissue",
  "model",
  "lead_snp",
  "lead_coding",
  "cs_snps",
  "min_pv",
  "mean_count",
  "distance_to_tss_kb"
)

missing_cs_columns <- setdiff(
  required_cs_columns,
  names(cs_res)
)

if (length(missing_cs_columns) > 0L) {
  stop(
    "CS summary is missing: ",
    paste(missing_cs_columns, collapse = ", "),
    ". Rerun the CS-centric summarizer."
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


parse_semicolon_set <- function(x) {

  if (
    length(x) != 1L ||
    is.na(x) ||
    !nzchar(trimws(x))
  ) {
    return(character())
  }

  values <- trimws(
    strsplit(
      x,
      ";",
      fixed = TRUE
    )[[1L]]
  )

  unique(
    values[
      !is.na(values) &
        nzchar(values)
    ]
  )
}


collapse_sorted_values <- function(x) {

  values <- trimws(
    as.character(x)
  )

  values <- values[
    !is.na(values) &
      nzchar(values)
  ]

  paste(
    sort(values),
    collapse = ";"
  )
}


summarize_tss_distribution <- function(
    x,
    plot_limit_kb,
    bin_width_kb) {

  model_levels <- c(
    "SuSiE",
    "SuSiE-mix"
  )

  bin_centers <- seq(
    -plot_limit_kb,
    plot_limit_kb,
    by = bin_width_kb
  )

  bin_breaks <- c(
    bin_centers - bin_width_kb / 2,
    tail(bin_centers, 1) + bin_width_kb / 2
  )

  rbindlist(
    lapply(
      model_levels,
      function(model_name) {

        model_distance <- x[
          model == model_name &
            is.finite(distance_to_tss_kb),
          distance_to_tss_kb
        ]

        in_window <- (
          model_distance >= min(bin_breaks) &
            model_distance <= max(bin_breaks)
        )

        plotted_distance <- model_distance[
          in_window
        ]

        bin_count <- hist(
          plotted_distance,
          breaks = bin_breaks,
          plot = FALSE,
          include.lowest = TRUE,
          right = FALSE
        )$counts

        total_in_window <- length(
          plotted_distance
        )

        data.table(
          model = model_name,
          distance_to_tss_kb = bin_centers,
          bin_lower_kb = head(
            bin_breaks,
            -1
          ),
          bin_upper_kb = tail(
            bin_breaks,
            -1
          ),
          count = bin_count,
          proportion = if (
            total_in_window > 0L
          ) {
            bin_count / total_in_window
          } else {
            NA_real_
          },
          percentage = if (
            total_in_window > 0L
          ) {
            100 * bin_count / total_in_window
          } else {
            NA_real_
          },
          n_cs_total = length(model_distance),
          n_cs_in_window = total_in_window
        )
      }
    )
  )
}


summarize_tissue_additivity <- function(x) {

  recognized_codings <- c(
    "additive",
    "dominant",
    "recessive"
  )

  mixed_cs <- x[
    model == "SuSiE-mix"
  ]

  tissue_additivity <- mixed_cs[
    ,
    .(
      total_mixed_cs = .N,
      classified_mixed_cs = sum(
        lead_coding %in% recognized_codings
      ),
      additive_cs = sum(
        lead_coding == "additive",
        na.rm = TRUE
      ),
      dominant_cs = sum(
        lead_coding == "dominant",
        na.rm = TRUE
      ),
      recessive_cs = sum(
        lead_coding == "recessive",
        na.rm = TRUE
      ),
      unclassified_cs = sum(
        is.na(lead_coding) |
          !lead_coding %in% recognized_codings
      )
    ),
    by = tissue
  ]

  if (nrow(tissue_additivity) == 0L) {
    return(data.table())
  }

  tissue_additivity[
    ,
    nonadditive_cs := dominant_cs + recessive_cs
  ]

  total_additive <- sum(
    tissue_additivity$additive_cs
  )

  total_nonadditive <- sum(
    tissue_additivity$nonadditive_cs
  )

  total_classified <- (
    total_additive + total_nonadditive
  )

  overall_additive_proportion <- if (
    total_classified > 0L
  ) {
    total_additive / total_classified
  } else {
    NA_real_
  }

  tissue_additivity[
    ,
    additive_proportion := safe_ratio_percent(
      additive_cs,
      classified_mixed_cs
    ) / 100
  ]

  z_value <- qnorm(0.975)

  tissue_additivity[
    ,
    `:=`(
      additive_ci_lower = (
        (
          additive_proportion +
            z_value^2 / (2 * classified_mixed_cs)
        ) /
          (1 + z_value^2 / classified_mixed_cs) -
          z_value * sqrt(
            additive_proportion *
              (1 - additive_proportion) /
              classified_mixed_cs +
              z_value^2 /
              (4 * classified_mixed_cs^2)
          ) /
          (1 + z_value^2 / classified_mixed_cs)
      ),
      additive_ci_upper = (
        (
          additive_proportion +
            z_value^2 / (2 * classified_mixed_cs)
        ) /
          (1 + z_value^2 / classified_mixed_cs) +
          z_value * sqrt(
            additive_proportion *
              (1 - additive_proportion) /
              classified_mixed_cs +
              z_value^2 /
              (4 * classified_mixed_cs^2)
          ) /
          (1 + z_value^2 / classified_mixed_cs)
      )
    )
  ]

  tissue_additivity[
    !is.finite(additive_ci_lower),
    additive_ci_lower := NA_real_
  ]

  tissue_additivity[
    !is.finite(additive_ci_upper),
    additive_ci_upper := NA_real_
  ]

  test_results <- rbindlist(
    lapply(
      seq_len(nrow(tissue_additivity)),
      function(i) {

        tissue_additive <- tissue_additivity$additive_cs[i]
        tissue_nonadditive <- tissue_additivity$nonadditive_cs[i]
        other_additive <- total_additive - tissue_additive
        other_nonadditive <- total_nonadditive - tissue_nonadditive

        contingency_table <- matrix(
          c(
            tissue_additive,
            tissue_nonadditive,
            other_additive,
            other_nonadditive
          ),
          nrow = 2,
          byrow = TRUE
        )

        fisher_result <- tryCatch(
          fisher.test(
            contingency_table,
            conf.int = TRUE
          ),
          error = function(e) NULL
        )

        if (is.null(fisher_result)) {
          return(data.table(
            additive_odds_ratio = NA_real_,
            odds_ratio_ci_lower = NA_real_,
            odds_ratio_ci_upper = NA_real_,
            fisher_p_value = NA_real_
          ))
        }

        fisher_estimate <- as.numeric(
          fisher_result$estimate
        )

        if (length(fisher_estimate) == 0L) {
          fisher_estimate <- NA_real_
        }

        fisher_confidence_interval <- as.numeric(
          fisher_result$conf.int
        )

        if (length(fisher_confidence_interval) < 2L) {
          fisher_confidence_interval <- c(
            NA_real_,
            NA_real_
          )
        }

        data.table(
          additive_odds_ratio = fisher_estimate[1],
          odds_ratio_ci_lower = fisher_confidence_interval[1],
          odds_ratio_ci_upper = fisher_confidence_interval[2],
          fisher_p_value = fisher_result$p.value
        )
      }
    )
  )

  tissue_additivity <- cbind(
    tissue_additivity,
    test_results
  )

  tissue_additivity[
    ,
    `:=`(
      overall_additive_proportion = overall_additive_proportion,
      additive_difference = (
        additive_proportion -
          overall_additive_proportion
      ),
      standardized_additive_residual = (
        (
          additive_cs -
            classified_mixed_cs *
            overall_additive_proportion
        ) /
          sqrt(
            classified_mixed_cs *
              overall_additive_proportion *
              (1 - overall_additive_proportion)
          )
      ),
      other_tissues_additive_proportion = safe_ratio_percent(
        total_additive - additive_cs,
        total_classified - classified_mixed_cs
      ) / 100,
      fisher_fdr = p.adjust(
        fisher_p_value,
        method = "BH"
      )
    )
  ]

  tissue_additivity[
    !is.finite(standardized_additive_residual),
    standardized_additive_residual := NA_real_
  ]

  tissue_additivity[
    ,
    direction_vs_overall := fifelse(
      additive_difference > 0,
      "over-additive",
      fifelse(
        additive_difference < 0,
        "under-additive",
        "at standard"
      )
    )
  ]

  tissue_additivity[
    ,
    additivity_classification := fifelse(
      is.finite(fisher_fdr) &
        fisher_fdr < 0.05,
      direction_vs_overall,
      "not significantly different"
    )
  ]

  tissue_additivity[
    ,
    `:=`(
      additive_percentage = 100 * additive_proportion,
      additive_ci_lower_percentage = 100 * additive_ci_lower,
      additive_ci_upper_percentage = 100 * additive_ci_upper,
      overall_additive_percentage = (
        100 * overall_additive_proportion
      ),
      additive_difference_percentage_points = (
        100 * additive_difference
      ),
      other_tissues_additive_percentage = (
        100 * other_tissues_additive_proportion
      )
    )
  ]

  setorder(
    tissue_additivity,
    -additive_difference
  )

  tissue_additivity
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

cs_idx <- cs_res[
  is.finite(min_pv) &
    min_pv < association_threshold &
    is.finite(mean_count) &
    mean_count >= minimum_mean_reads &
    model %in% c(
      "SuSiE",
      "SuSiE-mix"
    )
]

tss_distance_summary <- summarize_tss_distribution(
  x = cs_idx,
  plot_limit_kb = tss_plot_limit_kb,
  bin_width_kb = tss_bin_width_kb
)

finite_tss_distance <- cs_idx[
  is.finite(distance_to_tss_kb),
  distance_to_tss_kb
]

if (
  length(finite_tss_distance) > 0L &&
  all(finite_tss_distance >= 0)
) {
  warning(
    paste0(
      "All CS-to-TSS distances are nonnegative. ",
      "The workhorse may have saved absolute rather than signed ",
      "distances, so the TSS plot will be one-sided."
    )
  )
}

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
cat(
  "CS rows in primary set:",
  nrow(cs_idx),
  "\n"
)
cat(
  "CS rows with finite TSS distance:",
  length(finite_tss_distance),
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

additive_lead_summary <- cs_idx[
  model == "SuSiE",
  .(
    additive_cs_rows = .N,
    additive_valid_lead_count = sum(
      !is.na(lead_snp) &
        nzchar(trimws(as.character(lead_snp)))
    ),
    additive_lead_snp_set = collapse_sorted_values(
      lead_snp
    )
  ),
  by = .(
    gene,
    tissue
  )
]

mixed_lead_summary <- cs_idx[
  model == "SuSiE-mix",
  .(
    mixed_cs_rows = .N,
    mixed_valid_lead_count = sum(
      !is.na(lead_snp) &
        nzchar(trimws(as.character(lead_snp)))
    ),
    mixed_lead_snp_set = collapse_sorted_values(
      lead_snp
    )
  ),
  by = .(
    gene,
    tissue
  )
]

res_idx[
  ,
  `:=`(
    additive_cs_rows = 0L,
    additive_valid_lead_count = 0L,
    additive_lead_snp_set = "",
    mixed_cs_rows = 0L,
    mixed_valid_lead_count = 0L,
    mixed_lead_snp_set = ""
  )
]

res_idx[
  additive_lead_summary,
  on = .(
    gene,
    tissue
  ),
  `:=`(
    additive_cs_rows = i.additive_cs_rows,
    additive_valid_lead_count = i.additive_valid_lead_count,
    additive_lead_snp_set = i.additive_lead_snp_set
  )
]

res_idx[
  mixed_lead_summary,
  on = .(
    gene,
    tissue
  ),
  `:=`(
    mixed_cs_rows = i.mixed_cs_rows,
    mixed_valid_lead_count = i.mixed_valid_lead_count,
    mixed_lead_snp_set = i.mixed_lead_snp_set
  )
]

lead_summary_mismatch <- res_idx[
  additive_cs_rows != ncs_susie |
    mixed_cs_rows != ncs_susie_mix |
    additive_valid_lead_count != additive_cs_rows |
    mixed_valid_lead_count != mixed_cs_rows
]

if (nrow(lead_summary_mismatch) > 0L) {
  stop(
    paste0(
      "The CS-centric and gene-tissue summaries disagree for ",
      nrow(lead_summary_mismatch),
      " analysis rows, or a CS is missing its biological lead SNP. ",
      "Regenerate both summaries from the same result files."
    )
  )
}

res_idx[
  ,
  same_lead_snp_set := (
    ncs_susie > 0 &
      ncs_susie == ncs_susie_mix &
      additive_lead_snp_set == mixed_lead_snp_set
  )
]

agreement_levels <- c(
  "Same CS count; same lead SNP(s)",
  "Same CS count; different lead SNP(s); SNP overlap",
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
      overlap_snp > 0 &
      same_lead_snp_set,
    "Same CS count; same lead SNP(s)",

    ncs_susie == ncs_susie_mix &
      overlap_snp > 0 &
      !same_lead_snp_set,
    "Same CS count; different lead SNP(s); SNP overlap",

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


# ============================================================
# Detailed agreement when each model reported exactly one CS
# ============================================================

one_cs_keys <- unique(
  one_cs_each[
    ,
    .(
      gene,
      tissue
    )
  ]
)

one_cs_cs <- cs_idx[
  one_cs_keys,
  on = .(
    gene,
    tissue
  ),
  nomatch = 0L
]

additive_one_cs <- one_cs_cs[
  model == "SuSiE",
  .(
    gene,
    tissue,
    additive_lead_snp = as.character(lead_snp),
    additive_cs_snps = as.character(cs_snps)
  )
]

mixed_one_cs <- one_cs_cs[
  model == "SuSiE-mix",
  .(
    gene,
    tissue,
    mixed_lead_snp = as.character(lead_snp),
    mixed_cs_snps = as.character(cs_snps)
  )
]

if (
  anyDuplicated(
    additive_one_cs[
      ,
      .(
        gene,
        tissue
      )
    ]
  ) > 0L ||
  anyDuplicated(
    mixed_one_cs[
      ,
      .(
        gene,
        tissue
      )
    ]
  ) > 0L
) {
  stop(
    paste0(
      "The CS-centric summary contains more than one CS row for a ",
      "gene-tissue-model combination classified as having one CS."
    )
  )
}

one_cs_overlap_detail <- merge(
  additive_one_cs,
  mixed_one_cs,
  by = c(
    "gene",
    "tissue"
  ),
  all = FALSE
)

if (nrow(one_cs_overlap_detail) != nrow(one_cs_each)) {
  stop(
    paste0(
      "Could not match one additive and one mixed-model CS row for ",
      "every one-CS gene-tissue pair. Expected ",
      nrow(one_cs_each),
      " pairs but matched ",
      nrow(one_cs_overlap_detail),
      ". Rerun the CS-centric summarizer from the same result files."
    )
  )
}

one_cs_comparisons <- Map(
  function(
    additive_cs_string,
    mixed_cs_string,
    additive_lead,
    mixed_lead) {

    additive_snps <- parse_semicolon_set(
      additive_cs_string
    )

    mixed_snps <- parse_semicolon_set(
      mixed_cs_string
    )

    shared_snps <- intersect(
      additive_snps,
      mixed_snps
    )

    union_snps <- union(
      additive_snps,
      mixed_snps
    )

    full_overlap <- (
      length(additive_snps) > 0L &&
        length(mixed_snps) > 0L &&
        setequal(
          additive_snps,
          mixed_snps
        )
    )

    same_lead <- (
      !is.na(additive_lead) &&
        nzchar(additive_lead) &&
        !is.na(mixed_lead) &&
        nzchar(mixed_lead) &&
        identical(
          additive_lead,
          mixed_lead
        )
    )

    data.table(
      additive_cs_size_snps = length(additive_snps),
      mixed_cs_size_snps = length(mixed_snps),
      shared_snp_count = length(shared_snps),
      union_snp_count = length(union_snps),
      jaccard_similarity = if (
        length(union_snps) > 0L
      ) {
        length(shared_snps) / length(union_snps)
      } else {
        NA_real_
      },
      same_lead_snp = same_lead,
      full_snp_overlap = full_overlap,
      overlap_class = if (full_overlap) {
        "Full SNP overlap"
      } else if (length(shared_snps) > 0L) {
        "Partial SNP overlap"
      } else {
        "No SNP overlap"
      }
    )
  },
  one_cs_overlap_detail$additive_cs_snps,
  one_cs_overlap_detail$mixed_cs_snps,
  one_cs_overlap_detail$additive_lead_snp,
  one_cs_overlap_detail$mixed_lead_snp
)

one_cs_overlap_detail <- cbind(
  one_cs_overlap_detail,
  rbindlist(one_cs_comparisons)
)

one_cs_overlap_detail[
  ,
  lead_agreement := ifelse(
    same_lead_snp,
    "Same lead SNP",
    "Different lead SNP"
  )
]

one_cs_overlap_levels <- c(
  "Full SNP overlap",
  "Partial SNP overlap",
  "No SNP overlap"
)

one_cs_lead_levels <- c(
  "Same lead SNP",
  "Different lead SNP"
)

one_cs_overlap_summary <- CJ(
  overlap_class = one_cs_overlap_levels,
  lead_agreement = one_cs_lead_levels,
  unique = TRUE
)

one_cs_overlap_summary <- one_cs_overlap_detail[
  ,
  .(
    count = .N
  ),
  by = .(
    overlap_class,
    lead_agreement
  )
][
  one_cs_overlap_summary,
  on = .(
    overlap_class,
    lead_agreement
  )
]

one_cs_overlap_summary[
  is.na(count),
  count := 0L
]

one_cs_overlap_summary[
  ,
  percentage := safe_ratio_percent(
    count,
    nrow(one_cs_overlap_detail)
  )
]

one_cs_overlap_summary[
  ,
  `:=`(
    overlap_class = factor(
      overlap_class,
      levels = one_cs_overlap_levels
    ),
    lead_agreement = factor(
      lead_agreement,
      levels = one_cs_lead_levels
    )
  )
]

setorder(
  one_cs_overlap_summary,
  overlap_class,
  lead_agreement
)

print(one_cs_overlap_summary)


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
      "Same CS count and same biological lead SNP set",
      sum(
        as.character(res_idx$agreement_category) ==
          "Same CS count; same lead SNP(s)"
      ),
      nrow(res_idx)
    ),

    make_metric(
      paste0(
        "Same CS count, different biological lead SNP set, ",
        "and SNP overlap"
      ),
      sum(
        as.character(res_idx$agreement_category) ==
          paste0(
            "Same CS count; different lead SNP(s); ",
            "SNP overlap"
          )
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
    ),

    make_metric(
      "One CS each with the same biological lead SNP",
      sum(
        one_cs_overlap_detail$same_lead_snp
      ),
      nrow(one_cs_overlap_detail)
    ),

    make_metric(
      "One CS each with full biological-SNP overlap",
      sum(
        one_cs_overlap_detail$full_snp_overlap
      ),
      nrow(one_cs_overlap_detail)
    ),

    make_metric(
      "One CS each with partial biological-SNP overlap",
      sum(
        one_cs_overlap_detail$overlap_class ==
          "Partial SNP overlap"
      ),
      nrow(one_cs_overlap_detail)
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
# Tissue-specific coding patterns and additive enrichment
# ============================================================

tissues_for_plot <- sort(
  unique(res_idx$tissue)
)

tissue_coding_pattern_summary <- merge(
  CJ(
    tissue = tissues_for_plot,
    coding_pattern = coding_pattern_levels,
    unique = TRUE
  ),
  res_idx[
    ,
    .(count = .N),
    by = .(
      tissue,
      coding_pattern
    )
  ],
  by = c(
    "tissue",
    "coding_pattern"
  ),
  all.x = TRUE,
  sort = FALSE
)

tissue_coding_pattern_summary[
  is.na(count),
  count := 0L
]

tissue_coding_pattern_summary[
  ,
  `:=`(
    total_gene_tissue_pairs = sum(count),
    percentage = safe_ratio_percent(
      count,
      sum(count)
    ),
    pattern_order = match(
      coding_pattern,
      coding_pattern_levels
    )
  ),
  by = tissue
]

setorder(
  tissue_coding_pattern_summary,
  tissue,
  pattern_order
)

tissue_additivity_summary <- summarize_tissue_additivity(
  cs_idx
)

tissue_tss_distance_summary <- rbindlist(
  lapply(
    tissues_for_plot,
    function(tissue_name) {

      tissue_tss <- summarize_tss_distribution(
        x = cs_idx[
          tissue == tissue_name
        ],
        plot_limit_kb = tss_plot_limit_kb,
        bin_width_kb = tss_bin_width_kb
      )

      tissue_tss[
        ,
        tissue := tissue_name
      ]

      setcolorder(
        tissue_tss,
        c(
          "tissue",
          setdiff(
            names(tissue_tss),
            "tissue"
          )
        )
      )

      tissue_tss
    }
  )
)

print(tissue_additivity_summary)

if (nrow(tissue_additivity_summary) > 0L) {

  cat(
    "Overall additive lead-CS percentage:",
    unique(
      tissue_additivity_summary$
        overall_additive_percentage
    )[1],
    "\n"
  )

  cat("Tissues classified as over- or under-additive:\n")

  print(
    tissue_additivity_summary[
      additivity_classification %in%
        c(
          "over-additive",
          "under-additive"
        ),
      .(
        tissue,
        additive_cs,
        classified_mixed_cs,
        additive_percentage,
        additive_difference_percentage_points,
        additive_odds_ratio,
        fisher_fdr,
        additivity_classification
      )
    ]
  )
}

tissue_plot_order <- unique(
  c(
    tissue_additivity_summary$tissue,
    tissues_for_plot
  )
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


plot_tissue_coding_patterns <- function(
    tissue_name,
    pattern_summary,
    additivity_summary) {

  plot_data <- pattern_summary[
    tissue == tissue_name
  ]

  if (nrow(plot_data) == 0L) {
    plot.new()
    title(
      paste0(
        tissue_name,
        ": coding patterns"
      ),
      cex.main = 0.85
    )
    text(
      0.5,
      0.5,
      "No gene-tissue analyses"
    )
    return(invisible(NULL))
  }

  setorder(
    plot_data,
    pattern_order
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

  par(
    mar = c(3.8, 4.2, 2.8, 1),
    mgp = c(2.4, 0.7, 0)
  )

  bar_positions <- barplot(
    plot_data$percentage,
    names.arg = plot_data$coding_pattern,
    ylim = c(0, 100),
    ylab = "Gene-tissue pairs (%)",
    main = paste0(
      tissue_name,
      ": mixed-model coding patterns"
    ),
    col = unname(
      pattern_colors[
        plot_data$coding_pattern
      ]
    ),
    border = NA,
    cex.names = 0.7,
    cex.axis = 0.75,
    cex.lab = 0.8,
    cex.main = 0.85
  )

  nonzero_pattern <- plot_data$count > 0L

  text(
    x = bar_positions[nonzero_pattern],
    y = plot_data$percentage[nonzero_pattern],
    labels = plot_data$count[nonzero_pattern],
    pos = 3,
    cex = 0.6
  )

  additivity_row <- additivity_summary[
    tissue == tissue_name
  ]

  if (nrow(additivity_row) > 0L) {

    fdr_label <- if (
      is.finite(additivity_row$fisher_fdr[1])
    ) {
      format.pval(
        additivity_row$fisher_fdr[1],
        digits = 2,
        eps = 0.001
      )
    } else {
      "NA"
    }

    mtext(
      sprintf(
        paste0(
          "Additive lead CS: %.1f%%; overall: %.1f%%; ",
          "%s (FDR %s)"
        ),
        additivity_row$additive_percentage[1],
        additivity_row$overall_additive_percentage[1],
        additivity_row$direction_vs_overall[1],
        fdr_label
      ),
      side = 3,
      line = 0.15,
      adj = 1,
      cex = 0.58
    )
  }

  invisible(bar_positions)
}


plot_tss_distribution <- function(
    x,
    plot_limit_kb,
    main_title = "CS lead SNPs around the TSS",
    compact = FALSE) {

  model_levels <- c(
    "SuSiE",
    "SuSiE-mix"
  )

  model_colors <- c(
    SuSiE = "#E69F00",
    `SuSiE-mix` = "#0072B2"
  )

  valid_proportion <- is.finite(
    x$proportion
  )

  if (!any(valid_proportion)) {
    plot.new()
    title(
      main_title,
      cex.main = if (compact) 0.85 else 1
    )
    text(
      0.5,
      0.5,
      "No finite CS-to-TSS distances"
    )
    return(invisible(NULL))
  }

  y_max <- max(
    x$proportion[valid_proportion]
  )

  y_max <- max(
    0.01,
    1.12 * y_max
  )

  par(
    mar = if (compact) {
      c(3.8, 4.2, 2.8, 1)
    } else {
      c(5, 5, 4, 1)
    },
    mgp = if (compact) {
      c(2.4, 0.7, 0)
    } else {
      c(2.8, 0.8, 0)
    }
  )

  plot(
    NA_real_,
    NA_real_,
    xlim = c(
      -plot_limit_kb,
      plot_limit_kb
    ),
    ylim = c(0, y_max),
    xaxs = "i",
    yaxs = "i",
    xaxt = "n",
    bty = "l",
    xlab = "Distance to TSS (kb)",
    ylab = "Proportion of credible sets",
    main = main_title,
    cex.main = if (compact) 0.85 else 1,
    cex.lab = if (compact) 0.8 else 1,
    cex.axis = if (compact) 0.75 else 1
  )

  axis(
    side = 1,
    at = seq(
      -plot_limit_kb,
      plot_limit_kb,
      length.out = 5
    ),
    cex.axis = if (compact) 0.75 else 1
  )

  abline(
    v = 0,
    col = "gray65",
    lty = 3,
    lwd = 1.5
  )

  for (model_name in model_levels) {

    model_index <- (
      x[["model"]] == model_name
    )

    lines(
      x[["distance_to_tss_kb"]][model_index],
      x[["proportion"]][model_index],
      col = model_colors[[model_name]],
      lwd = 2.5
    )
  }

  legend_labels <- vapply(
    model_levels,
    function(model_name) {

      model_index <- (
        x[["model"]] == model_name
      )

      model_n <- unique(
        x[["n_cs_in_window"]][model_index]
      )

      if (length(model_n) == 0L) {
        model_n <- 0L
      }

      paste0(
        model_name,
        " (n = ",
        model_n[1],
        ")"
      )
    },
    character(1)
  )

  legend(
    "topright",
    legend = legend_labels,
    col = unname(
      model_colors[model_levels]
    ),
    lwd = 2.5,
    bty = "n",
    cex = if (compact) 0.65 else 0.8
  )

  invisible(NULL)
}


plot_tissue_specific_page <- function(
    tissue_names,
    pattern_summary,
    tss_summary,
    additivity_summary,
    plot_limit_kb) {

  old_par <- par(no.readonly = TRUE)
  on.exit(
    par(old_par),
    add = TRUE
  )

  par(
    mfrow = c(4, 2),
    mar = c(3.8, 4.2, 2.8, 1)
  )

  for (tissue_name in tissue_names) {

    plot_tissue_coding_patterns(
      tissue_name = tissue_name,
      pattern_summary = pattern_summary,
      additivity_summary = additivity_summary
    )

    plot_tss_distribution(
      x = tss_summary[
        tissue == tissue_name
      ],
      plot_limit_kb = plot_limit_kb,
      main_title = paste0(
        tissue_name,
        ": CS lead SNPs around the TSS"
      ),
      compact = TRUE
    )
  }

  missing_tissue_rows <- 4L - length(tissue_names)

  if (missing_tissue_rows > 0L) {
    for (i in seq_len(2L * missing_tissue_rows)) {
      plot.new()
    }
  }

  invisible(NULL)
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
  "Same count;\nsame lead",
  "Same count;\ndifferent lead;\nSNP overlap",
  "Different\nCS count",
  "Same count;\nno overlap",
  "Neither\nmodel",
  "Additive\nonly",
  "Mixed\nonly"
)

agreement_plot_colors <- c(
  "#0072B2",
  "#56B4E9",
  "#E69F00",
  "#D55E00",
  "#999999",
  "#CC79A7",
  "#009E73"
)

plot_agreement_summary <- function(summary_table) {

  valid_percentage <- is.finite(
    summary_table$percentage
  )

  if (!any(valid_percentage)) {
    plot.new()
    title("Agreement between fine-mapping models")
    text(
      0.5,
      0.5,
      "No gene-tissue analyses"
    )
    return(invisible(NULL))
  }

  bar_heights <- ifelse(
    valid_percentage,
    summary_table$percentage,
    0
  )

  y_max <- max(
    1,
    1.14 * max(bar_heights)
  )

  bar_midpoints <- barplot(
    bar_heights,
    axisnames = FALSE,
    ylab = "Percentage of gene-tissue pairs",
    main = "Agreement between fine-mapping models",
    col = agreement_plot_colors,
    border = NA,
    ylim = c(
      0,
      y_max
    )
  )

  axis(
    side = 1,
    at = bar_midpoints,
    labels = agreement_plot_labels,
    tick = FALSE,
    cex.axis = 0.54,
    gap.axis = -1
  )

  nonzero_category <- summary_table$count > 0L

  text(
    x = bar_midpoints[nonzero_category],
    y = bar_heights[nonzero_category],
    labels = format(
      summary_table$count[nonzero_category],
      big.mark = ","
    ),
    pos = 3,
    cex = 0.62
  )

  invisible(bar_midpoints)
}


par(mar = c(9, 5, 4, 1))

plot_agreement_summary(
  agreement_summary
)

par(mfrow = c(1, 1))

plot_coding_tss_2x2 <- function() {

  old_par <- par(no.readonly = TRUE)
  on.exit(
    par(old_par),
    add = TRUE
  )

  par(
    mfrow = c(2, 2),
    mar = c(5, 5, 4, 1)
  )

  par(mar = c(9, 5, 4, 1))

  plot_agreement_summary(
    agreement_summary
  )

  par(mar = c(8, 5, 4, 1))

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

  } else {

    plot.new()
    title("Mixed-model preferred coding")
    text(
      0.5,
      0.5,
      "No shared SNP occurrences"
    )
  }

  plot_coding_pattern_summary(
    coding_pattern_summary
  )

  plot_tss_distribution(
    tss_distance_summary,
    plot_limit_kb = tss_plot_limit_kb
  )

  invisible(NULL)
}


plot_coding_tss_2x2()


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

pdf(
  file.path(
    output_dir,
    "coding_and_tss_summary_2x2.pdf"
  ),
  width = 11,
  height = 8.5
)

plot_coding_tss_2x2()

invisible(dev.off())

if (length(tissue_plot_order) > 0L) {

  pdf(
    file.path(
      output_dir,
      "tissue_specific_coding_and_tss_4x2.pdf"
    ),
    width = 12,
    height = 16,
    onefile = TRUE
  )

  tissue_pages <- split(
    tissue_plot_order,
    ceiling(
      seq_along(tissue_plot_order) / 4
    )
  )

  for (tissue_page in tissue_pages) {
    plot_tissue_specific_page(
      tissue_names = tissue_page,
      pattern_summary = tissue_coding_pattern_summary,
      tss_summary = tissue_tss_distance_summary,
      additivity_summary = tissue_additivity_summary,
      plot_limit_kb = tss_plot_limit_kb
    )
  }

  invisible(dev.off())
}

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
  one_cs_overlap_summary,
  file.path(
    output_dir,
    "one_cs_overlap_and_lead_summary.csv"
  )
)

fwrite(
  one_cs_overlap_detail,
  file.path(
    output_dir,
    "one_cs_overlap_and_lead_detail.csv"
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
  tss_distance_summary,
  file.path(
    output_dir,
    "tss_distance_distribution.csv"
  )
)

fwrite(
  tissue_coding_pattern_summary[
    ,
    setdiff(
      names(tissue_coding_pattern_summary),
      "pattern_order"
    ),
    with = FALSE
  ],
  file.path(
    output_dir,
    "tissue_coding_patterns.csv"
  )
)

fwrite(
  tissue_tss_distance_summary,
  file.path(
    output_dir,
    "tissue_tss_distance_distribution.csv"
  )
)

fwrite(
  tissue_additivity_summary,
  file.path(
    output_dir,
    "tissue_additivity_summary.csv"
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

same_count_same_lead_percentage <- safe_percent(
  as.character(res_idx$agreement_category) ==
    "Same CS count; same lead SNP(s)"
)

same_count_different_lead_overlap_percentage <- safe_percent(
  as.character(res_idx$agreement_category) ==
    paste0(
      "Same CS count; different lead SNP(s); ",
      "SNP overlap"
    )
)

one_cs_same_lead_percentage <- safe_percent(
  one_cs_overlap_detail$same_lead_snp
)

one_cs_full_overlap_percentage <- safe_percent(
  one_cs_overlap_detail$full_snp_overlap
)

one_cs_partial_overlap_percentage <- safe_percent(
  one_cs_overlap_detail$overlap_class ==
    "Partial SNP overlap"
)

one_cs_no_overlap_percentage <- safe_percent(
  one_cs_overlap_detail$overlap_class ==
    "No SNP overlap"
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
      "of analyses; %.1f%% had the same CS count and the same set of ",
      "biological lead SNPs, whereas %.1f%% had the same CS count and ",
      "SNP overlap but different lead-SNP sets. Among the %s ",
      "gene-tissue pairs for which both ",
      "models inferred exactly one credible set, %.1f%% selected the ",
      "same biological lead SNP, %.1f%% had identical biological-SNP ",
      "membership (full overlap), %.1f%% overlapped partially, and ",
      "%.1f%% showed no SNP-level overlap. Across additive credible-set ",
      "SNP occurrences, %.1f%% ",
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
    same_count_same_lead_percentage,
    same_count_different_lead_overlap_percentage,
    format(nrow(one_cs_overlap_detail), big.mark = ","),
    one_cs_same_lead_percentage,
    one_cs_full_overlap_percentage,
    one_cs_partial_overlap_percentage,
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
