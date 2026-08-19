# summarize_susie_results.R

path_res <- "/project2/mstephens/wdenault/susie_mix/results/"

summary_file <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "res_summary.RData"
)

summary_csv <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "res_summary.csv"
)

error_csv <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "res_errors.csv"
)

failed_gene_file <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "failed_genes.txt"
)

result_files <- list.files(
  path_res,
  pattern = "\\.rds$",
  full.names = TRUE
)

cat("Found", length(result_files), "result files.\n")


# ============================================================
# Helper functions
# ============================================================

get_cs <- function(fit) {

  cs <- fit$sets$cs

  if (is.null(cs)) {
    return(list())
  }

  cs
}


max_elbo <- function(fit) {

  x <- fit$elbo
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  max(x)
}


get_log_lik_metric <- function(fit) {

  elbo <- fit$elbo
  elbo <- elbo[is.finite(elbo)]

  if (
    length(elbo) == 0L ||
    is.null(fit$KL)
  ) {
    return(NA_real_)
  }

  tail(elbo, 1) + sum(
    fit$KL,
    na.rm = TRUE
  )
}


value_or_na <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0L
  ) {
    return(NA)
  }

  x[1]
}


safe_percentage <- function(
    numerator,
    denominator) {

  if (
    length(denominator) == 0L ||
    is.na(denominator) ||
    denominator <= 0
  ) {
    return(NA_real_)
  }

  100 * numerator / denominator
}


clean_snp_vector <- function(x) {

  x <- as.character(x)

  unique(
    x[
      !is.na(x) &
        nzchar(x)
    ]
  )
}


validate_predictor_map <- function(
    fit,
    predictor_map,
    map_name) {

  if (!is.data.frame(predictor_map)) {
    stop(
      map_name,
      " is not a data frame."
    )
  }

  required_columns <- c(
    "predictor_index",
    "predictor_name",
    "snp",
    "coding"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(predictor_map)
  )

  if (length(missing_columns) > 0L) {
    stop(
      map_name,
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (nrow(predictor_map) != ncol(fit$alpha)) {
    stop(
      map_name,
      " has ",
      nrow(predictor_map),
      " rows, but the SuSiE fit has ",
      ncol(fit$alpha),
      " predictors."
    )
  }

  expected_index <- seq_len(
    nrow(predictor_map)
  )

  if (!identical(
    as.integer(predictor_map$predictor_index),
    expected_index
  )) {
    stop(
      map_name,
      "$predictor_index is not aligned with its row order."
    )
  }

  invisible(TRUE)
}


get_cs_snp_sets <- function(
    fit,
    predictor_map) {

  validate_predictor_map(
    fit,
    predictor_map,
    "predictor_map"
  )

  cs <- get_cs(fit)

  if (length(cs) == 0L) {
    return(list())
  }

  lapply(
    cs,
    function(index) {

      index <- as.integer(index)

      if (
        anyNA(index) ||
        any(index < 1L) ||
        any(index > nrow(predictor_map))
      ) {
        stop(
          "A credible set contains an invalid predictor index."
        )
      }

      clean_snp_vector(
        predictor_map$snp[index]
      )
    }
  )
}


get_cs_predictor_indices <- function(fit) {

  cs <- get_cs(fit)

  if (length(cs) == 0L) {
    return(integer(0))
  }

  index <- unique(
    as.integer(
      unlist(
        cs,
        use.names = FALSE
      )
    )
  )

  index[
    !is.na(index)
  ]
}


get_preferred_mix_coding <- function(
    fit,
    mix_predictor_map,
    snps,
    tie_tolerance = 1e-10) {

  if (length(snps) == 0L) {
    return(
      setNames(
        character(0),
        character(0)
      )
    )
  }

  if (
    is.null(fit$pip) ||
    length(fit$pip) != nrow(mix_predictor_map)
  ) {
    stop(
      "The mixed-model PIP vector does not match ",
      "the mixed predictor map."
    )
  }

  preferred_coding <- vapply(
    snps,
    function(snp) {

      index <- which(
        mix_predictor_map$snp == snp
      )

      if (length(index) == 0L) {
        return("unresolved")
      }

      snp_pip <- fit$pip[index]
      valid <- is.finite(snp_pip)

      if (!any(valid)) {
        return("unresolved")
      }

      index <- index[valid]
      snp_pip <- snp_pip[valid]

      maximum_pip <- max(snp_pip)

      winner <- index[
        snp_pip >= maximum_pip - tie_tolerance
      ]

      winner_coding <- unique(
        as.character(
          mix_predictor_map$coding[winner]
        )
      )

      winner_coding <- winner_coding[
        !is.na(winner_coding) &
          nzchar(winner_coding)
      ]

      if (length(winner_coding) == 1L) {
        return(winner_coding)
      }

      "ambiguous"
    },
    character(1L)
  )

  names(preferred_coding) <- snps

  preferred_coding
}


calculate_cs_overlap <- function(
    add_fit,
    mix_fit,
    add_predictor_map,
    mix_predictor_map) {

  validate_predictor_map(
    add_fit,
    add_predictor_map,
    "add_predictor_map"
  )

  validate_predictor_map(
    mix_fit,
    mix_predictor_map,
    "mix_predictor_map"
  )

  add_cs <- get_cs_snp_sets(
    add_fit,
    add_predictor_map
  )

  mix_cs <- get_cs_snp_sets(
    mix_fit,
    mix_predictor_map
  )

  # Pairwise biological-SNP overlap between credible sets.
  overlap_matrix <- matrix(
    0L,
    nrow = length(add_cs),
    ncol = length(mix_cs),
    dimnames = list(
      names(add_cs),
      names(mix_cs)
    )
  )

  if (
    length(add_cs) > 0L &&
    length(mix_cs) > 0L
  ) {

    for (i in seq_along(add_cs)) {
      for (j in seq_along(mix_cs)) {

        overlap_matrix[i, j] <- length(
          intersect(
            add_cs[[i]],
            mix_cs[[j]]
          )
        )
      }
    }
  }

  # Collect all predictor indices appearing in reported CSs.
  add_cs_index <- get_cs_predictor_indices(
    add_fit
  )

  mix_cs_index <- get_cs_predictor_indices(
    mix_fit
  )

  if (
    any(add_cs_index < 1L) ||
    any(add_cs_index > nrow(add_predictor_map))
  ) {
    stop(
      "The additive fit contains an invalid CS predictor index."
    )
  }

  if (
    any(mix_cs_index < 1L) ||
    any(mix_cs_index > nrow(mix_predictor_map))
  ) {
    stop(
      "The mixed fit contains an invalid CS predictor index."
    )
  }

  add_cs_members <- add_predictor_map[
    add_cs_index,
    ,
    drop = FALSE
  ]

  mix_cs_members <- mix_predictor_map[
    mix_cs_index,
    ,
    drop = FALSE
  ]

  # Distinct biological SNPs across all reported CSs.
  add_snps <- clean_snp_vector(
    add_cs_members$snp
  )

  mix_snps <- clean_snp_vector(
    mix_cs_members$snp
  )

  # Biological SNPs represented by each coding in mixed CSs.
  mix_additive_snps <- clean_snp_vector(
    mix_cs_members$snp[
      mix_cs_members$coding == "additive"
    ]
  )

  mix_nonadditive_snps <- clean_snp_vector(
    mix_cs_members$snp[
      mix_cs_members$coding %in%
        c("dominant", "recessive")
    ]
  )

  # SNP overlap ignores coding.
  shared_snps <- intersect(
    add_snps,
    mix_snps
  )

  add_snps_not_retained <- setdiff(
    add_snps,
    mix_snps
  )

  # Coding overlap requires the additive coding itself to occur
  # in a mixed-model credible set.
  additive_coding_overlap <- intersect(
    add_snps,
    mix_additive_snps
  )

  # Additive-CS SNPs retained in mixed CSs only through a
  # dominant or recessive predictor.
  other_coding_only <- intersect(
    add_snps,
    setdiff(
      mix_nonadditive_snps,
      mix_additive_snps
    )
  )

  # Among shared biological SNPs, identify the coding-specific
  # predictor with the largest mixed-model PIP.
  preferred_shared_coding <- get_preferred_mix_coding(
    fit = mix_fit,
    mix_predictor_map = mix_predictor_map,
    snps = shared_snps
  )

  n_preferred_additive <- sum(
    preferred_shared_coding == "additive"
  )

  n_preferred_dominant <- sum(
    preferred_shared_coding == "dominant"
  )

  n_preferred_recessive <- sum(
    preferred_shared_coding == "recessive"
  )

  n_preferred_nonadditive <- sum(
    preferred_shared_coding %in%
      c("dominant", "recessive")
  )

  n_preference_ambiguous <- sum(
    preferred_shared_coding %in%
      c("ambiguous", "unresolved")
  )

  list(
    overlap_matrix = overlap_matrix,

    n_add_cs_snps = length(add_snps),
    n_mix_cs_snps = length(mix_snps),

    # Biological-SNP overlap, irrespective of coding.
    overlap_snp = length(shared_snps),

    pct_add_cs_snps_retained = safe_percentage(
      length(shared_snps),
      length(add_snps)
    ),

    n_add_cs_snps_not_retained = length(
      add_snps_not_retained
    ),

    pct_add_cs_snps_not_retained = safe_percentage(
      length(add_snps_not_retained),
      length(add_snps)
    ),

    # Additive-coding overlap.
    overlap_coding = length(
      additive_coding_overlap
    ),

    pct_add_cs_snps_retained_additive_coding = (
      safe_percentage(
        length(additive_coding_overlap),
        length(add_snps)
      )
    ),

    # Additive-CS SNPs retained only under nonadditive coding.
    n_add_cs_snps_other_coding_only = length(
      other_coding_only
    ),

    pct_add_cs_snps_other_coding_only = (
      safe_percentage(
        length(other_coding_only),
        length(add_snps)
      )
    ),

    pct_shared_snps_other_coding_only = (
      safe_percentage(
        length(other_coding_only),
        length(shared_snps)
      )
    ),

    # Highest-PIP coding among shared SNPs.
    n_shared_snps_preferred_additive = (
      n_preferred_additive
    ),

    n_shared_snps_preferred_dominant = (
      n_preferred_dominant
    ),

    n_shared_snps_preferred_recessive = (
      n_preferred_recessive
    ),

    n_shared_snps_preferred_nonadditive = (
      n_preferred_nonadditive
    ),

    n_shared_snps_preference_ambiguous = (
      n_preference_ambiguous
    ),

    pct_shared_snps_preferred_nonadditive = (
      safe_percentage(
        n_preferred_nonadditive,
        length(shared_snps)
      )
    ),

    # Existing CS-pair overlap summaries.
    n_overlapping_cs_pairs = sum(
      overlap_matrix > 0L
    ),

    pairwise_overlap_sum = sum(
      overlap_matrix
    )
  )
}


get_cs_lead_predictor <- function(
    fit,
    cs,
    cs_number) {

  index <- as.integer(
    cs[[cs_number]]
  )

  if (length(index) == 0L) {
    return(NA_integer_)
  }

  cs_index <- fit$sets$cs_index

  if (
    !is.null(cs_index) &&
    length(cs_index) >= cs_number &&
    !is.na(cs_index[cs_number]) &&
    cs_index[cs_number] >= 1L &&
    cs_index[cs_number] <= nrow(fit$alpha)
  ) {

    probability <- fit$alpha[
      cs_index[cs_number],
      index
    ]

  } else {

    probability <- fit$pip[index]
  }

  valid <- is.finite(probability)

  if (!any(valid)) {
    return(index[1])
  }

  index <- index[valid]
  probability <- probability[valid]

  index[
    which.max(probability)
  ]
}


count_mix_cs_types <- function(
    fit,
    mix_predictor_map) {

  validate_predictor_map(
    fit,
    mix_predictor_map,
    "mix_predictor_map"
  )

  type_count <- c(
    additive = 0L,
    recessive = 0L,
    dominant = 0L
  )

  cs <- get_cs(fit)

  if (length(cs) == 0L) {
    return(type_count)
  }

  for (i in seq_along(cs)) {

    lead_index <- get_cs_lead_predictor(
      fit = fit,
      cs = cs,
      cs_number = i
    )

    if (is.na(lead_index)) {
      next
    }

    coding <- as.character(
      mix_predictor_map$coding[lead_index]
    )

    if (
      !is.na(coding) &&
      coding %in% names(type_count)
    ) {
      type_count[coding] <- (
        type_count[coding] + 1L
      )
    }
  }

  type_count
}


is_successful_tissue_result <- function(x) {

  if (!is.list(x)) {
    return(FALSE)
  }

  required_fields <- c(
    "susie_add",
    "susie_add_perm",
    "susie_mix",
    "susie_mix_perm",
    "add_predictor_map",
    "mix_predictor_map"
  )

  all(
    required_fields %in% names(x)
  )
}


summarize_tissue <- function(
    x,
    gene,
    tissue,
    result_file) {

  overlap_real <- calculate_cs_overlap(
    add_fit = x$susie_add,
    mix_fit = x$susie_mix,
    add_predictor_map = x$add_predictor_map,
    mix_predictor_map = x$mix_predictor_map
  )

  overlap_perm <- calculate_cs_overlap(
    add_fit = x$susie_add_perm,
    mix_fit = x$susie_mix_perm,
    add_predictor_map = x$add_predictor_map,
    mix_predictor_map = x$mix_predictor_map
  )

  type_real <- count_mix_cs_types(
    x$susie_mix,
    x$mix_predictor_map
  )

  type_perm <- count_mix_cs_types(
    x$susie_mix_perm,
    x$mix_predictor_map
  )

  data.frame(
    gene = gene,
    tissue = tissue,
    result_file = basename(result_file),

    # Sample and predictor counts.
    n_ind = as.numeric(
      value_or_na(x$n_ind)
    ),
    n_SNP = as.numeric(
      value_or_na(x$n_SNP)
    ),
    n_mix_predictor = as.numeric(
      value_or_na(x$n_mix_predictor)
    ),
    n_add_rm = as.numeric(
      value_or_na(x$n_add_rm)
    ),
    n_rec_rm = as.numeric(
      value_or_na(x$n_rec_rm)
    ),
    n_dom_rm = as.numeric(
      value_or_na(x$n_dom_rm)
    ),

    # Expression and marginal association summaries.
    mean_count = as.numeric(
      value_or_na(x$mean_read)
    ),
    median_count = as.numeric(
      value_or_na(x$median_read)
    ),
    min_pv = as.numeric(
      value_or_na(x$min_pv)
    ),

    # ELBO differences: positive values favor the mixed model.
    dif_elbo = (
      max_elbo(x$susie_mix) -
        max_elbo(x$susie_add)
    ),
    dif_elbo_perm = (
      max_elbo(x$susie_mix_perm) -
        max_elbo(x$susie_add_perm)
    ),

    # Previous log-likelihood-like summary.
    log_lik_add = get_log_lik_metric(
      x$susie_add
    ),
    log_lik_mix = get_log_lik_metric(
      x$susie_mix
    ),
    log_lik_add_perm = get_log_lik_metric(
      x$susie_add_perm
    ),
    log_lik_mix_perm = get_log_lik_metric(
      x$susie_mix_perm
    ),

    # Credible-set counts.
    ncs_susie = length(
      get_cs(x$susie_add)
    ),
    ncs_susie_mix = length(
      get_cs(x$susie_mix)
    ),
    perm_cs_susie = length(
      get_cs(x$susie_add_perm)
    ),
    perm_cs_susie_mix = length(
      get_cs(x$susie_mix_perm)
    ),

    # Biological-SNP overlap, ignoring coding.
    n_add_cs_snps = overlap_real$n_add_cs_snps,
    n_mix_cs_snps = overlap_real$n_mix_cs_snps,
    overlap_snp = overlap_real$overlap_snp,
    pct_add_cs_snps_retained = (
      overlap_real$pct_add_cs_snps_retained
    ),
    n_add_cs_snps_not_retained = (
      overlap_real$n_add_cs_snps_not_retained
    ),
    pct_add_cs_snps_not_retained = (
      overlap_real$pct_add_cs_snps_not_retained
    ),

    # Backward-compatible alias for older descriptive code.
    overlap = overlap_real$overlap_snp,

    overlap_cs_pairs = (
      overlap_real$n_overlapping_cs_pairs
    ),
    overlap_pairwise_sum = (
      overlap_real$pairwise_overlap_sum
    ),

    # Coding-aware overlap.
    overlap_coding = overlap_real$overlap_coding,
    pct_add_cs_snps_retained_additive_coding = (
      overlap_real$
        pct_add_cs_snps_retained_additive_coding
    ),
    n_add_cs_snps_other_coding_only = (
      overlap_real$n_add_cs_snps_other_coding_only
    ),
    pct_add_cs_snps_other_coding_only = (
      overlap_real$pct_add_cs_snps_other_coding_only
    ),
    pct_shared_snps_other_coding_only = (
      overlap_real$pct_shared_snps_other_coding_only
    ),

    # Highest-PIP coding among shared biological SNPs.
    n_shared_snps_preferred_additive = (
      overlap_real$n_shared_snps_preferred_additive
    ),
    n_shared_snps_preferred_dominant = (
      overlap_real$n_shared_snps_preferred_dominant
    ),
    n_shared_snps_preferred_recessive = (
      overlap_real$n_shared_snps_preferred_recessive
    ),
    n_shared_snps_preferred_nonadditive = (
      overlap_real$n_shared_snps_preferred_nonadditive
    ),
    n_shared_snps_preference_ambiguous = (
      overlap_real$n_shared_snps_preference_ambiguous
    ),
    pct_shared_snps_preferred_nonadditive = (
      overlap_real$pct_shared_snps_preferred_nonadditive
    ),

    # Permuted biological-SNP overlap.
    n_add_cs_snps_perm = overlap_perm$n_add_cs_snps,
    n_mix_cs_snps_perm = overlap_perm$n_mix_cs_snps,
    overlap_snp_perm = overlap_perm$overlap_snp,
    overlap_perm = overlap_perm$overlap_snp,
    pct_add_cs_snps_retained_perm = (
      overlap_perm$pct_add_cs_snps_retained
    ),
    n_add_cs_snps_not_retained_perm = (
      overlap_perm$n_add_cs_snps_not_retained
    ),
    pct_add_cs_snps_not_retained_perm = (
      overlap_perm$pct_add_cs_snps_not_retained
    ),
    overlap_cs_pairs_perm = (
      overlap_perm$n_overlapping_cs_pairs
    ),
    overlap_pairwise_sum_perm = (
      overlap_perm$pairwise_overlap_sum
    ),

    # Permuted coding-aware overlap and preference.
    overlap_coding_perm = overlap_perm$overlap_coding,
    pct_add_cs_snps_retained_additive_coding_perm = (
      overlap_perm$
        pct_add_cs_snps_retained_additive_coding
    ),
    n_add_cs_snps_other_coding_only_perm = (
      overlap_perm$n_add_cs_snps_other_coding_only
    ),
    pct_add_cs_snps_other_coding_only_perm = (
      overlap_perm$pct_add_cs_snps_other_coding_only
    ),
    n_shared_snps_preferred_additive_perm = (
      overlap_perm$n_shared_snps_preferred_additive
    ),
    n_shared_snps_preferred_dominant_perm = (
      overlap_perm$n_shared_snps_preferred_dominant
    ),
    n_shared_snps_preferred_recessive_perm = (
      overlap_perm$n_shared_snps_preferred_recessive
    ),
    n_shared_snps_preferred_nonadditive_perm = (
      overlap_perm$n_shared_snps_preferred_nonadditive
    ),
    n_shared_snps_preference_ambiguous_perm = (
      overlap_perm$n_shared_snps_preference_ambiguous
    ),
    pct_shared_snps_preferred_nonadditive_perm = (
      overlap_perm$pct_shared_snps_preferred_nonadditive
    ),

    # Coding of the lead predictor in each mixed-model CS.
    n_add = unname(
      type_real["additive"]
    ),
    n_rec = unname(
      type_real["recessive"]
    ),
    n_dom = unname(
      type_real["dominant"]
    ),
    n_add_perm = unname(
      type_perm["additive"]
    ),
    n_rec_perm = unname(
      type_perm["recessive"]
    ),
    n_dom_perm = unname(
      type_perm["dominant"]
    ),

    # Convergence diagnostics.
    converged_add = isTRUE(
      x$susie_add$converged
    ),
    converged_mix = isTRUE(
      x$susie_mix$converged
    ),
    converged_add_perm = isTRUE(
      x$susie_add_perm$converged
    ),
    converged_mix_perm = isTRUE(
      x$susie_mix_perm$converged
    ),

    stringsAsFactors = FALSE
  )
}


# ============================================================
# Read and summarize result files
# ============================================================

result_rows <- list()
error_rows <- list()

result_counter <- 0L
error_counter <- 0L


record_error <- function(
    file,
    gene,
    stage,
    message) {

  error_counter <<- error_counter + 1L

  error_rows[[error_counter]] <<- data.frame(
    result_file = basename(file),
    gene = gene,
    stage = stage,
    error = as.character(message),
    stringsAsFactors = FALSE
  )
}


for (file_index in seq_along(result_files)) {

  result_file <- result_files[file_index]

  gene_from_file <- tools::file_path_sans_ext(
    basename(result_file)
  )

  out <- tryCatch(
    readRDS(result_file),
    error = function(e) e
  )

  if (inherits(out, "error")) {

    record_error(
      file = result_file,
      gene = gene_from_file,
      stage = "readRDS",
      message = conditionMessage(out)
    )

    next
  }

  # Explicit job-failure object:
  # list(gene = "...", error = "...")
  if (
    is.list(out) &&
    all(c("gene", "error") %in% names(out))
  ) {

    failed_gene <- as.character(
      value_or_na(out$gene)
    )

    if (
      is.na(failed_gene) ||
      !nzchar(failed_gene)
    ) {
      failed_gene <- gene_from_file
    }

    record_error(
      file = result_file,
      gene = failed_gene,
      stage = "gene_analysis",
      message = value_or_na(out$error)
    )

    next
  }

  if (!is.list(out)) {

    record_error(
      file = result_file,
      gene = gene_from_file,
      stage = "invalid_result",
      message = paste(
        "Expected a list of tissue results, but found",
        paste(class(out), collapse = "/")
      )
    )

    next
  }

  if (length(out) == 0L) {

    record_error(
      file = result_file,
      gene = gene_from_file,
      stage = "no_tissues",
      message = paste(
        "No tissues passed the sample and predictor filters."
      )
    )

    next
  }

  tissue_names <- names(out)

  if (is.null(tissue_names)) {
    tissue_names <- paste0(
      "unnamed_tissue_",
      seq_along(out)
    )
  }

  for (tissue_index in seq_along(out)) {

    tissue_result <- out[[tissue_index]]
    tissue <- tissue_names[tissue_index]

    if (
      is.na(tissue) ||
      !nzchar(tissue)
    ) {
      tissue <- paste0(
        "unnamed_tissue_",
        tissue_index
      )
    }

    if (!is_successful_tissue_result(
      tissue_result
    )) {

      record_error(
        file = result_file,
        gene = gene_from_file,
        stage = "invalid_tissue_result",
        message = paste(
          "Tissue",
          tissue,
          "does not contain the expected SuSiE fits and predictor maps."
        )
      )

      next
    }

    temp_row <- tryCatch(
      summarize_tissue(
        x = tissue_result,
        gene = gene_from_file,
        tissue = tissue,
        result_file = result_file
      ),
      error = function(e) e
    )

    if (inherits(temp_row, "error")) {

      record_error(
        file = result_file,
        gene = gene_from_file,
        stage = paste0(
          "summarize_tissue:",
          tissue
        ),
        message = conditionMessage(temp_row)
      )

      next
    }

    result_counter <- result_counter + 1L
    result_rows[[result_counter]] <- temp_row
  }

  if (
    file_index %% 100L == 0L ||
    file_index == length(result_files)
  ) {
    cat(
      "Processed",
      file_index,
      "of",
      length(result_files),
      "files;",
      result_counter,
      "tissue results;",
      error_counter,
      "errors.\n"
    )
  }
}


# ============================================================
# Combine and save
# ============================================================

if (length(result_rows) > 0L) {

  res_summary <- do.call(
    rbind,
    result_rows
  )

  rownames(res_summary) <- NULL

} else {

  res_summary <- data.frame()
}


if (length(error_rows) > 0L) {

  res_errors <- do.call(
    rbind,
    error_rows
  )

  rownames(res_errors) <- NULL

} else {

  res_errors <- data.frame(
    result_file = character(),
    gene = character(),
    stage = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
}


save(
  res_summary,
  res_errors,
  file = summary_file
)

write.csv(
  res_summary,
  summary_csv,
  row.names = FALSE
)

write.csv(
  res_errors,
  error_csv,
  row.names = FALSE
)


failed_genes <- sort(
  unique(
    res_errors$gene[
      !is.na(res_errors$gene) &
        nzchar(res_errors$gene)
    ]
  )
)

writeLines(
  failed_genes,
  failed_gene_file
)


cat("\nSummary complete.\n")
cat(
  "Successful genes:",
  length(unique(res_summary$gene)),
  "\n"
)
cat(
  "Successful gene-tissue analyses:",
  nrow(res_summary),
  "\n"
)
cat(
  "Files or tissues with errors:",
  nrow(res_errors),
  "\n"
)
cat(
  "Genes to inspect or rerun:",
  length(failed_genes),
  "\n"
)
cat("Saved:", summary_file, "\n")
cat("Saved:", summary_csv, "\n")
cat("Saved:", error_csv, "\n")
cat("Saved:", failed_gene_file, "\n")
