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

      snps <- predictor_map$snp[index]

      unique(
        snps[
          !is.na(snps)
        ]
      )
    }
  )
}


calculate_cs_overlap <- function(
    add_fit,
    mix_fit,
    add_predictor_map,
    mix_predictor_map) {

  add_cs <- get_cs_snp_sets(
    add_fit,
    add_predictor_map
  )

  mix_cs <- get_cs_snp_sets(
    mix_fit,
    mix_predictor_map
  )

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

  add_union <- unique(
    unlist(
      add_cs,
      use.names = FALSE
    )
  )

  mix_union <- unique(
    unlist(
      mix_cs,
      use.names = FALSE
    )
  )

  list(
    overlap_matrix = overlap_matrix,

    # Number of distinct biological SNPs in both analyses.
    n_shared_snps = length(
      intersect(
        add_union,
        mix_union
      )
    ),

    # Number of CS pairs sharing at least one biological SNP.
    n_overlapping_cs_pairs = sum(
      overlap_matrix > 0L
    ),

    # Total pairwise number of shared SNPs.
    pairwise_overlap_sum = sum(
      overlap_matrix
    )
  )
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

  for (index in cs) {

    if (length(index) == 0L) {
      next
    }

    # Preserve the original definition: classify each CS using
    # its first, highest-priority predictor.
    lead_index <- as.integer(index[1])

    coding <- mix_predictor_map$coding[
      lead_index
    ]

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

    # SNP-level overlap.
    overlap = overlap_real$n_shared_snps,
    overlap_cs_pairs = (
      overlap_real$n_overlapping_cs_pairs
    ),
    overlap_pairwise_sum = (
      overlap_real$pairwise_overlap_sum
    ),

    # Permuted SNP-level overlap.
    overlap_perm = overlap_perm$n_shared_snps,
    overlap_cs_pairs_perm = (
      overlap_perm$n_overlapping_cs_pairs
    ),
    overlap_pairwise_sum_perm = (
      overlap_perm$pairwise_overlap_sum
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

  # The RDS file itself could not be read.
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

  # Successful workhorse output must be a list of tissues.
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

    if (
      !is_successful_tissue_result(
        tissue_result
      )
    ) {

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
cat("Saved:", error_csv, "\n")
cat("Saved:", failed_gene_file, "\n")


load(summary_file)
