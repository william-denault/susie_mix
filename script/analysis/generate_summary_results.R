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

cs_summary_file <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "res_cs_summary.RData"
)

cs_summary_csv <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "res_cs_summary.csv"
)

cs_error_csv <- paste0(
  "/project2/mstephens/wdenault/susie_mix/",
  "res_cs_errors.csv"
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


find_metadata_scalar <- function(
    objects,
    candidate_names,
    numeric = FALSE) {

  for (object in objects) {

    if (
      is.null(object) ||
      is.null(names(object))
    ) {
      next
    }

    for (candidate_name in candidate_names) {

      normalized_object_names <- tolower(
        gsub(
          "[^a-z0-9]",
          "",
          names(object)
        )
      )

      normalized_candidate_name <- tolower(
        gsub(
          "[^a-z0-9]",
          "",
          candidate_name
        )
      )

      matched_name_index <- match(
        normalized_candidate_name,
        normalized_object_names
      )

      if (is.na(matched_name_index)) {
        next
      }

      matched_name <- names(object)[matched_name_index]
      value <- object[[matched_name]]

      if (numeric) {

        value <- suppressWarnings(
          as.numeric(
            as.character(value)
          )
        )

        value <- value[is.finite(value)]

      } else {

        value <- as.character(value)
        value <- value[
          !is.na(value) &
            nzchar(value)
        ]
      }

      if (length(value) > 0L) {
        return(list(
          value = value[1],
          source = matched_name
        ))
      }
    }
  }

  list(
    value = if (numeric) NA_real_ else NA_character_,
    source = NA_character_
  )
}


get_tss_metadata <- function(
    x,
    predictor_map) {

  nested_metadata_names <- intersect(
    c(
      "gene_info",
      "gene_metadata",
      "gene_annotation",
      "annotation"
    ),
    names(x)
  )

  metadata_objects <- c(
    list(x),
    unname(x[nested_metadata_names]),
    list(predictor_map)
  )

  strand_info <- find_metadata_scalar(
    metadata_objects,
    c(
      "strand",
      "gene_strand",
      "tx_strand"
    )
  )

  strand <- strand_info$value

  if (!is.na(strand)) {

    strand_lower <- tolower(strand)

    if (strand_lower %in% c("-", "-1", "minus", "reverse")) {
      strand <- "-"
    } else if (
      strand_lower %in% c("+", "1", "plus", "forward")
    ) {
      strand <- "+"
    } else {
      strand <- NA_character_
    }
  }

  tss_info <- find_metadata_scalar(
    metadata_objects,
    c(
      "tss",
      "TSS",
      "tss_position",
      "tss_pos",
      "gene_tss",
      "transcription_start_site"
    ),
    numeric = TRUE
  )

  if (!is.finite(tss_info$value)) {

    gene_start_info <- find_metadata_scalar(
      metadata_objects,
      c(
        "gene_start",
        "tx_start",
        "transcription_start"
      ),
      numeric = TRUE
    )

    gene_end_info <- find_metadata_scalar(
      metadata_objects,
      c(
        "gene_end",
        "tx_end",
        "transcription_end"
      ),
      numeric = TRUE
    )

    if (
      identical(strand, "-") &&
      is.finite(gene_end_info$value)
    ) {

      tss_info <- list(
        value = gene_end_info$value,
        source = paste0(
          "derived_from_",
          gene_end_info$source
        )
      )

    } else if (is.finite(gene_start_info$value)) {

      tss_info <- list(
        value = gene_start_info$value,
        source = paste0(
          "derived_from_",
          gene_start_info$source
        )
      )
    }
  }

  chromosome_info <- find_metadata_scalar(
    metadata_objects,
    c(
      "chromosome",
      "chrom",
      "chr",
      "gene_chromosome"
    )
  )

  list(
    chromosome = chromosome_info$value,
    tss_position = as.numeric(tss_info$value),
    strand = strand,
    tss_source = tss_info$source
  )
}


parse_variant_position <- function(snp) {

  snp <- as.character(snp)

  matches <- regexec(
    paste0(
      "^(?:chr)?(?:[0-9]+|X|Y|M|MT)",
      "[:_]([0-9]+)"
    ),
    snp,
    ignore.case = TRUE,
    perl = TRUE
  )

  match_values <- regmatches(
    snp,
    matches
  )

  vapply(
    match_values,
    function(value) {

      if (length(value) < 2L) {
        return(NA_real_)
      }

      suppressWarnings(
        as.numeric(value[2])
      )
    },
    numeric(1)
  )
}


parse_variant_chromosome <- function(snp) {

  snp <- as.character(snp)

  matches <- regexec(
    "^((?:chr)?(?:[0-9]+|X|Y|M|MT))[:_]",
    snp,
    ignore.case = TRUE,
    perl = TRUE
  )

  match_values <- regmatches(
    snp,
    matches
  )

  vapply(
    match_values,
    function(value) {

      if (length(value) < 2L) {
        return(NA_character_)
      }

      value[2]
    },
    character(1)
  )
}


get_predictor_coordinates <- function(predictor_map) {

  position_candidates <- c(
    "position",
    "pos",
    "bp",
    "snp_position",
    "variant_position"
  )

  position_column <- intersect(
    position_candidates,
    names(predictor_map)
  )

  position <- rep(
    NA_real_,
    nrow(predictor_map)
  )

  position_source <- rep(
    NA_character_,
    nrow(predictor_map)
  )

  if (length(position_column) > 0L) {

    position_column <- position_column[1]

    position <- suppressWarnings(
      as.numeric(
        as.character(
          predictor_map[[position_column]]
        )
      )
    )

    position_source[is.finite(position)] <- position_column
  }

  parsed_position <- parse_variant_position(
    predictor_map$snp
  )

  use_parsed_position <- (
    !is.finite(position) &
      is.finite(parsed_position)
  )

  position[use_parsed_position] <- (
    parsed_position[use_parsed_position]
  )

  position_source[use_parsed_position] <- "parsed_from_snp"

  chromosome_candidates <- c(
    "chromosome",
    "chrom",
    "chr"
  )

  chromosome_column <- intersect(
    chromosome_candidates,
    names(predictor_map)
  )

  chromosome <- rep(
    NA_character_,
    nrow(predictor_map)
  )

  if (length(chromosome_column) > 0L) {

    chromosome_column <- chromosome_column[1]
    chromosome <- as.character(
      predictor_map[[chromosome_column]]
    )
  }

  parsed_chromosome <- parse_variant_chromosome(
    predictor_map$snp
  )

  use_parsed_chromosome <- (
    (is.na(chromosome) | !nzchar(chromosome)) &
      !is.na(parsed_chromosome)
  )

  chromosome[use_parsed_chromosome] <- (
    parsed_chromosome[use_parsed_chromosome]
  )

  list(
    position = position,
    position_source = position_source,
    chromosome = chromosome
  )
}


get_record_scalar <- function(
    record,
    candidate_names,
    numeric = FALSE) {

  if (
    is.null(record) ||
    length(record) == 0L
  ) {
    return(list(
      value = if (numeric) NA_real_ else NA_character_,
      source = NA_character_
    ))
  }

  find_metadata_scalar(
    objects = list(record),
    candidate_names = candidate_names,
    numeric = numeric
  )
}


get_saved_cs_tss_record <- function(
    saved_tss_summary,
    cs_number,
    cs_name,
    component) {

  if (
    is.null(saved_tss_summary) ||
    length(saved_tss_summary) == 0L
  ) {
    return(list())
  }

  if (
    is.atomic(saved_tss_summary) &&
    is.null(dim(saved_tss_summary))
  ) {

    if (length(saved_tss_summary) < cs_number) {
      return(list())
    }

    return(list(
      distance_to_tss = saved_tss_summary[cs_number]
    ))
  }

  if (
    is.list(saved_tss_summary) &&
    !is.data.frame(saved_tss_summary)
  ) {

    field_names <- c(
      "cs",
      "cs_number",
      "cs_name",
      "component",
      "lead_snp",
      "distance",
      "distance_to_tss",
      "tss_distance"
    )

    if (!any(field_names %in% names(saved_tss_summary))) {

      if (
        !is.null(names(saved_tss_summary)) &&
        cs_name %in% names(saved_tss_summary)
      ) {
        return(saved_tss_summary[[cs_name]])
      }

      if (length(saved_tss_summary) >= cs_number) {
        return(saved_tss_summary[[cs_number]])
      }

      return(list())
    }
  }

  saved_tss_summary <- tryCatch(
    as.data.frame(
      saved_tss_summary,
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )

  if (
    is.null(saved_tss_summary) ||
    nrow(saved_tss_summary) == 0L
  ) {
    return(list())
  }

  selected_row <- NA_integer_

  name_columns <- intersect(
    c(
      "cs_name",
      "cs",
      "credible_set"
    ),
    names(saved_tss_summary)
  )

  for (name_column in name_columns) {

    values <- as.character(
      saved_tss_summary[[name_column]]
    )

    matched_row <- which(
      values == cs_name
    )

    if (length(matched_row) > 0L) {
      selected_row <- matched_row[1]
      break
    }
  }

  if (is.na(selected_row)) {

    number_columns <- intersect(
      c(
        "cs_number",
        "cs_index",
        "component",
        "L"
      ),
      names(saved_tss_summary)
    )

    for (number_column in number_columns) {

      values <- suppressWarnings(
        as.numeric(
          as.character(
            saved_tss_summary[[number_column]]
          )
        )
      )

      target_value <- if (
        number_column %in% c("component", "L") &&
        is.finite(component)
      ) {
        component
      } else {
        cs_number
      }

      matched_row <- which(
        values == target_value
      )

      if (length(matched_row) > 0L) {
        selected_row <- matched_row[1]
        break
      }
    }
  }

  if (
    is.na(selected_row) &&
    nrow(saved_tss_summary) >= cs_number
  ) {
    selected_row <- cs_number
  }

  if (is.na(selected_row)) {
    return(list())
  }

  as.list(
    saved_tss_summary[
      selected_row,
      ,
      drop = FALSE
    ]
  )
}


get_cs_vector_value <- function(
    values,
    cs_number,
    cs_name) {

  if (
    is.null(values) ||
    length(values) == 0L
  ) {
    return(NA_real_)
  }

  if (
    !is.null(names(values)) &&
    cs_name %in% names(values)
  ) {
    return(as.numeric(values[[cs_name]]))
  }

  if (length(values) < cs_number) {
    return(NA_real_)
  }

  as.numeric(values[[cs_number]])
}


get_cs_purity_value <- function(
    fit,
    cs_number,
    cs_name,
    candidate_columns) {

  purity <- fit$sets$purity

  if (
    is.null(purity) ||
    length(purity) == 0L
  ) {
    return(NA_real_)
  }

  purity <- as.data.frame(
    purity,
    stringsAsFactors = FALSE
  )

  purity_column <- intersect(
    candidate_columns,
    names(purity)
  )

  if (length(purity_column) == 0L) {
    return(NA_real_)
  }

  purity_row <- cs_number

  if (
    !is.null(rownames(purity)) &&
    cs_name %in% rownames(purity)
  ) {
    purity_row <- match(
      cs_name,
      rownames(purity)
    )
  }

  if (
    is.na(purity_row) ||
    purity_row > nrow(purity)
  ) {
    return(NA_real_)
  }

  as.numeric(
    purity[
      purity_row,
      purity_column[1]
    ]
  )
}


summarize_credible_sets <- function(
    x,
    gene,
    tissue,
    result_file,
    tissue_summary) {

  model_inputs <- list(
    list(
      model = "SuSiE",
      model_key = "susie_add",
      fit = x$susie_add,
      predictor_map = x$add_predictor_map,
      saved_tss_summary = (
        x$susie_add_lead_snp_tss_distance
      )
    ),
    list(
      model = "SuSiE-mix",
      model_key = "susie_mix",
      fit = x$susie_mix,
      predictor_map = x$mix_predictor_map,
      saved_tss_summary = (
        x$susie_mix_lead_snp_tss_distance
      )
    )
  )

  cs_rows <- list()
  cs_row_counter <- 0L

  repeated_columns <- setdiff(
    names(tissue_summary),
    c(
      "gene",
      "tissue",
      "result_file"
    )
  )

  repeated_summary <- tissue_summary[
    1,
    repeated_columns,
    drop = FALSE
  ]

  for (model_input in model_inputs) {

    fit <- model_input$fit
    predictor_map <- model_input$predictor_map

    validate_predictor_map(
      fit,
      predictor_map,
      paste0(
        model_input$model_key,
        "_predictor_map"
      )
    )

    cs <- get_cs(fit)

    if (length(cs) == 0L) {
      next
    }

    coordinates <- get_predictor_coordinates(
      predictor_map
    )

    tss_metadata <- get_tss_metadata(
      x,
      predictor_map
    )

    cs_names <- names(cs)

    for (cs_number in seq_along(cs)) {

      predictor_index <- as.integer(
        cs[[cs_number]]
      )

      if (
        length(predictor_index) == 0L ||
        anyNA(predictor_index) ||
        any(predictor_index < 1L) ||
        any(predictor_index > nrow(predictor_map))
      ) {
        stop(
          "A credible set contains invalid predictor indices."
        )
      }

      cs_name <- paste0(
        "CS",
        cs_number
      )

      if (
        !is.null(cs_names) &&
        length(cs_names) >= cs_number &&
        !is.na(cs_names[cs_number]) &&
        nzchar(cs_names[cs_number])
      ) {
        cs_name <- cs_names[cs_number]
      }

      lead_index <- get_cs_lead_predictor(
        fit = fit,
        cs = cs,
        cs_number = cs_number
      )

      component <- get_cs_vector_value(
        fit$sets$cs_index,
        cs_number,
        cs_name
      )

      saved_cs_tss_record <- get_saved_cs_tss_record(
        saved_tss_summary = model_input$saved_tss_summary,
        cs_number = cs_number,
        cs_name = cs_name,
        component = component
      )

      saved_lead_snp_info <- get_record_scalar(
        saved_cs_tss_record,
        c(
          "lead_snp",
          "snp",
          "variant",
          "lead_variant"
        )
      )

      saved_lead_position_info <- get_record_scalar(
        saved_cs_tss_record,
        c(
          "lead_position",
          "lead_snp_position",
          "snp_position",
          "variant_position",
          "position",
          "pos",
          "bp"
        ),
        numeric = TRUE
      )

      saved_tss_position_info <- get_record_scalar(
        saved_cs_tss_record,
        c(
          "tss_position",
          "tss_pos",
          "tss",
          "TSS"
        ),
        numeric = TRUE
      )

      saved_distance_info <- get_record_scalar(
        saved_cs_tss_record,
        c(
          "distance_to_tss_bp",
          "distance_to_tss",
          "distance_tss",
          "tss_distance",
          "distance"
        ),
        numeric = TRUE
      )

      saved_distance_kb_info <- get_record_scalar(
        saved_cs_tss_record,
        c(
          "distance_to_tss_kb",
          "tss_distance_kb"
        ),
        numeric = TRUE
      )

      if (
        !is.finite(saved_distance_info$value) &&
        is.finite(saved_distance_kb_info$value)
      ) {
        saved_distance_info <- list(
          value = 1000 * saved_distance_kb_info$value,
          source = saved_distance_kb_info$source
        )
      }

      lead_alpha <- NA_real_

      if (
        is.finite(component) &&
        component >= 1L &&
        component <= nrow(fit$alpha) &&
        !is.na(lead_index)
      ) {
        lead_alpha <- fit$alpha[
          as.integer(component),
          lead_index
        ]
      }

      lead_pip <- NA_real_

      if (
        !is.null(fit$pip) &&
        !is.na(lead_index) &&
        length(fit$pip) >= lead_index
      ) {
        lead_pip <- fit$pip[lead_index]
      }

      member_snps <- clean_snp_vector(
        predictor_map$snp[predictor_index]
      )

      member_codings <- sort(
        clean_snp_vector(
          predictor_map$coding[predictor_index]
        )
      )

      member_positions <- unique(
        coordinates$position[predictor_index]
      )

      member_positions <- member_positions[
        is.finite(member_positions)
      ]

      saved_chromosome_info <- get_record_scalar(
        saved_cs_tss_record,
        c(
          "lead_chromosome",
          "chromosome",
          "chrom",
          "chr"
        )
      )

      saved_strand_info <- get_record_scalar(
        saved_cs_tss_record,
        c(
          "gene_strand",
          "strand"
        )
      )

      lead_snp <- as.character(
        predictor_map$snp[lead_index]
      )

      tss_lead_snp <- lead_snp

      if (!is.na(saved_lead_snp_info$value)) {
        tss_lead_snp <- saved_lead_snp_info$value
      }

      lead_position <- coordinates$position[lead_index]
      lead_position_source <- coordinates$position_source[lead_index]

      if (
        !is.finite(lead_position) &&
        is.finite(saved_lead_position_info$value)
      ) {
        lead_position <- saved_lead_position_info$value
        lead_position_source <- paste0(
          "workhorse:",
          saved_lead_position_info$source
        )
      }

      if (!is.finite(lead_position)) {

        parsed_lead_position <- parse_variant_position(
          tss_lead_snp
        )

        if (is.finite(parsed_lead_position)) {
          lead_position <- parsed_lead_position
          lead_position_source <- "parsed_from_workhorse_lead_snp"
        }
      }

      lead_chromosome <- coordinates$chromosome[lead_index]

      if (
        (is.na(lead_chromosome) || !nzchar(lead_chromosome)) &&
        !is.na(saved_chromosome_info$value)
      ) {
        lead_chromosome <- saved_chromosome_info$value
      }

      if (is.na(lead_chromosome) || !nzchar(lead_chromosome)) {
        lead_chromosome <- parse_variant_chromosome(
          tss_lead_snp
        )
      }

      if (is.finite(saved_tss_position_info$value)) {
        tss_metadata$tss_position <- saved_tss_position_info$value
        tss_metadata$tss_source <- paste0(
          "workhorse:",
          saved_tss_position_info$source
        )
      }

      if (!is.na(saved_strand_info$value)) {

        saved_strand <- tolower(
          saved_strand_info$value
        )

        if (saved_strand %in% c("-", "-1", "minus", "reverse")) {
          tss_metadata$strand <- "-"
        } else if (
          saved_strand %in% c("+", "1", "plus", "forward")
        ) {
          tss_metadata$strand <- "+"
        }
      }

      tss_position <- tss_metadata$tss_position
      genomic_offset <- lead_position - tss_position

      orientation_multiplier <- 1
      distance_orientation <- "genomic_coordinate"

      if (identical(tss_metadata$strand, "-")) {
        orientation_multiplier <- -1
        distance_orientation <- "transcription_direction"
      } else if (identical(tss_metadata$strand, "+")) {
        distance_orientation <- "transcription_direction"
      }

      distance_to_tss_bp <- (
        genomic_offset * orientation_multiplier
      )

      distance_source <- "reconstructed_from_lead_position_and_tss"

      if (
        !is.finite(distance_to_tss_bp) &&
        is.finite(saved_distance_info$value)
      ) {
        distance_to_tss_bp <- saved_distance_info$value
        distance_orientation <- "workhorse_saved_unspecified"
        distance_source <- paste0(
          "workhorse:",
          saved_distance_info$source
        )
      }

      member_distance_to_tss_bp <- (
        (member_positions - tss_position) *
          orientation_multiplier
      )

      chromosome_matches_tss <- NA

      if (
        !is.na(lead_chromosome) &&
        nzchar(lead_chromosome) &&
        !is.na(tss_metadata$chromosome) &&
        nzchar(tss_metadata$chromosome)
      ) {

        normalize_chromosome <- function(value) {
          tolower(
            sub(
              "^chr",
              "",
              value,
              ignore.case = TRUE
            )
          )
        }

        chromosome_matches_tss <- identical(
          normalize_chromosome(lead_chromosome),
          normalize_chromosome(tss_metadata$chromosome)
        )
      }

      cs_row <- data.frame(
        gene = gene,
        tissue = tissue,
        result_file = basename(result_file),
        model = model_input$model,
        model_key = model_input$model_key,
        cs_number = cs_number,
        cs_name = cs_name,
        cs_id = paste(
          gene,
          tissue,
          model_input$model_key,
          cs_number,
          sep = "|"
        ),
        component = component,
        cs_size_predictors = length(predictor_index),
        cs_size_snps = length(member_snps),
        cs_predictor_indices = paste(
          predictor_index,
          collapse = ";"
        ),
        cs_predictors = paste(
          predictor_map$predictor_name[predictor_index],
          collapse = ";"
        ),
        cs_snps = paste(
          member_snps,
          collapse = ";"
        ),
        cs_coding_types = paste(
          member_codings,
          collapse = ";"
        ),
        cs_has_multiple_codings = length(member_codings) > 1L,
        n_additive_predictors = sum(
          predictor_map$coding[predictor_index] == "additive",
          na.rm = TRUE
        ),
        n_recessive_predictors = sum(
          predictor_map$coding[predictor_index] == "recessive",
          na.rm = TRUE
        ),
        n_dominant_predictors = sum(
          predictor_map$coding[predictor_index] == "dominant",
          na.rm = TRUE
        ),
        cs_coverage = get_cs_vector_value(
          fit$sets$coverage,
          cs_number,
          cs_name
        ),
        cs_min_abs_corr = get_cs_purity_value(
          fit,
          cs_number,
          cs_name,
          c(
            "min.abs.corr",
            "min_abs_corr"
          )
        ),
        cs_mean_abs_corr = get_cs_purity_value(
          fit,
          cs_number,
          cs_name,
          c(
            "mean.abs.corr",
            "mean_abs_corr"
          )
        ),
        cs_median_abs_corr = get_cs_purity_value(
          fit,
          cs_number,
          cs_name,
          c(
            "median.abs.corr",
            "median_abs_corr"
          )
        ),
        lead_predictor_index = lead_index,
        lead_predictor = as.character(
          predictor_map$predictor_name[lead_index]
        ),
        lead_snp = lead_snp,
        workhorse_lead_snp = saved_lead_snp_info$value,
        tss_lead_snp = tss_lead_snp,
        lead_coding = as.character(
          predictor_map$coding[lead_index]
        ),
        lead_alpha = as.numeric(lead_alpha),
        lead_pip = as.numeric(lead_pip),
        lead_chromosome = lead_chromosome,
        lead_position = lead_position,
        lead_position_source = lead_position_source,
        gene_chromosome = tss_metadata$chromosome,
        tss_position = tss_position,
        gene_strand = tss_metadata$strand,
        tss_source = tss_metadata$tss_source,
        chromosome_matches_tss = chromosome_matches_tss,
        genomic_offset_from_tss_bp = genomic_offset,
        distance_to_tss_bp = distance_to_tss_bp,
        distance_to_tss_kb = distance_to_tss_bp / 1000,
        absolute_distance_to_tss_bp = abs(distance_to_tss_bp),
        distance_orientation = distance_orientation,
        distance_source = distance_source,
        workhorse_distance_to_tss_bp = saved_distance_info$value,
        workhorse_distance_source = saved_distance_info$source,
        nearest_cs_snp_abs_distance_to_tss_bp = if (
          length(member_distance_to_tss_bp) > 0L
        ) {
          min(abs(member_distance_to_tss_bp))
        } else {
          NA_real_
        },
        min_cs_snp_distance_to_tss_bp = if (
          length(member_distance_to_tss_bp) > 0L
        ) {
          min(member_distance_to_tss_bp)
        } else {
          NA_real_
        },
        max_cs_snp_distance_to_tss_bp = if (
          length(member_distance_to_tss_bp) > 0L
        ) {
          max(member_distance_to_tss_bp)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      cs_row_counter <- cs_row_counter + 1L
      cs_rows[[cs_row_counter]] <- cbind(
        cs_row,
        repeated_summary
      )
    }
  }

  if (length(cs_rows) == 0L) {
    return(data.frame())
  }

  cs_summary <- do.call(
    rbind,
    cs_rows
  )

  rownames(cs_summary) <- NULL
  cs_summary
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
cs_result_rows <- list()
error_rows <- list()
cs_error_rows <- list()

result_counter <- 0L
cs_result_counter <- 0L
cs_row_counter <- 0L
error_counter <- 0L
cs_error_counter <- 0L


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


record_cs_error <- function(
    file,
    gene,
    tissue,
    message) {

  cs_error_counter <<- cs_error_counter + 1L

  cs_error_rows[[cs_error_counter]] <<- data.frame(
    result_file = basename(file),
    gene = gene,
    tissue = tissue,
    stage = "summarize_credible_sets",
    error = as.character(message),
    stringsAsFactors = FALSE
  )
}


for (file_index in  seq_along(result_files)) {

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

    temp_cs_rows <- tryCatch(
      summarize_credible_sets(
        x = tissue_result,
        gene = gene_from_file,
        tissue = tissue,
        result_file = result_file,
        tissue_summary = temp_row
      ),
      error = function(e) e
    )

    if (inherits(temp_cs_rows, "error")) {

      record_cs_error(
        file = result_file,
        gene = gene_from_file,
        tissue = tissue,
        message = conditionMessage(temp_cs_rows)
      )

    } else if (nrow(temp_cs_rows) > 0L) {

      cs_result_counter <- cs_result_counter + 1L
      cs_result_rows[[cs_result_counter]] <- temp_cs_rows
      cs_row_counter <- cs_row_counter + nrow(temp_cs_rows)
    }
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
      cs_row_counter,
      "credible sets;",
      error_counter,
      "original-summary errors;",
      cs_error_counter,
      "CS-summary errors.\n"
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


if (length(cs_result_rows) > 0L) {

  res_cs_summary <- do.call(
    rbind,
    cs_result_rows
  )

  rownames(res_cs_summary) <- NULL

} else {

  res_cs_summary <- data.frame()
}


if (length(cs_error_rows) > 0L) {

  cs_errors <- do.call(
    rbind,
    cs_error_rows
  )

  rownames(cs_errors) <- NULL

} else {

  cs_errors <- data.frame(
    result_file = character(),
    gene = character(),
    tissue = character(),
    stage = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
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

save(
  res_cs_summary,
  cs_errors,
  file = cs_summary_file
)

write.csv(
  res_summary,
  summary_csv,
  row.names = FALSE
)

write.csv(
  res_cs_summary,
  cs_summary_csv,
  row.names = FALSE
)

write.csv(
  cs_errors,
  cs_error_csv,
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
  "Credible-set summary rows:",
  nrow(res_cs_summary),
  "\n"
)
if (
  nrow(res_cs_summary) > 0L &&
  "distance_to_tss_bp" %in% names(res_cs_summary)
) {
  cat(
    "CS rows with a TSS distance:",
    sum(
      is.finite(res_cs_summary$distance_to_tss_bp)
    ),
    "\n"
  )
}
cat(
  "Files or tissues with errors:",
  nrow(res_errors),
  "\n"
)
cat(
  "CS-summary errors:",
  nrow(cs_errors),
  "\n"
)
cat(
  "Genes to inspect or rerun:",
  length(failed_genes),
  "\n"
)
cat("Saved:", summary_file, "\n")
cat("Saved:", summary_csv, "\n")
cat("Saved:", cs_summary_file, "\n")
cat("Saved:", cs_summary_csv, "\n")
cat("Saved:", cs_error_csv, "\n")
cat("Saved:", error_csv, "\n")
cat("Saved:", failed_gene_file, "\n")
