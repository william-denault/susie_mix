# Shared utilities for plotting normalized GTEx expression against the lead
# biological SNP from a saved SuSiE-mix credible set. Data reconstruction
# intentionally matches scan_tissue_attempt/workhorse.R.

library(data.table)
library(matrixStats)


inverse_normal_transform <- function(x) {
  qnorm(
    (rank(x, na.last = "keep") - 0.5) /
      sum(!is.na(x))
  )
}


sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}


qc_filter_geno_for_expression_plot <- function(
    X,
    gene_name,
    hwe_threshold,
    maf_threshold) {

  X <- as.matrix(X)
  storage.mode(X) <- "integer"

  allele_frequency <- colMeans(X) / 2
  flip <- is.finite(allele_frequency) & allele_frequency > 0.5

  if (any(flip)) {
    X[, flip] <- 2L - X[, flip, drop = FALSE]
  }

  maf <- colMeans(X) / 2
  keep_maf <- is.finite(maf) & maf > maf_threshold
  X <- X[, keep_maf, drop = FALSE]

  if (ncol(X) == 0L) {
    stop("No SNPs remained after MAF filtering for ", gene_name, ".")
  }

  n_samples <- nrow(X)
  maf <- colMeans(X) / 2
  count0 <- colSums(X == 0L)
  count1 <- colSums(X == 1L)
  count2 <- colSums(X == 2L)

  expected0 <- n_samples * (1 - maf)^2
  expected1 <- n_samples * 2 * maf * (1 - maf)
  expected2 <- n_samples * maf^2

  valid_hwe <- expected0 > 0 & expected1 > 0 & expected2 > 0
  hwe_statistic <- rep(NA_real_, ncol(X))

  hwe_statistic[valid_hwe] <- (
    (count0[valid_hwe] - expected0[valid_hwe])^2 /
      expected0[valid_hwe] +
      (count1[valid_hwe] - expected1[valid_hwe])^2 /
      expected1[valid_hwe] +
      (count2[valid_hwe] - expected2[valid_hwe])^2 /
      expected2[valid_hwe]
  )

  hwe_p <- pchisq(
    hwe_statistic,
    df = 1,
    lower.tail = FALSE
  )

  fail_hwe <- !is.na(hwe_p) & hwe_p < hwe_threshold
  X <- X[, !fail_hwe, drop = FALSE]

  if (ncol(X) == 0L) {
    stop("No SNPs remained after HWE filtering for ", gene_name, ".")
  }

  X
}


recode_fit_mix_predictors <- function(X) {

  X <- as.matrix(X)
  storage.mode(X) <- "integer"

  snp_names <- colnames(X)
  additive <- X
  recessive <- (X == 2L) * 1L
  dominant <- (X >= 1L) * 1L

  colnames(additive) <- paste0(snp_names, "__additive")
  colnames(recessive) <- paste0(snp_names, "__recessive")
  colnames(dominant) <- paste0(snp_names, "__dominant")

  mixed <- cbind(additive, recessive, dominant)
  storage.mode(mixed) <- "double"
  rownames(mixed) <- rownames(X)
  mixed
}


get_fit_mix_cs_leads <- function(fit_mix, predictor_map) {

  cs_list <- fit_mix$sets[["cs"]]

  if (is.null(cs_list) || length(cs_list) == 0L) {
    return(data.table())
  }

  required_columns <- c(
    "predictor_index",
    "predictor_name",
    "snp",
    "coding"
  )

  missing_columns <- setdiff(required_columns, names(predictor_map))

  if (length(missing_columns) > 0L) {
    stop(
      "mix_predictor_map is missing: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  cs_index <- fit_mix$sets$cs_index
  cs_names <- names(cs_list)

  lead_rows <- lapply(
    seq_along(cs_list),
    function(cs_number) {

      members <- as.integer(cs_list[[cs_number]])
      members <- members[
        members >= 1L & members <= ncol(fit_mix$alpha)
      ]

      if (length(members) == 0L) {
        stop("fit_mix CS ", cs_number, " has no valid predictors.")
      }

      component_index <- NA_integer_

      if (
        !is.null(cs_index) &&
          length(cs_index) >= cs_number &&
          is.finite(cs_index[cs_number]) &&
          cs_index[cs_number] >= 1L &&
          cs_index[cs_number] <= nrow(fit_mix$alpha)
      ) {
        component_index <- as.integer(cs_index[cs_number])
      }

      if (is.finite(component_index)) {
        lead_score <- fit_mix$alpha[component_index, members]
      } else if (!is.null(fit_mix$pip)) {
        lead_score <- fit_mix$pip[members]
      } else {
        lead_score <- rep(NA_real_, length(members))
      }

      lead_score[!is.finite(lead_score)] <- -Inf

      lead_predictor_index <- if (all(lead_score == -Inf)) {
        members[1]
      } else {
        members[which.max(lead_score)]
      }

      map_row <- match(
        lead_predictor_index,
        predictor_map$predictor_index
      )

      if (is.na(map_row)) {
        stop(
          "Lead predictor ",
          lead_predictor_index,
          " is absent from mix_predictor_map."
        )
      }

      cs_name <- if (
        !is.null(cs_names) &&
          length(cs_names) >= cs_number &&
          nzchar(cs_names[cs_number])
      ) {
        cs_names[cs_number]
      } else {
        paste0("CS", cs_number)
      }

      data.table(
        cs_number = cs_number,
        cs_name = cs_name,
        component_index = component_index,
        lead_predictor_index = lead_predictor_index,
        lead_predictor_name = as.character(
          predictor_map$predictor_name[map_row]
        ),
        lead_snp = as.character(predictor_map$snp[map_row]),
        lead_coding = as.character(predictor_map$coding[map_row]),
        lead_pip = if (!is.null(fit_mix$pip)) {
          as.numeric(fit_mix$pip[lead_predictor_index])
        } else {
          NA_real_
        },
        cs_size = length(members)
      )
    }
  )

  rbindlist(lead_rows)
}


add_mix_coding_boundaries <- function(predictor_map) {

  required_columns <- c("predictor_index", "coding")
  missing_columns <- setdiff(required_columns, names(predictor_map))

  if (length(missing_columns) > 0L) {
    stop(
      "mix_predictor_map is missing: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  ordered_map <- predictor_map[
    order(predictor_map$predictor_index),
    ,
    drop = FALSE
  ]

  coding <- as.character(ordered_map$coding)

  if (length(coding) == 0L) {
    return(invisible(NULL))
  }

  coding_runs <- rle(coding)
  run_ends <- cumsum(coding_runs$lengths)
  run_starts <- c(1L, head(run_ends, -1L) + 1L)
  run_midpoints <- (run_starts + run_ends) / 2

  if (length(run_ends) > 1L) {
    abline(
      v = head(run_ends, -1L) + 0.5,
      col = "#B22222",
      lty = 2,
      lwd = 1.4
    )
  }

  plot_limits <- par("usr")

  text(
    x = run_midpoints,
    y = plot_limits[4],
    labels = tools::toTitleCase(coding_runs$values),
    pos = 3,
    offset = 0.20,
    xpd = NA,
    cex = 0.68,
    col = "#7A1F1F"
  )

  invisible(NULL)
}


finite_plot_limits <- function(...) {

  values <- unlist(list(...), use.names = FALSE)
  values <- values[is.finite(values)]

  if (length(values) == 0L) {
    stop("The expression values are all non-finite.")
  }

  limits <- range(values)
  spread <- diff(limits)
  padding <- if (spread > 0) 0.06 * spread else 0.25

  limits + c(-padding, padding)
}


open_four_panel_png <- function(filename) {

  png_arguments <- list(
    filename = filename,
    width = 3200,
    height = 2800,
    res = 250
  )

  if (capabilities("cairo")) {
    png_arguments$type <- "cairo-png"
  }

  do.call(png, png_arguments)
}


open_expression_png <- function(filename) {

  png_arguments <- list(
    filename = filename,
    width = 2200,
    height = 1800,
    res = 250
  )

  if (capabilities("cairo")) {
    png_arguments$type <- "cairo-png"
  }

  do.call(png, png_arguments)
}


plot_fit_mix_lead_expression <- function(
    gene_name,
    tissue_name,
    lead,
    raw_genotype,
    normalized_expression,
    output_file,
    genotype_axis_ticks = FALSE) {

  genotype_factor <- factor(raw_genotype, levels = 0:2)
  genotype_counts <- table(genotype_factor)
  expression_groups <- split(
    normalized_expression,
    genotype_factor,
    drop = FALSE
  )

  genotype_labels <- sprintf(
    "%d\n(n = %d)",
    0:2,
    as.integer(genotype_counts)
  )

  plot_title <- paste0(
    gene_name,
    " in ",
    tissue_name,
    " — fit_mix ",
    lead$cs_name,
    "\n",
    lead$lead_snp,
    " (",
    lead$lead_coding,
    " lead predictor; PIP = ",
    if (is.finite(lead$lead_pip)) {
      sprintf("%.3f", lead$lead_pip)
    } else {
      "NA"
    },
    ")"
  )

  open_expression_png(output_file)
  on.exit(invisible(dev.off()), add = TRUE)

  par(
    mar = c(5, 5, 5, 2) + 0.1,
    las = 1
  )

  boxplot(
    expression_groups,
    names = genotype_labels,
    xaxt = "n",
    col = c("#DCEAF7", "#91BCE2", "#4C78A8"),
    border = "#254E70",
    outline = FALSE,
    ylab = "Normalized gene expression",
    xlab = paste0(
      "Minor-allele dosage for ",
      lead$lead_snp
    ),
    main = plot_title,
    cex.main = 0.9
  )

  # Draw the genotype labels without vertical tick marks. The default
  # ticks run through the centered 0/1/2 labels on high-resolution plots.
  axis(
    side = 1,
    at = seq_along(genotype_labels),
    labels = genotype_labels,
    tick = genotype_axis_ticks,
    las = 1
  )

  set.seed(1L + lead$cs_number)
  stripchart(
    expression_groups,
    vertical = TRUE,
    method = "jitter",
    jitter = 0.18,
    add = TRUE,
    pch = 21,
    cex = 0.75,
    col = adjustcolor("#1F2937", alpha.f = 0.55),
    bg = adjustcolor("white", alpha.f = 0.50)
  )

  group_means <- vapply(
    expression_groups,
    function(x) {
      if (length(x) > 0L) mean(x) else NA_real_
    },
    numeric(1L)
  )

  valid_means <- is.finite(group_means)

  points(
    which(valid_means),
    group_means[valid_means],
    pch = 23,
    cex = 1.3,
    lwd = 1.2,
    col = "#7F2704",
    bg = "#F28E2B"
  )

  legend(
    "topright",
    legend = "Group mean",
    pch = 23,
    pt.bg = "#F28E2B",
    col = "#7F2704",
    bty = "n",
    cex = 0.8
  )

  invisible(NULL)
}


conditional_expression_for_cs <- function(
    y,
    geno_mix,
    lead_indices,
    focal_cs) {

  if (
    focal_cs < 1L ||
      focal_cs > length(lead_indices)
  ) {
    stop("focal_cs is outside the lead-index vector.")
  }

  focal_index <- lead_indices[focal_cs]

  other_indices <- setdiff(
    unique(lead_indices[-focal_cs]),
    focal_index
  )

  if (length(other_indices) == 0L) {
    return(y)
  }

  focal_predictor <- geno_mix[
    ,
    focal_index,
    drop = FALSE
  ]

  other_predictors <- geno_mix[
    ,
    other_indices,
    drop = FALSE
  ]

  # Estimate the other-CS effects jointly with the focal predictor. This
  # prevents LD between the leads from making the other predictor absorb
  # part of the focal effect.
  conditional_fit <- lm.fit(
    x = cbind(
      intercept = 1,
      focal_predictor,
      other_predictors
    ),
    y = y
  )

  other_coefficient_positions <- seq.int(
    from = 3L,
    length.out = ncol(other_predictors)
  )
  other_coefficients <- conditional_fit$coefficients[
    other_coefficient_positions
  ]

  if (any(!is.finite(other_coefficients))) {
    stop(
      "The focal and other CS lead predictors are collinear; ",
      "their conditional effects cannot be separated."
    )
  }

  # Remove only the slopes for the other CS leads. Keeping the intercept
  # leaves the adjusted phenotype on the same location as the original y.
  y - as.numeric(
    other_predictors %*% other_coefficients
  )
}


plot_fit_mix_expression_panel <- function(
    lead,
    raw_genotype,
    expression_to_plot,
    main,
    ylab,
    ylim = NULL,
    show_mean_legend = TRUE,
    jitter_seed = 1L,
    genotype_axis_ticks = FALSE) {

  genotype_factor <- factor(raw_genotype, levels = 0:2)
  genotype_counts <- table(genotype_factor)
  expression_groups <- split(
    expression_to_plot,
    genotype_factor,
    drop = FALSE
  )

  genotype_labels <- sprintf(
    "%d\n(n = %d)",
    0:2,
    as.integer(genotype_counts)
  )

  boxplot(
    expression_groups,
    names = genotype_labels,
    xaxt = "n",
    col = c("#DCEAF7", "#91BCE2", "#4C78A8"),
    border = "#254E70",
    outline = FALSE,
    ylim = ylim,
    ylab = ylab,
    xlab = paste0(
      "Minor-allele dosage for ",
      lead$lead_snp
    ),
    main = main,
    cex.main = 0.82
  )

  axis(
    side = 1,
    at = seq_along(genotype_labels),
    labels = genotype_labels,
    tick = genotype_axis_ticks,
    las = 1
  )

  set.seed(jitter_seed)
  stripchart(
    expression_groups,
    vertical = TRUE,
    method = "jitter",
    jitter = 0.18,
    add = TRUE,
    pch = 21,
    cex = 0.65,
    col = adjustcolor("#1F2937", alpha.f = 0.55),
    bg = adjustcolor("white", alpha.f = 0.50)
  )

  group_means <- vapply(
    expression_groups,
    function(x) {
      if (length(x) > 0L) mean(x) else NA_real_
    },
    numeric(1L)
  )

  valid_means <- is.finite(group_means)

  points(
    which(valid_means),
    group_means[valid_means],
    pch = 23,
    cex = 1.2,
    lwd = 1.1,
    col = "#7F2704",
    bg = "#F28E2B"
  )

  abline(
    h = 0,
    col = "gray75",
    lty = 3
  )

  if (show_mean_legend) {
    legend(
      "topright",
      legend = "Group mean",
      pch = 23,
      pt.bg = "#F28E2B",
      col = "#7F2704",
      bty = "n",
      cex = 0.72
    )
  }

  invisible(genotype_counts)
}


load_fit_mix_expression_inputs <- function(
    subject_pheno_file,
    sample_attr_file,
    expr_file,
    gtf_file,
    gene_annotation_function) {

  source(gene_annotation_function)

  cat("Importing covariate data.\n")

  cov_subject <- read.table(
    subject_pheno_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  cov_sample <- read.table(
    sample_attr_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )

  cov_sample <- transform(
    cov_sample,
    SUBJID = substr(SAMPID, 1, 10)
  )

  cov <- merge(cov_subject, cov_sample, by = "SUBJID")

  cov <- cov[
    c(
      "SUBJID",
      "SAMPID",
      "SEX",
      "AGE",
      "SMTS",
      "SMTSD",
      "SMGEBTCHT",
      "SMAFRZE"
    )
  ]

  cov <- subset(cov, SMAFRZE == "RNASEQ")

  cov <- transform(
    cov,
    SEX = SEX - 1,
    AGE = factor(AGE),
    SMTS = factor(SMTS),
    SMTSD = factor(SMTSD),
    SMGEBTCHT = factor(SMGEBTCHT),
    SMAFRZE = factor(SMAFRZE)
  )

  rownames(cov) <- cov$SAMPID

  cat("Importing gene-expression data.\n")

  pheno_all <- fread(
    expr_file,
    sep = "\t",
    skip = 2,
    header = TRUE,
    showProgress = TRUE
  )

  class(pheno_all) <- "data.frame"
  gene_info <- pheno_all[1:2]
  pheno_all <- pheno_all[-(1:2)]
  pheno_all <- as.matrix(pheno_all)
  pheno_all <- t(pheno_all)
  storage.mode(pheno_all) <- "double"
  colnames(pheno_all) <- gene_info$Name

  sample_ids <- intersect(cov$SAMPID, rownames(pheno_all))

  cov <- cov[
    match(sample_ids, cov$SAMPID),
    ,
    drop = FALSE
  ]

  pheno_all <- pheno_all[
    match(sample_ids, rownames(pheno_all)),
    ,
    drop = FALSE
  ]

  list(
    cov = cov,
    pheno_all = pheno_all,
    gene_info = gene_info,
    total_read_count = rowSums(pheno_all, na.rm = TRUE),
    gene_annotations = get_gene_annotations(gtf_file)
  )
}


prepare_fit_mix_gene_data <- function(
    gene_name,
    shared_inputs,
    results_dir,
    temp_dir,
    plink_exec,
    geno_file,
    cis_window,
    min_maf_plink,
    min_maf,
    hwe_thresh) {

  result_file <- file.path(results_dir, paste0(gene_name, ".rds"))

  if (!file.exists(result_file)) {
    stop("Result file does not exist: ", result_file)
  }

  gene_results <- readRDS(result_file)

  expression_column <- which(
    shared_inputs$gene_info$Description == gene_name
  )

  if (length(expression_column) == 0L) {
    stop("Gene not found in expression data: ", gene_name)
  }

  if (length(expression_column) > 1L) {
    warning(
      "More than one expression column found for ",
      gene_name,
      "; using the first."
    )
    expression_column <- expression_column[1]
  }

  pheno_gene <- cbind(
    shared_inputs$cov,
    data.frame(
      count = shared_inputs$pheno_all[, expression_column],
      total_read_count = shared_inputs$total_read_count
    )
  )

  gene_annotation <- shared_inputs$gene_annotations[
    shared_inputs$gene_annotations$gene_name == gene_name,
    ,
    drop = FALSE
  ]

  if (nrow(gene_annotation) == 0L) {
    stop("Gene annotation not found: ", gene_name)
  }

  chromosome <- sub(
    "^chr",
    "",
    as.character(gene_annotation$chromosome[1])
  )

  if (chromosome == "M") {
    chromosome <- "MT"
  }

  tss <- with(
    gene_annotation[1, ],
    ifelse(strand == "+", start, end)
  )

  region_start <- max(0, tss - cis_window)
  region_end <- tss + cis_window

  plink_prefix <- tempfile(
    pattern = paste0(
      "plot_",
      sanitize_filename(gene_name),
      "_"
    ),
    tmpdir = temp_dir
  )

  on.exit(
    unlink(
      Sys.glob(paste0(plink_prefix, ".*")),
      force = TRUE
    ),
    add = TRUE
  )

  cat("Extracting cis-genotypes for ", gene_name, ".\n", sep = "")

  plink_status <- system2(
    plink_exec,
    args = c(
      "--bfile",
      shQuote(geno_file),
      "--chr",
      chromosome,
      "--from-bp",
      format(region_start, scientific = FALSE),
      "--to-bp",
      format(region_end, scientific = FALSE),
      "--snps-only",
      "--max-alleles",
      "2",
      "--rm-dup",
      "exclude-all",
      "--threads",
      "2",
      "--memory",
      "8000",
      "--maf",
      format(min_maf_plink, scientific = FALSE),
      "--recode",
      "A",
      "--out",
      shQuote(plink_prefix)
    )
  )

  raw_file <- paste0(plink_prefix, ".raw")

  if (plink_status != 0L || !file.exists(raw_file)) {
    stop("PLINK failed with exit status ", plink_status, ".")
  }

  geno_all <- fread(raw_file, sep = "\t", header = TRUE)
  class(geno_all) <- "data.frame"
  rownames(geno_all) <- geno_all$IID
  geno_all <- as.matrix(geno_all[, -(1:6), drop = FALSE])
  storage.mode(geno_all) <- "double"

  keep_complete <- colSums(is.na(geno_all)) == 0L
  geno_all <- geno_all[, keep_complete, drop = FALSE]

  keep_variable <- matrixStats::colSds(geno_all) > 0
  geno_all <- geno_all[, keep_variable, drop = FALSE]

  geno_all <- qc_filter_geno_for_expression_plot(
    X = geno_all,
    gene_name = gene_name,
    hwe_threshold = hwe_thresh,
    maf_threshold = min_maf
  )

  list(
    gene_results = gene_results,
    pheno_gene = pheno_gene,
    geno_all = geno_all,
    geno_mix_all = recode_fit_mix_predictors(geno_all)
  )
}


prepare_fit_mix_tissue_data <- function(
    gene_data,
    gene_name,
    tissue_name,
    tissue_result,
    min_samples,
    min_n_rec) {

  pheno <- subset(
    gene_data$pheno_gene,
    SMTS == tissue_name
  )

  matched_ids <- intersect(
    pheno$SUBJID,
    rownames(gene_data$geno_all)
  )

  pheno <- pheno[
    match(matched_ids, pheno$SUBJID),
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(
      as.character(pheno$SUBJID),
      as.character(matched_ids)
    )
  )

  keep_sample <- (
    is.finite(pheno$count) &
      is.finite(pheno$total_read_count) &
      pheno$total_read_count > 0 &
      !is.na(pheno$SEX)
  )

  pheno <- pheno[keep_sample, , drop = FALSE]
  matched_ids <- as.character(pheno$SUBJID)

  if (nrow(pheno) < min_samples) {
    stop(
      "Only ",
      nrow(pheno),
      " samples remained for ",
      gene_name,
      " / ",
      tissue_name,
      "."
    )
  }

  if (
    !is.null(tissue_result$n_ind) &&
      is.finite(tissue_result$n_ind) &&
      nrow(pheno) != tissue_result$n_ind
  ) {
    stop(
      "Reconstructed sample count (",
      nrow(pheno),
      ") does not match the saved fit (",
      tissue_result$n_ind,
      ")."
    )
  }

  pheno$library_size_factor <- (
    pheno$total_read_count /
      mean(pheno$total_read_count)
  )

  pheno$normalized_expression <- log1p(
    pheno$count /
      pheno$library_size_factor
  )

  pheno$normalized_expression <- inverse_normal_transform(
    pheno$normalized_expression
  )

  pheno$y <- resid(
    lm(
      normalized_expression ~ SEX,
      data = pheno
    )
  )

  pheno$y <- inverse_normal_transform(pheno$y)

  geno_for_counts <- gene_data$geno_all[
    matched_ids,
    ,
    drop = FALSE
  ]

  geno_mix <- gene_data$geno_mix_all[
    matched_ids,
    ,
    drop = FALSE
  ]

  keep_mix <- (
    colSums(geno_mix) >= min_n_rec &
      matrixStats::colSds(geno_mix) > 0
  )

  geno_mix <- geno_mix[, keep_mix, drop = FALSE]

  if (is.null(tissue_result$mix_predictor_map)) {
    stop(
      "mix_predictor_map is absent; rerun ",
      gene_name,
      " with the current workhorse."
    )
  }

  predictor_map <- as.data.frame(tissue_result$mix_predictor_map)

  expected_names <- as.character(
    predictor_map$predictor_name[
      order(predictor_map$predictor_index)
    ]
  )

  fit_mix <- tissue_result$susie_mix

  if (length(expected_names) != ncol(fit_mix$alpha)) {
    stop("mix_predictor_map length does not match fit_mix$alpha.")
  }

  predictor_order <- match(expected_names, colnames(geno_mix))

  if (anyNA(predictor_order)) {
    stop(
      "Could not reconstruct saved fit_mix predictors: ",
      paste(expected_names[is.na(predictor_order)], collapse = ", ")
    )
  }

  geno_mix <- geno_mix[, predictor_order, drop = FALSE]

  if (!identical(colnames(geno_mix), expected_names)) {
    stop("Reconstructed fit_mix predictor order is incorrect.")
  }

  list(
    y = pheno$y,
    geno_for_counts = geno_for_counts,
    geno_mix = geno_mix,
    predictor_map = predictor_map
  )
}


get_fit_add_plot_lead <- function(fit_add, predictor_map) {

  required_columns <- c("predictor_index", "predictor_name", "snp", "coding")
  if (!is.data.frame(predictor_map) ||
      !all(required_columns %in% names(predictor_map))) {
    stop("add_predictor_map is absent or missing required columns.")
  }
  if (nrow(predictor_map) != ncol(fit_add$alpha) ||
      !identical(sort(as.integer(predictor_map$predictor_index)),
                 seq_len(ncol(fit_add$alpha))) ||
      anyNA(predictor_map$snp) || any(!nzchar(predictor_map$snp)) ||
      anyNA(predictor_map$coding) || any(predictor_map$coding != "additive")) {
    stop("add_predictor_map does not match the saved additive fit.")
  }
  if (length(fit_add$pip) != ncol(fit_add$alpha) ||
      !any(is.finite(fit_add$pip))) {
    stop("The saved additive fit has no usable PIPs.")
  }

  # A single CS uses its component-alpha lead, matching the mixed model.
  # Keep the original case selection, which only restricts the mixed CS count:
  # for multiple additive CSs, plot the CS lead with the highest PIP.
  leads <- get_fit_mix_cs_leads(fit_add, predictor_map)
  if (nrow(leads) > 0L) {
    scores <- leads$lead_pip
    scores[!is.finite(scores)] <- -Inf
    lead <- leads[which.max(scores)]
    lead[, lead_selection := if (nrow(leads) == 1L) {
      "Single credible-set lead"
    } else {
      "Highest-PIP credible-set lead"
    }]
    return(lead)
  }

  # With no additive CS, explicitly label the highest-PIP SNP as a fallback.
  scores <- fit_add$pip
  scores[!is.finite(scores)] <- -Inf
  lead_index <- which.max(scores)
  map_row <- match(lead_index, predictor_map$predictor_index)
  data.table(
    cs_number = NA_integer_,
    cs_name = NA_character_,
    component_index = NA_integer_,
    lead_predictor_index = lead_index,
    lead_predictor_name = as.character(predictor_map$predictor_name[map_row]),
    lead_snp = as.character(predictor_map$snp[map_row]),
    lead_coding = "additive",
    lead_pip = as.numeric(fit_add$pip[lead_index]),
    cs_size = NA_integer_,
    lead_selection = "Highest-PIP SNP (no credible set)"
  )
}


get_one_cs_likelihood_comparison <- function(fit_add, fit_mix) {

  # Match get_log_lik_metric() in generate_summary_results.R exactly.
  # This is the existing ELBO + KL log-likelihood-like metric, not a pair
  # of maximized nested-model likelihoods with a calibrated chi-square test.
  log_lik_metric <- function(fit) {
    elbo <- fit$elbo
    elbo <- elbo[is.finite(elbo)]
    if (length(elbo) == 0L || is.null(fit$KL)) {
      return(NA_real_)
    }
    tail(elbo, 1) + sum(fit$KL, na.rm = TRUE)
  }

  log_lik_add <- log_lik_metric(fit_add)
  log_lik_mix <- log_lik_metric(fit_mix)
  list(
    log_lik_add = log_lik_add,
    log_lik_mix = log_lik_mix,
    likelihood_statistic = -2 * (log_lik_add - log_lik_mix)
  )
}


describe_one_cs_lead_change <- function(add_lead, mix_lead) {
  paste0(
    if (identical(add_lead$lead_snp, mix_lead$lead_snp)) {
      "Same lead SNP; coding change only: "
    } else {
      "Different lead SNP; coding change: "
    },
    "additive -> ", mix_lead$lead_coding
  )
}


plot_one_cs_four_panel_case <- function(
    gene_name,
    tissue_name,
    fit_add,
    fit_mix,
    predictor_map,
    add_lead,
    mix_lead,
    raw_genotype_add,
    raw_genotype_mix,
    normalized_expression,
    likelihood_comparison,
    output_file,
    genotype_axis_ticks = FALSE) {

  expression_limits <- finite_plot_limits(normalized_expression)
  n_add_cs <- length(fit_add$sets[["cs"]])
  lead_comparison <- describe_one_cs_lead_change(add_lead, mix_lead)
  statistic <- likelihood_comparison$likelihood_statistic
  statistic_label <- if (is.finite(statistic)) {
    sprintf(
      "Likelihood-ratio statistic: 2(log lik mix - log lik additive) = %.2f",
      statistic
    )
  } else {
    "Likelihood-ratio statistic unavailable (missing finite ELBO/KL metric)"
  }

  open_four_panel_png(output_file)
  on.exit(invisible(dev.off()), add = TRUE)
  par(
    mfrow = c(2, 2),
    oma = c(1, 1, 6.5, 1),
    mar = c(5, 5, 4.5, 2) + 0.1,
    las = 1
  )

  susieR::susie_plot(
    fit_add,
    y = "PIP",
    main = paste0("Additive SuSiE: ", n_add_cs, " credible set",
                  if (n_add_cs == 1L) "" else "s")
  )
  susieR::susie_plot(
    fit_mix,
    y = "PIP",
    main = "SuSiE-mix: 1 credible set"
  )
  add_mix_coding_boundaries(predictor_map)

  add_title <- paste0(
    "SuSiE lead: ", add_lead$lead_snp,
    "\nAdditive; PIP = ", sprintf("%.3f", add_lead$lead_pip),
    if (n_add_cs == 0L) {
      " (highest PIP; no CS)"
    } else if (n_add_cs > 1L) {
      paste0(" (highest-PIP CS lead; ", add_lead$cs_name, ")")
    } else {
      paste0(" (", add_lead$cs_name, ")")
    }
  )
  mix_title <- paste0(
    "SuSiE-mix lead: ", mix_lead$lead_snp,
    "\n", tools::toTitleCase(mix_lead$lead_coding),
    "; PIP = ", sprintf("%.3f", mix_lead$lead_pip),
    " (", mix_lead$cs_name, ")"
  )

  # Both panels show the same workhorse-normalized phenotype, grouped by
  # each biological SNP's raw 0/1/2 dosage. Do not condition one model on the
  # other: these are alternative explanations of a single mixed-model CS.
  # Use identical jitter so a coding-only change gives identical boxplots.
  plot_fit_mix_expression_panel(
    lead = add_lead,
    raw_genotype = raw_genotype_add,
    expression_to_plot = normalized_expression,
    main = add_title,
    ylab = "Normalized gene expression",
    ylim = expression_limits,
    show_mean_legend = FALSE,
    jitter_seed = 101L,
    genotype_axis_ticks = genotype_axis_ticks
  )
  plot_fit_mix_expression_panel(
    lead = mix_lead,
    raw_genotype = raw_genotype_mix,
    expression_to_plot = normalized_expression,
    main = mix_title,
    ylab = "Normalized gene expression",
    ylim = expression_limits,
    show_mean_legend = TRUE,
    jitter_seed = 101L,
    genotype_axis_ticks = genotype_axis_ticks
  )

  mtext(paste(gene_name, tissue_name, sep = " - "), side = 3,
        outer = TRUE, line = 4.7, font = 2, cex = 1.25)
  mtext(lead_comparison, side = 3, outer = TRUE, line = 3.2, cex = 0.95)
  mtext(statistic_label, side = 3, outer = TRUE, line = 1.8, cex = 0.85)
  mtext("ELBO + KL summary metric; positive favors SuSiE-mix; no calibrated p-value",
        side = 3, outer = TRUE, line = 0.6, cex = 0.70, col = "gray35")

  invisible(NULL)
}


run_one_cs_fit_mix_plots <- function(
    cases,
    expected_coding,
    plot_dir,
    summary_filename,
    project_dir = "/project2/mstephens/wdenault/susie_mix",
    datadir = "/project2/mstephens/gtex",
    min_maf_plink = 0.00,
    min_maf = 0.05,
    hwe_thresh = 1e-8,
    min_n_rec = 5,
    cis_window = 5e5,
    min_samples = 50,
    genotype_axis_ticks = FALSE) {

  expected_coding <- match.arg(
    expected_coding,
    c("dominant", "recessive")
  )

  cases <- unique(
    as.data.table(cases)[, .(gene, tissue)],
    by = c("gene", "tissue")
  )

  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(cases) == 0L) {
    message("No one-CS ", expected_coding, " cases met the filters.")
    empty_summary <- data.table(
      gene = character(),
      tissue = character(),
      status = character(),
      message = character()
    )
    fwrite(empty_summary, file.path(plot_dir, summary_filename))
    return(invisible(empty_summary))
  }

  results_dir <- file.path(project_dir, "results")
  temp_dir <- file.path(project_dir, "temp_plink")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  plink_exec <- file.path(datadir, "plink2")
  geno_file <- file.path(
    datadir,
    "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"
  )
  subject_pheno_file <- file.path(
    datadir,
    "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz"
  )
  sample_attr_file <- file.path(
    datadir,
    "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz"
  )
  expr_file <- file.path(
    datadir,
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
  )
  gtf_file <- file.path(
    datadir,
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"
  )
  gene_annotation_function <- file.path(
    project_dir,
    "script",
    "scan_tissue_attempt",
    "get_gene_annotations.R"
  )

  shared_inputs <- load_fit_mix_expression_inputs(
    subject_pheno_file = subject_pheno_file,
    sample_attr_file = sample_attr_file,
    expr_file = expr_file,
    gtf_file = gtf_file,
    gene_annotation_function = gene_annotation_function
  )

  result_rows <- list()
  result_index <- 0L

  add_result <- function(result) {
    result_index <<- result_index + 1L
    result_rows[[result_index]] <<- result
  }

  make_error_row <- function(gene_name, tissue_name, message) {
    data.table(
      gene = gene_name,
      tissue = tissue_name,
      expected_coding = expected_coding,
      status = "error",
      message = message
    )
  }

  for (gene_name in unique(cases$gene)) {

    gene_cases <- cases[gene == gene_name]

    gene_data <- tryCatch(
      prepare_fit_mix_gene_data(
        gene_name = gene_name,
        shared_inputs = shared_inputs,
        results_dir = results_dir,
        temp_dir = temp_dir,
        plink_exec = plink_exec,
        geno_file = geno_file,
        cis_window = cis_window,
        min_maf_plink = min_maf_plink,
        min_maf = min_maf,
        hwe_thresh = hwe_thresh
      ),
      error = function(e) e
    )

    if (inherits(gene_data, "error")) {
      error_message <- conditionMessage(gene_data)

      for (tissue_name in gene_cases$tissue) {
        message(
          "ERROR [",
          gene_name,
          " / ",
          tissue_name,
          "]: ",
          error_message
        )
        add_result(make_error_row(gene_name, tissue_name, error_message))
      }

      next
    }

    for (tissue_name in gene_cases$tissue) {

      case_result <- tryCatch(
        {
          if (!tissue_name %in% names(gene_data$gene_results)) {
            stop("Tissue is absent from the saved gene result.")
          }

          tissue_result <- gene_data$gene_results[[tissue_name]]
          fit_add <- tissue_result$susie_add
          fit_mix <- tissue_result$susie_mix

          if (is.null(fit_add)) {
            stop("The saved tissue result has no susie_add fit.")
          }

          if (is.null(fit_mix)) {
            stop("The saved tissue result has no susie_mix fit.")
          }

          cs_list <- fit_mix$sets[["cs"]]
          n_cs <- if (is.null(cs_list)) 0L else length(cs_list)

          if (n_cs != 1L) {
            stop(
              "Expected exactly one fit_mix CS, but the saved fit has ",
              n_cs,
              "."
            )
          }

          tissue_data <- prepare_fit_mix_tissue_data(
            gene_data = gene_data,
            gene_name = gene_name,
            tissue_name = tissue_name,
            tissue_result = tissue_result,
            min_samples = min_samples,
            min_n_rec = min_n_rec
          )

          lead <- get_fit_mix_cs_leads(
            fit_mix,
            tissue_data$predictor_map
          )[1]

          add_lead <- get_fit_add_plot_lead(
            fit_add,
            tissue_result$add_predictor_map
          )
          likelihood_comparison <- get_one_cs_likelihood_comparison(
            fit_add,
            fit_mix
          )

          if (!identical(lead$lead_coding, expected_coding)) {
            stop(
              "Summary selected this as ",
              expected_coding,
              ", but the saved CS lead coding is ",
              lead$lead_coding,
              "."
            )
          }

          missing_lead_snps <- setdiff(
            c(add_lead$lead_snp, lead$lead_snp),
            colnames(tissue_data$geno_for_counts)
          )
          if (length(missing_lead_snps) > 0L) {
            stop(
              "Lead SNPs are absent from reconstructed genotypes: ",
              paste(missing_lead_snps, collapse = ", ")
            )
          }

          raw_genotype <- tissue_data$geno_for_counts[, lead$lead_snp]
          raw_genotype_add <- tissue_data$geno_for_counts[, add_lead$lead_snp]
          output_file <- file.path(
            plot_dir,
            paste0(
              sanitize_filename(gene_name),
              "_",
              sanitize_filename(tissue_name),
              "_add_vs_mix.png"
            )
          )

          plot_one_cs_four_panel_case(
            gene_name = gene_name,
            tissue_name = tissue_name,
            fit_add = fit_add,
            fit_mix = fit_mix,
            predictor_map = tissue_data$predictor_map,
            add_lead = add_lead,
            mix_lead = lead,
            raw_genotype_add = raw_genotype_add,
            raw_genotype_mix = raw_genotype,
            normalized_expression = tissue_data$y,
            likelihood_comparison = likelihood_comparison,
            output_file = output_file,
            genotype_axis_ticks = genotype_axis_ticks
          )

          data.table(
            gene = gene_name,
            tissue = tissue_name,
            expected_coding = expected_coding,
            status = "plotted",
            message = "Four-panel SuSiE / SuSiE-mix comparison created",
            additive_n_cs = length(fit_add$sets[["cs"]]),
            additive_cs_name = add_lead$cs_name,
            additive_component_index = add_lead$component_index,
            additive_lead_predictor_index = add_lead$lead_predictor_index,
            additive_lead_predictor = add_lead$lead_predictor_name,
            additive_lead_snp = add_lead$lead_snp,
            additive_lead_coding = add_lead$lead_coding,
            additive_lead_pip = add_lead$lead_pip,
            additive_cs_size = add_lead$cs_size,
            additive_lead_selection = add_lead$lead_selection,
            additive_n_genotype_0 = sum(raw_genotype_add == 0L),
            additive_n_genotype_1 = sum(raw_genotype_add == 1L),
            additive_n_genotype_2 = sum(raw_genotype_add == 2L),
            same_lead_snp = identical(add_lead$lead_snp, lead$lead_snp),
            lead_comparison = describe_one_cs_lead_change(add_lead, lead),
            log_lik_add = likelihood_comparison$log_lik_add,
            log_lik_mix = likelihood_comparison$log_lik_mix,
            likelihood_statistic = likelihood_comparison$likelihood_statistic,
            n_cs = 1L,
            cs_name = lead$cs_name,
            component_index = lead$component_index,
            lead_predictor_index = lead$lead_predictor_index,
            lead_predictor = lead$lead_predictor_name,
            lead_snp = lead$lead_snp,
            lead_coding = lead$lead_coding,
            lead_pip = lead$lead_pip,
            cs_size = lead$cs_size,
            n_genotype_0 = sum(raw_genotype == 0L),
            n_genotype_1 = sum(raw_genotype == 1L),
            n_genotype_2 = sum(raw_genotype == 2L),
            output_file = output_file
          )
        },
        error = function(e) {
          error_message <- conditionMessage(e)
          message(
            "ERROR [",
            gene_name,
            " / ",
            tissue_name,
            "]: ",
            error_message
          )
          make_error_row(gene_name, tissue_name, error_message)
        }
      )

      add_result(case_result)
    }
  }

  plot_summary <- rbindlist(result_rows, fill = TRUE)
  summary_file <- file.path(plot_dir, summary_filename)
  fwrite(plot_summary, summary_file)

  cat(
    "\nFinished ",
    expected_coding,
    " one-CS plots. Plotted ",
    sum(plot_summary$status == "plotted"),
    " of ",
    nrow(plot_summary),
    " candidate gene-tissue pairs.\n",
    "Summary: ",
    summary_file,
    "\n",
    sep = ""
  )

  invisible(plot_summary)
}
