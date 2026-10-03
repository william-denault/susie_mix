# Fit only discrete-prior slider SuSiE with frozen tissue-specific weights.
# Data preparation matches workhorse_em.R; no CS/P/read/N>=300 EM screen.
run_slide_prior_gene <- function(
    target_gene = "GTF2H2",
    tissue_priors,

    # --- paths ---
    project_dir = Sys.getenv(
      "SUSIE_MIX_PROJECT_DIR", "/project2/mstephens/wdenault/susie_mix"
    ),
    datadir = "/project2/mstephens/gtex",
    plink_exec = file.path(datadir, "plink2"),
    gene_annot_fun = file.path(
      project_dir, "script/scan_tissue_attempt/get_gene_annotations.R"
    ),
    gtf_file = file.path(
      datadir,
      "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"
    ),
    geno_file = file.path(
      datadir,
      "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"
    ),
    subject_pheno_file = file.path(
      datadir,
      "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz"
    ),
    sample_attr_file = file.path(
      datadir,
      "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz"
    ),
    expr_file = file.path(
      datadir,
      "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
    ),

    # --- analysis parameters ---
    min_maf_plink = 0.00,
    min_maf = 0.05,
    cis_window = 5e5,
    min_samples = 50,
    seed = 1,
    hwe_thresh = 1e-8,
    min_n_rec = 5,
    min_obs = 5,
    delta_grid = seq(-1, 1, length.out = 17),

    # --- SuSiE parameters ---
    L = 10,
    standardize = FALSE,
    estimate_prior_method = "EM",
    min_abs_corr = 0.5,
    verbose = FALSE,
    max_iter = 1000,
    tol = 1e-5,

    # --- misc ---
    temp_dir = file.path(project_dir, "temp_plink_slide_prior_em"),
    previous_result = NULL
) {


  library(data.table)
  library(matrixStats)
  library(susieRSlidePrior)

  source(file.path(project_dir, "script/scan_tissue_attempt/workhorse_utils.R"),
         local = TRUE)
  source(file.path(project_dir, "script/scan_tissue_attempt/em_utils.R"),
         local = TRUE)
  source(file.path(project_dir, "script/scan_tissue_attempt/slide_prior_em_utils.R"),
         local = TRUE)
  spe_validate_priors(tissue_priors)
  source(gene_annot_fun, local = TRUE)

  set.seed(seed)

  dir.create(
    temp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  plink_out_prefix <- tempfile(paste0("plink_", target_gene, "_"), tmpdir = temp_dir)

  # Always remove gene-specific PLINK intermediates, including when an
  # error occurs before the end of the gene analysis.
  on.exit(
    unlink(
      Sys.glob(
        paste0(plink_out_prefix, ".*")
      ),
      force = TRUE
    ),
    add = TRUE
  )

  # ------------------------------------------------------------
  # Import covariates
  # ------------------------------------------------------------

  cat("Importing covariate data.\n")

  cov1 <- read.table(
    subject_pheno_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  cov2 <- read.table(
    sample_attr_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )

  cov2 <- transform(
    cov2,
    SUBJID = substr(SAMPID, 1, 10)
  )

  cov <- merge(
    cov1,
    cov2,
    by = "SUBJID"
  )

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

  cov <- subset(
    cov,
    SMAFRZE == "RNASEQ"
  )

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

  # ------------------------------------------------------------
  # Import gene-expression data
  # ------------------------------------------------------------

  cat("Importing gene expression data.\n")

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

  # Align expression and covariate data.
  ids <- intersect(
    cov$SAMPID,
    rownames(pheno_all)
  )

  rows1 <- match(
    ids,
    cov$SAMPID
  )

  rows2 <- match(
    ids,
    rownames(pheno_all)
  )

  cov <- cov[rows1, , drop = FALSE]
  pheno_all <- pheno_all[rows2, , drop = FALSE]

  # Extract expression for the target gene.
  j <- which(
    gene_info$Description == target_gene
  )

  if (length(j) == 0L) {
    stop("Target gene was not found: ", target_gene)
  }

  if (length(j) > 1L) {
    warning(
      "More than one expression column found for ",
      target_gene,
      "; using the first."
    )

    j <- j[1]
  }

  pheno_gene <- cbind(
    cov,
    data.frame(count = pheno_all[, j])
  )

  read_count <- rowSums(
    pheno_all,
    na.rm = TRUE
  )

  pheno_gene <- cbind(
    cov,
    data.frame(
      count = pheno_all[, j],
      total_read_count = read_count
    )
  )
  all_tissues <- levels(
    droplevels(pheno_gene$SMTS)
  )

  cat(
    "Found",
    length(all_tissues),
    "tissues for gene",
    target_gene,
    ":\n"
  )

  print(all_tissues)

  # ------------------------------------------------------------
  # Identify the cis-region
  # ------------------------------------------------------------

  genes <- get_gene_annotations(gtf_file)

  genes <- subset(
    genes,
    gene_name == target_gene
  )

  if (nrow(genes) == 0L) {
    stop(
      "No gene annotation found for: ",
      target_gene
    )
  }

  chr <- sub(
    "^chr",
    "",
    as.character(genes$chromosome[1])
  )

  if (chr == "M") {
    chr <- "MT"
  }

  if (
    length(chr) != 1L ||
    is.na(chr) ||
    !nzchar(chr)
  ) {
    stop(
      "Could not determine the chromosome for: ",
      target_gene
    )
  }

  tss <- with(
    genes[1, ],
    ifelse(strand == "+", start, end)
  )

  pos0 <- max(0, tss - cis_window)
  pos1 <- tss + cis_window

  # ------------------------------------------------------------
  # Extract cis-genotypes
  # ------------------------------------------------------------

  geno_file_raw <- paste0(
    plink_out_prefix,
    ".raw"
  )

  if (!file.exists(geno_file_raw)) {

    cat("Extracting genotype data from PLINK file.\n")

    plink_call <- sprintf(
      paste(
        "%s --bfile %s --chr %s --from-bp %d --to-bp %d",
        "--snps-only --max-alleles 2 --rm-dup exclude-all",
        "--threads 2 --memory 8000 --maf %g",
        "--recode A --out %s"
      ),
      plink_exec,
      geno_file,
      chr,
      pos0,
      pos1,
      min_maf_plink,
      plink_out_prefix
    )

    plink_status <- system(plink_call)

    if (
      plink_status != 0L ||
      !file.exists(geno_file_raw)
    ) {
      stop(
        "PLINK failed for ",
        target_gene,
        " on chromosome ",
        chr,
        " with exit status ",
        plink_status
      )
    }
  }

  geno_all <- fread(
    geno_file_raw,
    sep = "\t",
    header = TRUE
  )

  class(geno_all) <- "data.frame"

  ids <- geno_all$IID
  rownames(geno_all) <- ids

  geno_all <- geno_all[, -(1:6)]
  geno_all <- as.matrix(geno_all)
  dim(geno_all)
  storage.mode(geno_all) <- "double"

  # Remove SNPs with missing genotypes.
  keep <- colSums(is.na(geno_all)) == 0

  geno_all <- geno_all[
    ,
    keep,
    drop = FALSE
  ]

  # Remove SNPs that do not vary.
  keep <- colSds(geno_all) > 0

  geno_all <- geno_all[
    ,
    keep,
    drop = FALSE
  ]

  # ------------------------------------------------------------
  # Genotype QC functions
  # ------------------------------------------------------------




  geno_all<- qc_filter_geno(
    X = geno_all,
    hwe_thresh = hwe_thresh,
    maf_min = min_maf
  )$X

  # Freeze the successful baseline cohort for subsequent EM updates.
  if (is.null(previous_result)) {
    all_tissues <- intersect(all_tissues, tissue_priors$tissue)
  } else {
    all_tissues <- names(previous_result)
  }

  analyze_tissue <- function(target_tissue) {

    cat("\n=====================================\n")
    cat("Tissue:", target_tissue, "\n")
    cat("=====================================\n")
    pheno <- subset(
      pheno_gene,
      SMTS == target_tissue
    )

    ids <- intersect(
      pheno$SUBJID,
      rownames(geno_all)
    )

    rows <- match(
      ids,
      pheno$SUBJID
    )

    # Subset and order the phenotype data before normalization.
    pheno <- pheno[
      rows,
      ,
      drop = FALSE
    ]

    stopifnot(
      identical(
        as.character(pheno$SUBJID),
        as.character(ids)
      )
    )

    # Remove samples with invalid expression or library-size values.
    keep_sample <- (
      is.finite(pheno$count) &
        is.finite(pheno$total_read_count) &
        pheno$total_read_count > 0 &
        !is.na(pheno$SEX)
    )

    pheno <- pheno[
      keep_sample,
      ,
      drop = FALSE
    ]

    # Update IDs after sample filtering so genotypes remain aligned.
    ids <- as.character(
      pheno$SUBJID
    )

    # Skip before normalization and matrix operations if too few matched
    # samples remain in this tissue.
    if (nrow(pheno) < min_samples) {

      cat(
        "Skipping",
        target_tissue,
        "- only",
        nrow(pheno),
        "samples.\n"
      )

      return(NULL)
    }

    median_read <- median(
      pheno$count
    )

    mean_read <- mean(
      pheno$count
    )

    pheno <- transform(
      pheno,
      SMGEBTCHT = factor(SMGEBTCHT)
    )

    # Library-size normalization within the tissue.
    pheno$library_size_factor <- (
      pheno$total_read_count /
        mean(pheno$total_read_count)
    )

    pheno$normalized_expression <- log1p(
      pheno$count /
        pheno$library_size_factor
    )
    x=pheno$normalized_expression
    pheno$normalized_expression =qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
    # Residualize normalized expression on sex.
    pheno$y <- resid(
      lm(
        normalized_expression ~ SEX,
        data = pheno
      )
    )
    x=pheno$y
    pheno$y =qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))

    y_variance <- var(
      pheno$y
    )

    if (
      !is.finite(y_variance) ||
      y_variance <= 1e-12
    ) {

      cat(
        "Skipping",
        target_tissue,
        "- residual phenotype variance is zero.\n"
      )

      return(NULL)
    }


    # Retain original hard calls. Rare genotype classes are handled by
    # min_obs inside the fitter, which forces delta=0 for that SNP.
    geno <- geno_all[ids, , drop = FALSE]
    keep <- colSums(geno) >= min_n_rec & matrixStats::colSds(geno) > 0
    geno <- geno[, keep, drop = FALSE]
    if (!ncol(geno)) return(NULL)
    stopifnot(identical(as.character(pheno$SUBJID), rownames(geno)))

    prior <- spe_prior(tissue_priors, target_tissue)
    old <- previous_result[[target_tissue]]
    warm <- old$fit_slide_prior
    if (!is.null(previous_result)) {
      if (is.null(warm) || !identical(old$sample_ids, ids) ||
          !identical(colnames(warm$alpha), colnames(geno)))
        em_stop_fit(paste("Slider warm-start samples/SNPs changed:", target_tissue))
    }
    args <- list(X = geno, y = pheno$y, delta_grid = delta_grid,
                 delta_prior = prior, min_obs = min_obs, L = L,
                 standardize = standardize, estimate_prior_method = estimate_prior_method,
                 max_iter = max_iter, tol = tol, min_abs_corr = min_abs_corr,
                 verbose = verbose)
    if (!is.null(warm)) args$model_init <- warm
    fit <- do.call(susieRSlidePrior::susie, args)
    tryCatch(spe_fit_counts(fit, prior, require_converged = FALSE),
             error = function(e) em_stop_fit(conditionMessage(e)))
    predictor_map <- data.frame(predictor_index = seq_len(ncol(geno)),
                                 predictor_name = colnames(geno), snp = colnames(geno),
                                 coding = "slider")
    # Save an additive marginal association screen for later summaries only.
    # It never determines which fits/components enter the prior update.
    r <- as.numeric(cor(geno, pheno$y))
    stat <- abs(r) * sqrt((length(pheno$y) - 2) / pmax(0, 1 - r^2))
    min_pv <- min(2 * pt(stat, df = length(pheno$y) - 2, lower.tail = FALSE))
    list(fit_slide_prior = fit, delta_prior = prior,
         slide_prior_lead_snp_tss_distance =
           get_cs_lead_tss_distance(fit, predictor_map, tss),
         predictor_map = predictor_map, sample_ids = ids,
         em_warm_started = !is.null(warm), n_SNP = ncol(geno),
         n_ind = length(pheno$y), mean_phe = mean(pheno$y),
         median_phe = median(pheno$y), median_read = median_read,
         mean_read = mean_read, min_pv = min_pv)
  }

  fits <- list()
  tissue_errors <- if (is.null(previous_result)) list() else
    attr(previous_result, "tissue_errors", exact = TRUE)

  for (target_tissue in all_tissues) {

    tissue_result <- tryCatch(
      analyze_tissue(
        target_tissue
      ),
      error = function(e) e
    )

    if (inherits(tissue_result, "error")) {

      if (inherits(tissue_result, "em_fit_error")) stop(tissue_result)
      if (!is.null(previous_result))
        em_stop_fit(paste("Previously fitted tissue failed:", target_tissue,
                          conditionMessage(tissue_result)))

      error_message <- conditionMessage(
        tissue_result
      )

      message(
        "ERROR [",
        target_gene,
        " / ",
        target_tissue,
        "]: ",
        error_message
      )

      tissue_errors[[target_tissue]] <- error_message
      next
    }

    if (is.null(tissue_result)) {
      if (!is.null(previous_result))
        em_stop_fit(paste("Previously fitted tissue no longer passes data QC:", target_tissue))
      next
    }

    fits[[target_tissue]] <- tissue_result
  }

  # Preserve tissue-level failures for auditing without changing the
  # existing list-of-successful-tissues result structure.
  attr(fits, "tissue_errors") <- tissue_errors

  fits
}
