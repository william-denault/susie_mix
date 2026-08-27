# run_susie_gene.R
#
# Run additive and mixed-coding SuSiE fine-mapping for all tissues
# available for one gene.


source("/project2/mstephens/wdenault/susie_mix/script/scan_tissue_attempt/workhorse_utils.R")
run_susie_gene <- function(
    target_gene = "GTF2H2",

    # --- paths ---
    datadir = "/project2/mstephens/gtex",
    plink_exec = file.path(datadir, "plink2"),
    gene_annot_fun = paste0(
      "/project2/mstephens/wdenault/susie_mix/",
      "script/scan_tissue_attempt/get_gene_annotations.R"
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

    # --- SuSiE parameters ---
    L = 10,
    standardize = FALSE,
    estimate_prior_method = "EM",
    min_abs_corr = 0.0,
    verbose = FALSE,

    # --- misc ---
    temp_dir = "/project2/mstephens/wdenault/susie_mix/temp_plink/"
) {


  library(tools)
  library(data.table)
  library(matrixStats)
  library(susieR)

  source(gene_annot_fun)

  set.seed(seed)

  plink_out_prefix <- paste0(
    temp_dir,
    "plink_",
    target_gene
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

  chr <- as.numeric(
    substr(genes$chromosome[1], 4, 5)
  )

  tss <- with(
    genes[1, ],
    ifelse(strand == "+", start, end)
  )

  pos0 <- max(0, tss - cis_window)
  pos1 <- tss + cis_window

  # ------------------------------------------------------------
  # Extract cis-genotypes
  # ------------------------------------------------------------

  if (!file.exists(paste0(plink_out_prefix, ".raw"))) {

    cat("Extracting genotype data from PLINK file.\n")

    plink_call <- sprintf(
      paste(
        "%s --bfile %s --chr %d --from-bp %d --to-bp %d",
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

    system(plink_call)
  }

  geno_file_raw <- paste0(
    plink_out_prefix,
    ".raw"
  )
  geno_file_raw <- paste0(
    plink_out_prefix,
    ".raw"
  )
  hwe_thresh = 1e-8
  maf_min = 0.05
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

  # ------------------------------------------------------------
  # Construct additive, recessive and dominant predictors
  # ------------------------------------------------------------



  geno_mix_parts <- recode_snp_matrix(
    geno_all
  )

  base_snp_names <- colnames(geno_all)

  # Give every mixed predictor a unique name.
  colnames(geno_mix_parts$additive) <- paste0(
    base_snp_names,
    "__additive"
  )

  colnames(geno_mix_parts$recessive) <- paste0(
    base_snp_names,
    "__recessive"
  )

  colnames(geno_mix_parts$dominant) <- paste0(
    base_snp_names,
    "__dominant"
  )

  geno_mix_all <- cbind(
    geno_mix_parts$additive,
    geno_mix_parts$recessive,
    geno_mix_parts$dominant
  )

  # Map every mixed predictor to its biological SNP and coding.
  geno_mix_snp_all <- rep(
    base_snp_names,
    times = 3L
  )

  geno_mix_coding_all <- rep(
    c(
      "additive",
      "recessive",
      "dominant"
    ),
    each = length(base_snp_names)
  )

  # Convert to double without losing names.
  storage.mode(geno_all) <- "double"
  storage.mode(geno_mix_all) <- "double"

  fits <- list()
  for (target_tissue in all_tissues) {

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

    qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
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


    geno <- geno_all[
      ids,
      ,
      drop = FALSE
    ]

    geno_mix <- geno_mix_all[
      ids,
      ,
      drop = FALSE
    ]

    # ----------------------------------------------------------
    # Tissue-specific predictor filtering
    # ----------------------------------------------------------

    # Filter the ordinary additive matrix.
    keep_add <- colSums(geno) >= min_n_rec

    n_add_rm <- sum(!keep_add)

    geno <- geno[
      ,
      keep_add,
      drop = FALSE
    ]

    add_snp_names <- colnames(geno)

    # Filter the mixed matrix and its predictor maps together.
    keep_mix <- colSums(geno_mix) >= min_n_rec

    n_rec_rm <- sum(
      !keep_mix &
        geno_mix_coding_all == "recessive"
    )

    n_dom_rm <- sum(
      !keep_mix &
        geno_mix_coding_all == "dominant"
    )

    geno_mix <- geno_mix[
      ,
      keep_mix,
      drop = FALSE
    ]

    mix_snp_names <- geno_mix_snp_all[
      keep_mix
    ]

    mix_coding <- geno_mix_coding_all[
      keep_mix
    ]

    mix_predictor_names <- colnames(
      geno_mix
    )

    # Convenient data frames for downstream processing.
    add_predictor_map <- data.frame(
      predictor_index = seq_len(ncol(geno)),
      predictor_name = colnames(geno),
      snp = add_snp_names,
      coding = rep("additive", ncol(geno)),
      stringsAsFactors = FALSE
    )

    mix_predictor_map <- data.frame(
      predictor_index = seq_len(ncol(geno_mix)),
      predictor_name = mix_predictor_names,
      snp = mix_snp_names,
      coding = mix_coding,
      stringsAsFactors = FALSE
    )

    print(
      all(
        pheno$SUBJID == rownames(geno)
      )
    )

    # Skip tissues with too few samples.
    if (nrow(geno) < min_samples) {

      cat(
        "Skipping",
        target_tissue,
        "- only",
        nrow(geno),
        "samples.\n"
      )

      next
    }

    # Skip if filtering removed all predictors.
    if (
      ncol(geno) == 0L ||
      ncol(geno_mix) == 0L
    ) {

      cat(
        "Skipping",
        target_tissue,
        "- no predictors remained after filtering.\n"
      )

      next
    }

    # ----------------------------------------------------------
    # Marginal association analysis
    # ----------------------------------------------------------

    cat("Running SuSiE.\n")

    set.seed(seed)

    perm_y <- sample(pheno$y)

    pv <- rep(
      1,
      ncol(geno)
    )

    for (k in seq_len(ncol(geno))) {

      if (var(geno[, k]) > 0) {

        pv[k] <- summary(
          lm(pheno$y ~ geno[, k])
        )$coefficients[2, 4]
      }
    }

    # ----------------------------------------------------------
    # SuSiE fits
    # ----------------------------------------------------------

    fit <- susie(
      geno,
      pheno$y,
      L = L,
      standardize = standardize,
      estimate_prior_method = estimate_prior_method,
      min_abs_corr = min_abs_corr,
      verbose = verbose
    )

    fit_perm <- susie(
      geno,
      perm_y,
      L = L,
      standardize = standardize,
      estimate_prior_method = estimate_prior_method,
      min_abs_corr = min_abs_corr,
      verbose = verbose
    )

    fit_mix <- susie(
      geno_mix,
      pheno$y,
      L = L,
      standardize = standardize,
      estimate_prior_method = estimate_prior_method,
      min_abs_corr = min_abs_corr,
      verbose = verbose
    )

    fit_mix_perm <- susie(
      geno_mix,
      perm_y,
      L = L,
      standardize = standardize,
      estimate_prior_method = estimate_prior_method,
      min_abs_corr = min_abs_corr,
      verbose = verbose
    )



    # ----------------------------------------------------------
    # Lead SNP and distance to the TSS for each credible set
    # ----------------------------------------------------------

    susie_add_lead_snp_tss_distance <- (
      get_cs_lead_tss_distance(
        fit = fit,
        predictor_map = add_predictor_map,
        tss = tss
      )
    )

    susie_mix_lead_snp_tss_distance <- (
      get_cs_lead_tss_distance(
        fit = fit_mix,
        predictor_map = mix_predictor_map,
        tss = tss
      )
    )
    # ----------------------------------------------------------
    # Store results and predictor mappings
    # ----------------------------------------------------------

    fits[[target_tissue]] <- list(
      susie_add = fit,
      susie_add_perm = fit_perm,
      susie_mix = fit_mix,
      susie_mix_perm = fit_mix_perm,

      # Lead biological SNP and absolute distance to the TSS.
      susie_add_lead_snp_tss_distance =
        susie_add_lead_snp_tss_distance,

      susie_mix_lead_snp_tss_distance =
        susie_mix_lead_snp_tss_distance,

      # Predictor maps.
      add_predictor_map = add_predictor_map,
      mix_predictor_map = mix_predictor_map,

      add_snp_names = add_snp_names,
      mix_snp_names = mix_snp_names,
      mix_coding = mix_coding,
      mix_predictor_names = mix_predictor_names,

      n_SNP = ncol(geno),
      n_mix_predictor = ncol(geno_mix),
      n_add_rm = n_add_rm,
      n_rec_rm = n_rec_rm,
      n_dom_rm = n_dom_rm,

      n_ind = length(perm_y),
      mean_phe = mean(pheno$y),
      median_phe = median(pheno$y),
      min_pv = min(pv),
      median_read = median_read,
      mean_read = mean_read
    )
  }


  # Clean up temporary PLINK files.
  file.remove(
    Sys.glob(
      paste0(plink_out_prefix, ".*")
    )
  )

  fits
}
