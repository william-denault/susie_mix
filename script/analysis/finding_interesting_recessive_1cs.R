
load("/project2/mstephens/wdenault/susie_mix/res_summary.RData")

idx= which(res_summary$min_pv< 5e-8)


plot_folder= "/project2/mstephens/wdenault/susie_mix/plot/one_cs_recessive"
res_idx= res_summary[idx,]
res_idx= res_idx[which(res_idx$mean_count>100),]
idx= which(res_summary$ncs_susie>0)

res_1cs=  res_idx[which(res_idx$ncs_susie==1 & res_idx$ncs_susie_mix==1  ),]
res_1cs[ which(res_1cs$n_dom== 1),]

int_df= res_1cs[ which(res_1cs$n_rec== 1),]



min_n_rec <- 5


source("/project2/mstephens/wdenault/susie_mix/script/scan_tissue_attempt/workhorse_utils.R")
for (l in 1:nrow(int_df)) {

  target_gene <- int_df$gene[l]
  tissue <- int_df$tissue[l]
  target_tissue=tissue
  # --- paths ---
  datadir <- "/project2/mstephens/gtex"
  plink_exec <- file.path(datadir, "plink2")
  gene_annot_fun <- paste0(
    "/project2/mstephens/wdenault/susie_mix/",
    "script/scan_tissue_attempt/get_gene_annotations.R"
  )

  gtf_file <- file.path(
    datadir,
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"
  )

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

  # --- analysis parameters ---
  min_maf_plink <- 0.00
  min_maf <- 0.05
  hwe_thresh <- 1e-5
  cis_window <- 5e5
  min_samples <- 50
  seed <- 1

  # --- susie parameters ---
  L <- 10
  standardize <- FALSE
  estimate_prior_method <- "EM"
  min_abs_corr <- 0.0
  verbose <- FALSE

  # --- misc ---
  temp_dir <- "/project2/mstephens/wdenault/susie_mix/temp_plink/"

  set.seed(seed)

  plink_out_prefix <- paste0(temp_dir, "plink_", target_gene)
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

  pheno$normalized_expression <- log1p(
    pheno$count /
      pheno$library_size_factor
  )

  # Residualize normalized expression on sex.
  pheno$y <- resid(
    lm(
      normalized_expression ~ SEX,
      data = pheno
    )
  )



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


  fit_mix <- susie(
    geno_mix,
    pheno$y,
    L = L,
    standardize = standardize,
    estimate_prior_method = estimate_prior_method,
    min_abs_corr = min_abs_corr,
    verbose = verbose
  )

  out <- readRDS(
    paste0(
      "/project2/mstephens/wdenault/susie_mix/results/",
      target_gene,
      ".rds"
    )
  )

  k <- which(names(out) == target_tissue)

  fit$sets



  fit_mix$sets

  # Select the predictor from the MIXED-CODING credible set.
  selected_mix_column <- fit_mix$sets$cs[[1]][1]
  selected_mix_name <- colnames(geno_mix)[selected_mix_column]

  # Recover the original SNP and its coding.
  selected_snp_name <- sub(
    "__(additive|recessive|dominant)$",
    "",
    selected_mix_name
  )

  selected_coding <- sub(
    "^.*__",
    "",
    selected_mix_name
  )

  # Use the original 0/1/2 genotype to calculate the displayed counts.
  SNP <- geno[, selected_snp_name]

  number <- table(
    factor(
      SNP,
      levels = 0:2
    )
  )

  my_main <- paste0(
    target_gene,
    " expression level in ",
    tissue,
    "\n",
    selected_snp_name,
    " (",
    selected_coding,
    "):\n n0 = ",
    number["0"],
    ", n1 = ",
    number["1"],
    ", n2 = ",
    number["2"]
  )

  print(
    boxplot(
      pheno$y ~ SNP,
      main = my_main
    )
  )

  png(
    file.path(
      plot_folder,
      paste0(
        target_gene,
        "_",
        tissue,
        "_boxplot.png"
      )
    ),
    width = 2000,
    height = 2000,
    res = 300
  )

  boxplot(
    pheno$y ~ SNP,
    main = my_main
  )

  dev.off()
}
