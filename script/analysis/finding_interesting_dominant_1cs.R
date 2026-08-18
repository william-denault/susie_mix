
load("/project2/mstephens/wdenault/susie_mix/res_summary.RData")

idx= which(res_summary$min_pv< 5e-8)


plot_folder= "/project2/mstephens/wdenault/susie_mix/plot/one_cs_dominant"
res_idx= res_summary[idx,]
res_idx= res_idx[which(res_idx$mean_count>100),]
idx= which(res_summary$ncs_susie>0)

res_1cs=  res_idx[which(res_idx$ncs_susie==1 & res_idx$ncs_susie_mix==1 & res_idx$overlap==0),]
res_1cs[ which(res_1cs$n_dom== 1),]

int_df= res_1cs[ which(res_1cs$n_dom== 1),]


min_n_rec=5



min_n_rec <- 5

for (l in 1:nrow(int_df)) {

  target_gene <- int_df$gene[l]
  tissue <- int_df$tissue[l]

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

  # Read in the covariate data.
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

  cov2 <- transform(cov2, SUBJID = substr(SAMPID, 1, 10))
  cov <- merge(cov1, cov2, by = "SUBJID")

  cov <- cov[
    c(
      "SUBJID", "SAMPID", "SEX", "AGE", "SMTS", "SMTSD",
      "SMGEBTCHT", "SMAFRZE"
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

  # Read in the gene expression data.
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

  # Align the gene expression and covariate data.
  ids <- intersect(cov$SAMPID, rownames(pheno_all))
  rows1 <- match(ids, cov$SAMPID)
  rows2 <- match(ids, rownames(pheno_all))

  cov <- cov[rows1, ]
  pheno_all <- pheno_all[rows2, ]

  # Extract the gene expression data for the target gene.
  j <- which(gene_info$Description == target_gene)
  pheno_gene <- cbind(cov, data.frame(count = pheno_all[, j]))

  # All tissues available for this gene.
  all_tissues <- levels(droplevels(pheno_gene$SMTS))

  cat(
    "Found", length(all_tissues),
    "tissues for gene", target_gene, ":\n"
  )

  print(all_tissues)

  # Select the SNPs for the target gene.
  genes <- get_gene_annotations(gtf_file)
  genes <- subset(genes, gene_name == target_gene)

  chr <- as.numeric(substr(genes$chromosome, 4, 5))
  tss <- with(genes, ifelse(strand == "+", start, end))
  pos0 <- tss - cis_window
  pos1 <- tss + cis_window

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

  geno_file_raw <- paste0(plink_out_prefix, ".raw")

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
  storage.mode(geno_all) <- "double"

  # Remove SNPs with missing genotypes.
  x <- colSums(is.na(geno_all))
  j <- which(x == 0)
  geno_all <- geno_all[, j, drop = FALSE]

  # Remove SNPs that do not vary.
  x <- colSds(geno_all)
  j <- which(x > 0)
  geno_all <- geno_all[, j, drop = FALSE]

  recode_matrix_by_freq <- function(X) {

    storage.mode(X) <- "integer"

    counts <- rbind(
      colSums(X == 0L),
      colSums(X == 1L),
      colSums(X == 2L)
    )

    new_code <- apply(
      counts,
      2,
      function(cnt) rank(-cnt, ties.method = "first") - 1L
    )

    out <- X

    for (v in 0:2) {
      idx <- which(X == v, arr.ind = TRUE)
      out[idx] <- new_code[v + 1L, idx[, "col"]]
    }

    dimnames(out) <- dimnames(X)
    out
  }

  qc_recode_geno <- function(
    X,
    hwe_thresh = 1e-5,
    maf_min = 0.05) {

    storage.mode(X) <- "integer"

    # 1. Orient to the minor allele.
    af <- colMeans(X) / 2
    flip <- af > 0.5

    if (any(flip)) {
      X[, flip] <- 2L - X[, flip, drop = FALSE]
    }

    # 2. Filter by MAF.
    X <- X[
      ,
      which(colMeans(X) / 2 > maf_min),
      drop = FALSE
    ]

    n <- nrow(X)
    af <- colMeans(X) / 2
    maf <- af

    # 3. HWE chi-square test.
    count0 <- colSums(X == 0L)
    count1 <- colSums(X == 1L)
    count2 <- colSums(X == 2L)

    exp0 <- n * (1 - maf)^2
    exp1 <- n * 2 * maf * (1 - maf)
    exp2 <- n * maf^2

    chisq <- (
      (count0 - exp0)^2 / exp0 +
        (count1 - exp1)^2 / exp1 +
        (count2 - exp2)^2 / exp2
    )

    hwe_p <- pchisq(
      chisq,
      df = 1,
      lower.tail = FALSE
    )

    # Recode columns that fail HWE.
    fail <- which(
      !is.na(hwe_p) &
        hwe_p < hwe_thresh
    )

    if (length(fail) > 0) {
      cat(
        length(fail),
        "SNP(s) failed HWE (p <",
        hwe_thresh,
        ") - recoding by frequency:\n"
      )

      print(colnames(X)[fail])

      X[, fail] <- recode_matrix_by_freq(
        X[, fail, drop = FALSE]
      )
    }

    list(
      X = X,
      maf = maf,
      hwe_p = hwe_p,
      flipped = which(flip),
      recoded = fail
    )
  }

  geno_all <- qc_recode_geno(
    X = geno_all,
    hwe_thresh = hwe_thresh,
    maf_min = min_maf
  )$X

  pos <- as.numeric(
    sapply(
      strsplit(colnames(geno_all), "_"),
      "[[",
      2
    )
  ) / 1e6

  recode_snp_matrix <- function(X, warn = TRUE) {

    X <- as.matrix(X)

    if (warn && !all(X %in% c(0L, 1L, 2L))) {
      warning(
        paste0(
          "recode_snp_matrix: entries are not all in {0,1,2}; ",
          "expecting additive 0/1/2 dosages."
        ),
        call. = FALSE
      )
    }

    storage.mode(X) <- "integer"

    dominant <- (X >= 1L) * 1L
    recessive <- (X == 2L) * 1L

    dimnames(dominant) <- dimnames(X)
    dimnames(recessive) <- dimnames(X)

    list(
      additive = X,
      dominant = dominant,
      recessive = recessive
    )
  }

  geno_mix_all <- recode_snp_matrix(geno_all)

  # Give each mixed predictor a unique name so that its original SNP
  # can still be identified after columns are removed.
  snp_names <- colnames(geno_all)

  colnames(geno_mix_all$additive) <- paste0(
    snp_names,
    "__additive"
  )

  colnames(geno_mix_all$recessive) <- paste0(
    snp_names,
    "__recessive"
  )

  colnames(geno_mix_all$dominant) <- paste0(
    snp_names,
    "__dominant"
  )

  geno_mix_all <- cbind(
    geno_mix_all$additive,
    geno_mix_all$recessive,
    geno_mix_all$dominant
  )

  # Convert to double without removing row names or column names.
  storage.mode(geno_all) <- "double"
  storage.mode(geno_mix_all) <- "double"

  target_tissue <- tissue

  cat("\n=====================================\n")
  cat("Tissue:", target_tissue, "\n")
  cat("=====================================\n")

  # Extract the data for this tissue.
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

  median_read <- median(pheno$count[rows])
  mean_read <- mean(pheno$count[rows])

  pheno <- transform(
    pheno,
    SMGEBTCHT = factor(SMGEBTCHT)
  )

  pheno$y <- resid(
    lm(count ~ SEX, pheno)
  )

  # Align genotype and phenotype data.
  pheno <- pheno[rows, ]

  # Keep the original 0/1/2 genotypes for displaying counts.
  geno_for_counts <- geno_all[ids, , drop = FALSE]

  geno <- geno_all[ids, , drop = FALSE]
  geno_mix <- geno_mix_all[ids, , drop = FALSE]

  # Remove additive predictors with too few minor alleles.
  pb_col_insample <- which(
    colSums(geno) < min_n_rec
  )

  length(pb_col_insample)
  dim(geno)

  if (length(pb_col_insample) > 0) {
    geno <- geno[
      ,
      -pb_col_insample,
      drop = FALSE
    ]
  }

  # Remove mixed predictors with too few observations.
  pb_col <- which(
    colSums(geno_mix) < min_n_rec
  )

  length(pb_col)

  if (length(pb_col) > 0) {

    removed_names <- colnames(geno_mix)[pb_col]

    n_rec_rm <- sum(
      grepl("__recessive$", removed_names)
    )

    n_dom_rm <- sum(
      grepl("__dominant$", removed_names)
    )

    geno_mix <- geno_mix[
      ,
      -pb_col,
      drop = FALSE
    ]

  } else {

    n_rec_rm <- 0
    n_dom_rm <- 0
  }

  out <- readRDS(
    paste0(
      "/project2/mstephens/wdenault/susie_mix/results/",
      target_gene,
      ".rds"
    )
  )

  k <- which(names(out) == target_tissue)

  fit <- susie(
    geno,
    pheno$y,
    L = L,
    standardize = standardize,
    estimate_prior_method = estimate_prior_method,
    min_abs_corr = min_abs_corr,
    verbose = verbose
  )

  fit$sets

  fit_mix <- susie(
    geno_mix,
    pheno$y,
    L = L,
    standardize = standardize,
    estimate_prior_method = estimate_prior_method,
    min_abs_corr = min_abs_corr,
    verbose = verbose
  )

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
  SNP <- geno_for_counts[, selected_snp_name]

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
    "):\n  n0 = ",
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
