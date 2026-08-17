# run_susie_gene.R
#
# Function version of the per-gene, all-tissue SuSiE fine-mapping script.
# Everything that was hardcoded in the original script is now a function
# argument (with the original value kept as the default), and target_gene
# is the primary argument.
#
# Returns: a named list of susie fit objects, one per tissue that had
# enough samples, e.g. fits[["Lung"]]$sets

run_susie_gene <- function(
    target_gene="GTF2H2",#,"CBX8",

    # --- paths ---
    datadir          = "/project2/mstephens/gtex",
    plink_exec       = file.path(datadir,"plink2"),
    gene_annot_fun   = "/project2/mstephens/wdenault/susie_mix/script/scan_tissue_attempt/get_gene_annotations.R",
    gtf_file         = file.path(datadir,
                                 "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"),
    geno_file        = file.path(datadir,
                                 "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"),
    subject_pheno_file = file.path(datadir,
                                   "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz"),
    sample_attr_file   = file.path(datadir,
                                   "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz"),
    expr_file          = file.path(datadir,
                                   "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"),

    # --- analysis parameters ---
    min_maf_plink         = 0.00,
    min_maf           = 0.05,# I found problem in the coding better load everything then filter
    cis_window        = 5e5,
    min_samples       = 50,
    seed              = 1,
    hwe_thresh= 1e-5,
    # --- susie parameters ---
    L                     = 10,
    standardize           = FALSE,
    estimate_prior_method = "EM",
    min_abs_corr          = 0.0,
    verbose               = FALSE,
    min_n_rec=5,

    # --- misc ---
    temp_dir= "/project2/mstephens/wdenault/susie_mix/temp_plink/"
    # avoids collisions across parallel calls
) {

  library(tools)
  library(data.table)
  library(matrixStats)
  library(susieR)
  source(gene_annot_fun)

  set.seed(seed)

  plink_out_prefix  = paste0( temp_dir, "plink_",target_gene)

  # Read in the covariate data.
  cat("Importing covariate data.\n")
  cov1 <- read.table(subject_pheno_file,
                     header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  cov2 <- read.table(sample_attr_file,
                     header = TRUE,sep = "\t",quote = "",
                     stringsAsFactors = FALSE)
  cov2 <- transform(cov2,SUBJID = substr(SAMPID,1,10))
  cov <- merge(cov1,cov2,by = "SUBJID")
  cov <- cov[c("SUBJID","SAMPID","SEX","AGE","SMTS","SMTSD","SMGEBTCHT",
               "SMAFRZE")]
  cov <- subset(cov,SMAFRZE == "RNASEQ")
  cov <- transform(cov,
                   SEX       = SEX - 1,
                   AGE       = factor(AGE),
                   SMTS      = factor(SMTS),
                   SMTSD     = factor(SMTSD),
                   SMGEBTCHT = factor(SMGEBTCHT),
                   SMAFRZE   = factor(SMAFRZE))
  rownames(cov) <- cov$SAMPID

  # Read in the gene expression data.
  cat("Importing gene expression data.\n")
  pheno_all <- fread(expr_file,sep = "\t",skip = 2,header = TRUE,
                     showProgress = TRUE)
  class(pheno_all) <- "data.frame"
  gene_info <- pheno_all[1:2]
  pheno_all <- pheno_all[-(1:2)]
  pheno_all <- as.matrix(pheno_all)
  pheno_all <- t(pheno_all)
  storage.mode(pheno_all) <- "double"
  colnames(pheno_all) <- gene_info$Name

  # Align the gene expression and covariate data.
  ids   <- intersect(cov$SAMPID,rownames(pheno_all))
  rows1 <- match(ids,cov$SAMPID)
  rows2 <- match(ids,rownames(pheno_all))
  cov       <- cov[rows1,]
  pheno_all <- pheno_all[rows2,]

  # Extract the gene expression data for the target gene.
  j <- which(gene_info$Description == target_gene)
  pheno_gene <- cbind(cov,data.frame(count = pheno_all[,j]))

  # All tissues available for this gene.
  all_tissues <- levels(droplevels(pheno_gene$SMTS))
  cat("Found",length(all_tissues),"tissues for gene",target_gene,":\n")
  print(all_tissues)

  # Select the SNPs for the target gene (same cis-window for every tissue).
  genes <- get_gene_annotations(gtf_file)
  genes <- subset(genes,gene_name == target_gene)
  chr   <- as.numeric(substr(genes$chromosome,4,5))
  tss   <- with(genes,ifelse(strand == "+",start,end))
  pos0  <- tss - cis_window
  pos1  <- tss + cis_window


  if (!file.exists (  paste0(plink_out_prefix,".raw"))){
    # Read in the genotype data for SNPs near the target gene (done once,
    # reused across tissues since the cis-window doesn't change).
    cat("Extracting genotype data from PLINK file.\n")
    plink_call <- sprintf(paste("%s --bfile %s --chr %d --from-bp %d --to-bp %d",
                                "--snps-only --max-alleles 2 --rm-dup exclude-all",
                                "--threads 2 --memory 8000 --maf %g",
                                "--recode A --out %s"),
                          plink_exec,geno_file,chr,pos0,pos1,min_maf_plink,
                          plink_out_prefix)
    system(plink_call)
  }


  geno_file_raw <- paste0(plink_out_prefix,".raw")
  geno_all <- fread(geno_file_raw,sep = "\t",header = TRUE)
  class(geno_all) <- "data.frame"
  ids <- geno_all$IID
  rownames(geno_all) <- ids
  geno_all <- geno_all[,-(1:6)]
  geno_all <- as.matrix(geno_all)
  storage.mode(geno_all) <- "double"
  row_name_geno= row.names(geno_all)
  # Remove SNPs with missing genotypes (across all individuals, done once).
  x <- colSums(is.na(geno_all))
  j <- which(x == 0)
  geno_all <- geno_all[,j]

  # Remove SNPs that do not vary (across all individuals, done once).
  x <- colSds(geno_all)
  j <- which(x > 0)

  geno_all <- geno_all[,j]


  recode_matrix_by_freq <- function(X) {
    storage.mode(X) <- "integer"
    counts <- rbind(colSums(X == 0L),
                    colSums(X == 1L),
                    colSums(X == 2L))                 # 3 x ncol(X)
    new_code <- apply(counts, 2, function(cnt) rank(-cnt, ties.method = "first") - 1L)
    out <- X
    for (v in 0:2) {
      idx <- which(X == v, arr.ind = TRUE)
      out[idx] <- new_code[v + 1L, idx[, "col"]]
    }
    dimnames(out) <- dimnames(X)
    out
  }

  qc_recode_geno <- function(X, hwe_thresh = 1e-5, maf_min=0.05) {
    storage.mode(X) <- "integer"

    # 1. orient to minor allele
    af   <- colMeans(X) / 2
    flip <- af > 0.5
    if (any(flip)) X[, flip] <- 2L - X[, flip, drop = FALSE]

    # 2. filter by MAF (drops monomorphic SNPs too, since the comparison is strict >)
    X  <- X[, which(colMeans(X) / 2 > maf_min), drop = FALSE]
    n  <- nrow(X)
    af <- colMeans(X) / 2       # already the MAF post-flip/post-filter — no pmin needed
    maf <- af

    # 3. HWE chi-square test
    count0 <- colSums(X == 0L); count1 <- colSums(X == 1L); count2 <- colSums(X == 2L)
    exp0 <- n * (1 - maf)^2; exp1 <- n * 2 * maf * (1 - maf); exp2 <- n * maf^2
    chisq <- (count0 - exp0)^2 / exp0 + (count1 - exp1)^2 / exp1 + (count2 - exp2)^2 / exp2
    hwe_p <- pchisq(chisq, df = 1, lower.tail = FALSE)

    # 3. frequency-recode only the columns that fail HWE
    fail <- which(!is.na(hwe_p) & hwe_p < hwe_thresh)
    if (length(fail) > 0) {
      cat(length(fail), "SNP(s) failed HWE (p <", hwe_thresh, ") - recoding by frequency:\n")
      print(colnames(X)[fail])
      X[, fail] <- recode_matrix_by_freq(X[, fail, drop = FALSE])
    }

    list(X = X, maf = maf, hwe_p = hwe_p, flipped = which(flip), recoded = fail)
  }

  geno_all <-   qc_recode_geno(X=geno_all,
                               hwe_thresh =  hwe_thresh,
                               maf_min=min_maf)$X

  ### 1 check for allele flip
  ### 2 remove SNP with maf below 5% then
  ### 3 test for HWE and recode by frequency


  pos <- as.numeric(sapply(strsplit(colnames(geno_all),"_"),"[[",2))/1e6
  recode_snp_matrix <- function(X, warn = TRUE) {
    X <- as.matrix(X)
    if (warn && !all(X %in% c(0L, 1L, 2L)))
      warning("recode_snp_matrix: entries are not all in {0,1,2}; ",
              "expecting additive 0/1/2 dosages.", call. = FALSE)

    storage.mode(X) <- "integer"          # additive, unchanged
    dominant  <- (X >= 1L) * 1L           # 1 if 1 or 2, else 0
    recessive <- (X == 2L) * 1L           # 1 if 2,      else 0
    dimnames(dominant) <- dimnames(recessive) <- dimnames(X)

    list(additive = X, dominant = dominant, recessive = recessive)
  }
  geno_mix_all= recode_snp_matrix( geno_all)
  #SNP_rec_enough_point =  which( apply(geno_mix_all$recessive,2,sum)>min_n_rec)
  #n_rec_SNP= length(SNP_rec_enough_point)

 # if( n_rec_SNP==0){
 #    n_rec_SNP= ncol(geno_mix_all$additive)
     geno_mix_all=cbind( geno_mix_all$additive,
                          geno_mix_all$recessive  ,
                          geno_mix_all$dominant)
  #}else{
  #  geno_mix_all=cbind( geno_mix_all$additive,
  #                      geno_mix_all$recessive[,SNP_rec_enough_point]  ,
  #                      geno_mix_all$dominant)
  #}


  geno_all      =  matrix(as.double(geno_all),
                         ncol=ncol(geno_all),
                         nrow = nrow(geno_all))
  row.names(geno_all ) = row_name_geno
  geno_mix_all =  matrix(as.double(geno_mix_all),
                       ncol=ncol(geno_mix_all),
                       nrow = nrow(geno_mix_all))
  row.names(geno_mix_all)= row.names(geno_all)
  # Store per-tissue results.
  fits <- list()
  #target_tissue = all_tissues[1]
  for (target_tissue in all_tissues) {

    cat("\n=====================================\n")
    cat("Tissue:",target_tissue,"\n")
    cat("=====================================\n")

    # Extract the data for this tissue.
    pheno <- subset(pheno_gene,SMTS == target_tissue)

    ids         <- intersect(pheno$SUBJID,rownames(geno_all))
    rows        <- match(ids,pheno$SUBJID)
    median_read <- median( pheno$count[ rows])
    mean_read   <- mean(pheno$count[ rows])

    pheno <- transform(pheno,SMGEBTCHT = factor(SMGEBTCHT))

    # From Sec. 4.1 of the GTEx supplement: "We concluded that 5 PCs is a
    # good choice that controls for population structure reasonably well
    # while avoiding reduction of power in smaller tissues. Additionally,
    # WGS sequencing platform (HiSeq 2000 or HiSeq X), WGS library
    # construction protocol (PCR-based or PCR-free) and donor sex were
    # included in the set of covariates used in the association analyses."
    pheno$y <- resid(lm(count ~ SEX,pheno))

    # Align genotype and phenotype data for this tissue's samples.

    pheno <- pheno[rows,]

    geno  <- geno_all[ids,]


    geno_mix  <- geno_mix_all[ids,]


    geno_mix  <- geno_mix_all[ids,]

    pb_col= which(apply(geno_mix,2,sum)< min_n_rec)
    if ( length(pb_col)>1){
      n_rec_rm= 0
      n_dom_rm=0

      if( length( which(pb_col <2*ncol(geno)+1))>1){
        n_rec_rm=length( which(pb_col <2*ncol(geno)+1))

      }
      if(length( which(pb_col >2*ncol(geno)+1))){
        n_dom_rm=length( which(pb_col <2*ncol(geno)+1))
      }



      geno_mix[ , -pb_col]= geno_mix[ ,- pb_col]
    }else{

      n_rec_rm= 0
      n_dom_rm=0
    }

    print(all(pheno$SUBJID == rownames(geno))) # Sanity check.

    # Skip tissues with too few samples for a meaningful fine-mapping fit.
    if (nrow(geno) < min_samples) {
      cat("Skipping",target_tissue,"- only",nrow(geno),"samples.\n")
      next
    }


    cat("Running susie.\n")
    set.seed(1)
    perm_y=sample(pheno$y)

    pv= rep(1, ncol(geno))
    for ( k in 1: ncol(geno)){

      if(var(geno[,k])>0){

        pv[k]= summary(lm( pheno$y~geno[,k]))$coefficients[2,4]
      }
    }

    fit <- susie(geno,pheno$y,L = L,standardize = standardize,
                 estimate_prior_method = estimate_prior_method,
                 min_abs_corr = min_abs_corr,verbose = verbose)

    fit_perm <- susie(geno,perm_y,L = L,standardize = standardize,
                 estimate_prior_method = estimate_prior_method,
                 min_abs_corr = min_abs_corr,verbose = verbose)

    fit_mix <- susie(geno_mix,pheno$y,L = L,standardize = standardize,
                     estimate_prior_method = estimate_prior_method,
                     min_abs_corr = min_abs_corr,verbose = verbose)

    fit_mix_perm <- susie(geno_mix,perm_y,L = L,standardize = standardize,
                     estimate_prior_method = estimate_prior_method,
                     min_abs_corr = min_abs_corr,verbose = verbose)

    fits[[target_tissue]] <- list(susie_add= fit,
                                  susie_add_perm= fit_perm,

                                  susie_mix=fit_mix,
                                  susie_mix_perm= fit_mix_perm,
                                  n_SNP= ncol(geno),
                                  n_rec_rm= n_rec_rm,
                                  n_dom_rm=n_dom_rm,
                                  n_ind= length(perm_y),
                                  mean_phe= mean(perm_y),

                                  median_phe= median(perm_y),
                                  min_pv=min(pv),
                                  median_read  =median_read,
                                  mean_read = mean_read
                                  )

  }

  # Clean up the temporary plink output files for this call.
  file.remove(Sys.glob(paste0(plink_out_prefix,".*")))

  # fits is a named list keyed by tissue, e.g. fits[["Lung"]]$sets
  fits
}
