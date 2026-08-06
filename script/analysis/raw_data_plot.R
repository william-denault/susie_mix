

res_1cs[order(res_1cs$n_ind, decreasing = TRUE), c(  "tissue",    "gene" , "n_ind", "n_add", "n_rec", "n_dom" )]


target_gene="AKR1C2"
tissue =  "Skin"
target_gene="ASB14"
tissue =  "Skin"
target_gene="ARPC2"
tissue =  "Muscle"


# --- paths ---
datadir          = "/project2/mstephens/gtex"
plink_exec       = file.path(datadir,"plink2")
gene_annot_fun   = "/project2/mstephens/wdenault/susie_mix/script/scan_tissue_attempt/get_gene_annotations.R"
gtf_file         = file.path(datadir,
                             "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
geno_file        = file.path(datadir,
                             "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv")
subject_pheno_file = file.path(datadir,
                               "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz")
sample_attr_file   = file.path(datadir,
                               "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz")
expr_file          = file.path(datadir,
                               "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")

# --- analysis parameters ---
min_maf           = 0.05 # I found problem in the coding better load everything then filter
cis_window        = 5e5
min_samples       = 50
seed              = 1

# --- susie parameters ---
L                     = 10
standardize           = FALSE
estimate_prior_method = "EM"
min_abs_corr          = 0.0
verbose               = FALSE

# --- misc ---
plink_out_prefix  = tempfile("plink_")  # avoids collisions across parallel calls


  library(tools)
  library(data.table)
  library(matrixStats)
  library(susieR)
  source(gene_annot_fun)

  set.seed(seed)

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

  # Read in the genotype data for SNPs near the target gene (done once,
  # reused across tissues since the cis-window doesn't change).
  cat("Extracting genotype data from PLINK file.\n")
  plink_call <- sprintf(paste("%s --bfile %s --chr %d --from-bp %d --to-bp %d",
                              "--snps-only --max-alleles 2 --rm-dup exclude-all",
                              "--threads 2 --memory 8000 --maf %g",
                              "--recode A --out %s"),
                        plink_exec,geno_file,chr,pos0,pos1,min_maf,
                        plink_out_prefix)
  system(plink_call)
  geno_file_raw <- paste0(plink_out_prefix,".raw")
  geno_all <- fread(geno_file_raw,sep = "\t",header = TRUE)
  class(geno_all) <- "data.frame"
  ids <- geno_all$IID
  rownames(geno_all) <- ids
  geno_all <- geno_all[,-(1:6)]
  geno_all <- as.matrix(geno_all)
  storage.mode(geno_all) <- "double"

  # Remove SNPs with missing genotypes (across all individuals, done once).
  x <- colSums(is.na(geno_all))
  j <- which(x == 0)
  geno_all <- geno_all[,j]

  # Remove SNPs that do not vary (across all individuals, done once).
  x <- colSds(geno_all)
  j <- which(x > 0)
  geno_all <- geno_all[,j]

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
  geno_mix_all=cbind( geno_mix_all$additive,
                      geno_mix_all$recessive,
                      geno_mix_all$dominant)

  geno_mix_all=  matrix(as.double(geno_mix_all),
                        ncol=ncol(geno_mix_all),
                        nrow = nrow(geno_mix_all))
  row.names(geno_mix_all)= row.names(geno_all)
  # Store per-tissue results.
  fits <- list()
  #target_tissue = all_tissues[1]
 target_tissue  =tissue

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
    print(all(pheno$SUBJID == rownames(geno))) # Sanity check.

    # Skip tissues with too few samples for a meaningful fine-mapping fit.

    out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",
                         target_gene,
                         ".rds"))
     k = which ( names(out)== target_tissue)


     fit <- susie(geno,pheno$y,L = L,standardize = standardize,
                  estimate_prior_method = estimate_prior_method,
                  min_abs_corr = min_abs_corr,verbose = verbose)

     fit$sets
     fit_mix <- susie(geno_mix,pheno$y,L = L,standardize = standardize,
                      estimate_prior_method = estimate_prior_method,
                      min_abs_corr = min_abs_corr,verbose = verbose)

     fit_mix$sets
     ncol(geno)+fit$sets$cs$L1[1]
SNP=geno[ , out[[k]]$susie_add$sets$cs$L1[1]]
     boxplot(pheno$y~ SNP, main= paste(target_gene, "Expression level vs SNP fine mapped\n n = " , length(pheno$y)))


     table(geno[ , out[[k]]$susie_add$sets$cs$L1[1]] )

