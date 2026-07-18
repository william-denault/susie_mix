# SuSiE analysis of GTEx data from this paper:
# https://doi.org/10.1126/science.aaz1776
#
# Note that Some of the files used in this analysis were downloaded
# from the GTEx Portal:
# https://gtexportal.org/home/downloads/adult-gtex/metadata
#
# sinteractive ...
# module load R/4.2.0
# R
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_20"
#
pt=proc.time()

## Gene names in the vcf file-----
#GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
library(tools)
library(data.table)
library(matrixStats)
library(susieR) # 0.14.2
library(ggplot2)
library(cowplot)
source("/project2/mstephens/wdenault/susie_mix/script/pc_script/get_gene_annotations.R")
set.seed(1)
 datadir     <- "/project2/mstephens/gtex"
#datadir       <- "/scratch/midway2/pcarbo/gtex"
plink_exec    <- file.path(datadir,"plink2")
target_gene   <- "CBX8"
target_tissue <- "Lung"#
min_maf       <- 0.05

# Read in the covariate data.
cat("Importing covariate data.\n")
cov1 <- read.table(file.path(datadir,
                             "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz"),
                   header = TRUE,sep = "\t",stringsAsFactors = FALSE)
cov2 <- read.table(file.path(datadir,
                             "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz"),
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
pheno <- fread(file.path(datadir,
                         "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"),
               sep = "\t",skip = 2,header = TRUE,showProgress = TRUE)
class(pheno) <- "data.frame"
gene_info <- pheno[1:2]
pheno <- pheno[-(1:2)]
pheno <- as.matrix(pheno)
pheno <- t(pheno)
storage.mode(pheno) <- "double"
colnames(pheno) <- gene_info$Name

# Align the gene expression and covariate data.
ids   <- intersect(cov$SAMPID,rownames(pheno))
rows1 <- match(ids,cov$SAMPID)
rows2 <- match(ids,rownames(pheno))
cov   <- cov[rows1,]
pheno <- pheno[rows2,]

# Extract the gene expression data for the target gene.
j <- which(gene_info$Description == target_gene)
pheno <- cbind(cov,data.frame(count = pheno[,j]))

# Extract the data for the target tissue.
pheno <- subset(pheno,SMTS == target_tissue)
#TODO TISSUE names -------
pheno <- transform(pheno,SMGEBTCHT = factor(SMGEBTCHT))

# From Sec. 4.1 of the GTEx supplement: "We concluded that 5 PCs is a
# good choice that controls for population structure reasonably well
# while avoiding reduction of power in smaller tissues. Additionally,
# WGS sequencing platform (HiSeq 2000 or HiSeq X), WGS library
# construction protocol (PCR-based or PCR-free) and donor sex were
# included in the set of covariates used in the association analyses."
pheno$y <- resid(lm(count ~ SEX,pheno))

# Select the SNPs for the target gene.
gene_file <-
  file.path(datadir,
            "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
genes <- get_gene_annotations(gene_file)
genes <- subset(genes,gene_name == target_gene)
chr   <- as.numeric(substr(genes$chromosome,4,5))
tss   <- with(genes,ifelse(strand == "+",start,end))
pos0  <- tss - 5e5
pos1  <- tss + 5e5

# Read in the genotype data for SNPs near the target gene.
cat("Extracting genotype data from PLINK file.\n")
geno_file <- file.path(datadir,
                       "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv")
plink_call <- sprintf(paste("%s --bfile %s --chr %d --from-bp %d --to-bp %d",
                            "--snps-only --max-alleles 2 --rm-dup exclude-all",
                            "--threads 2 --memory 8000 --maf %g",
                            "--recode A"),
                      plink_exec,geno_file,chr,pos0,pos1,min_maf)
system(plink_call)
geno_file_raw <- file.path(datadir,"plink2.raw")
geno <- fread("plink2.raw",sep = "\t",header = TRUE)
class(geno) <- "data.frame"
ids <- geno$IID
rownames(geno) <- ids
geno <- geno[,-(1:6)]
geno <- as.matrix(geno)
storage.mode(geno) <- "double"
ids   <- intersect(pheno$SUBJID,rownames(geno))
rows  <- match(ids,pheno$SUBJID)
pheno <- pheno[rows,]
geno  <- geno[ids,]
print(all(pheno$SUBJID == rownames(geno))) # Sanity check.

# Remove SNPs with missing genotypes.
x <- colSums(is.na(geno))
j <- which(x == 0)
geno <- geno[,j]

# Remove SNPs that do not vary.
x    <- colSds(geno)
j    <- which(x > 0)
geno <- geno[,j]

# Run susie.
cat("Running susie.\n")
fit <- susie(geno,pheno$y,L = 10,standardize = FALSE,
             estimate_prior_method = "EM",min_abs_corr = 0,
             verbose = TRUE)
fit$sets
x <- sapply(fit$sets$cs,length)
i <- which(x <= 10)
pos <- as.numeric(sapply(strsplit(colnames(geno),"_"),"[[",2))/1e6
print(fit$sets$cs[i])

# PIP plot.
# For CBX8 expression in lung, compare to Fig. 3 of the GTEx paper.
pdat <- data.frame(pos = pos,pip = fit$pip,cs = "none")
n <- length(fit$sets$cs)
for (i in names(fit$sets$cs)) {
  j <- fit$sets$cs[[i]]
  if (length(j) <= 10)
    pdat[j,"cs"] <- i
}
pdat <- transform(pdat,cs = factor(cs,c("none",names(fit$sets$cs))))
cs_colors <- c("black","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
               "#ffff33","#a65628","#f781bf")
p1 <- ggplot(pdat,aes(x = pos,y = pip,color = cs)) +
  geom_point() +
  scale_color_manual(values = cs_colors) +
  labs(x = "base-pair position (Mb)",y = "PIP") +
  theme_cowplot(font_size = 12)
print(p1)


proc.time()-pt
