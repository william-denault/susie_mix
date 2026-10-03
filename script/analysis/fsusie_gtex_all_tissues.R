# Minimal adaptation of the pasted script: pool all matched RNA tissue samples.
# A donor's genotype row is repeated for each of that donor's RNA samples.
# This does not add tissue covariates or account for within-donor correlation.
# Binning and model settings follow the supplied script.

rm(list=ls())
library(data.table)
library(vcfR)
library(fsusieR)

library(flashier)
library(mvsusieR)
ensembl_id= "ENSG00000112081"#SRSF3    "ENSG00000214941" #ZSWIM7
windows= 200000
n_bins=1024

#loading GTF3C6
lf= list.files("/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/")

id_gene= grep(ensembl_id,lf)
i = id_gene
desired_bins=n_bins
library(data.table)
library(vcfR)
library(fsusieR)
tt = fread(paste0("/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/", lf[i]),
           sep = ",",header = TRUE)
tt[1:10,1:10]



# Assuming 'tt' is your data frame, and the column names are in the format 'chrX:100627058'

# Get the column names of the dataframe
col_names <- colnames(tt)

# Extract the chromosome number and base pair using string manipulation
# This assumes the format is 'chrX:basepair'

# Define a function to extract chromosome and base pair



extract_chromosome_basepair <- function(column_name) {
  parts <- strsplit(column_name, ":")[[1]]  # Split by colon
  if(length(parts) == 2) {  # If it follows the expected format
    chromosome <- sub("^chr", "", parts[1])  # Remove 'chr' prefix from 'chrX' to get 'X'
    basepair <- parts[2]  # '100627058'
    return(c(chromosome, basepair))
  } else {
    return(NA)  # If it doesn't match the format
  }
}

# Apply the function to each column name
extracted_info <- lapply(col_names, extract_chromosome_basepair)

# Convert the list into a data frame for easier viewing
extracted_df <- do.call(rbind, extracted_info)
colnames(extracted_df) <- c("Chromosome", "BasePair")
# Print the extracted information
#print(extracted_df)
extracted_df =data.frame( extracted_df)
extracted_df$BasePair =as.numeric(extracted_df$BasePair)


# Find the minimum and maximum base pair values
min_basepair <- min(extracted_df$BasePair, na.rm = TRUE)
max_basepair <- max(extracted_df$BasePair, na.rm = TRUE)

# Define a 200kb window around the lowest and highest base pairs
lower_bound <- min_basepair - windows
upper_bound <- max_basepair + windows

# Define the output path in the folder /project2/mstephens/fsusie_gtex/
output_vcf <- paste0( "/project2/mstephens/fsusie_gtex/temp/extracted_snps_",i,".vcf")


datadir <- "/project2/mstephens/gtex"

plink_exec <- file.path(datadir, "plink2")
geno_prefix <- file.path(
  datadir,
  "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv"
)

# Check the actual inputs before launching PLINK.
required <- c(plink_exec, paste0(geno_prefix, c(".bed", ".bim", ".fam")))
missing <- required[file.access(required, 4) != 0L]
if (length(missing)) {
  stop("Missing or unreadable files:\n", paste(missing, collapse = "\n"))
}

outdir <- "/project2/mstephens/fsusie_gtex/temp"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Use the gene ID, rather than i, which your later loops overwrite.
output_prefix <- tempfile(
  pattern = paste0(ensembl_id, "_"),
  tmpdir = outdir
)
vcf_file <- paste0(output_prefix, ".vcf")

status <- system2(
  plink_exec,
  args = c(
    "--bfile", shQuote(geno_prefix),
    "--chr", sub("^chr", "", extracted_df$Chromosome[2]),
    "--from-bp", format(max(1, lower_bound), scientific = FALSE),
    "--to-bp", format(upper_bound, scientific = FALSE),
    "--snps-only",
    "--max-alleles", "2",
    "--rm-dup", "exclude-all",
    "--threads", "2",
    "--memory", "8000",
    "--export", "vcf-4.2", "id-paste=iid",
    "--out", shQuote(output_prefix)
  )
)

if (status != 0L || !file.exists(vcf_file)) {
  stop("Genotype extraction failed; inspect the PLINK output.")
}

vcf_data <- vcfR::read.vcfR(vcf_file)




vcf_meta <- getFIX(vcf_data)  # Extract fixed data (like CHROM, POS, etc.)
vcf_genotypes <- extract.gt(vcf_data)  # Extract genotypes


info_SNP =  data.frame(vcf_meta )


rm(vcf_data)
rm(vcf_meta)


# Function to transform genotypes
transform_genotype <- function(genotype) {
  if (is.na(genotype)) return(NA_real_)
  genotype <- chartr("|", "/", genotype)

  if (genotype == "0/0") {
    return(0)
  } else if (genotype == "0/1" || genotype == "1/0") {
    return(1)
  } else if (genotype == "1/1") {
    return(2)
  } else {
    return(NA)  # Handle missing or unexpected values
  }
}


# Assuming your matrix is called 'matrix_data'
vcf_genotypes[is.na(vcf_genotypes)] <- -9

# Apply the transformation across the whole dataframe
transformed_genotypes <- apply(as.matrix(vcf_genotypes) , 2, function(col) {
  sapply(col, transform_genotype)
})

# Convert the transformed data back into a data frame for easy handling
X <- as.data.frame(t(transformed_genotypes))
dim(X)
keep_snp <- colSums(is.na(X)) == 0
X <- X[, keep_snp, drop = FALSE]
info_SNP <- info_SNP[keep_snp, , drop = FALSE]
if (!ncol(X)) stop("No SNPs remain after removing missing genotypes.")



  ind_names=  tt[,1]
  data=tt[,-1]
  data=data.frame(data)

  # Assuming ind_names and data are already correctly defined
  base_pair <- extracted_df$BasePair[-1]

  # Calculate bin size
  total_bp <- ncol(data)
  if ( total_bp<n_bins){

    n_bins=  2^ max( which (2^(1:9)< total_bp))
    desired_bins=n_bins
  }
  bin_size <- total_bp / desired_bins

  base_pair =extracted_df$BasePair[-1]
  # Initialize the binned data frame
  binned_data <- (matrix(ncol = desired_bins, nrow = nrow(data)))
  start_bin <- rep(NA, desired_bins)

  for (i in 1:desired_bins) {
    # Calculate the start and end of the current bin
    start_col <- which.min(abs(base_pair - (min(base_pair) + (i - 1) * bin_size)))
    end_col <- which.min(abs(base_pair - (min(base_pair) + i * bin_size)))

    if (length(start_col:end_col) <floor( bin_size)  ) {
      start_bin[i] <- start_bin[i-1]+floor(bin_size)
      if(length(start_col:end_col)==1){
        binned_data[, i] =  data[, start_col]
      }else{
        binned_data[, i] = rowSums(data[, start_col:end_col], na.rm = TRUE)
      }

    } else {



      start_bin[i] <- base_pair[start_col]
      binned_data[, i] = rowSums(data[, start_col:end_col], na.rm = TRUE)
    }

  }

  Y=   binned_data

  pos= start_bin
  start_bin
  # View the binned data with positions

  # Match all tissue samples to donor genotypes (one X row per RNA sample).
  sample_ids <- as.character(tt[[1L]])
  info_ind_gen <- data.frame(
    ind = sub("^([^-]+-[^-]+).*", "\\1", sample_ids),
    tissue = sub(".*\\.", "", sample_ids),
    full_name = sample_ids,
    stringsAsFactors = FALSE
  )

  genotype_id <- sub("_.*", "", rownames(X))
  stopifnot(!anyDuplicated(genotype_id))
  rows_y <- which(info_ind_gen$ind %in% genotype_id)
  if (!length(rows_y)) stop("No coverage samples match the genotype donors.")

  sub_ind <- info_ind_gen[rows_y, , drop = FALSE]
  sub_ind$extracted <- chartr("-", ".", sub("\\..*", "", sub_ind$full_name))
  rows_x <- match(sub_ind$ind, genotype_id)

  # Calculate MAF once per matched donor, before repeating rows across tissues.
  af <- colMeans(as.matrix(X[unique(rows_x), , drop = FALSE])) / 2
  maf_filter <- pmin(af, 1 - af)
  keep_snp <- is.finite(maf_filter) & maf_filter >= 0.05 & maf_filter > 0
  if (!any(keep_snp)) stop("No SNPs remain after MAF filtering.")
  info_SNP <- info_SNP[keep_snp, , drop = FALSE]
  X <- as.matrix(X[rows_x, keep_snp, drop = FALSE])
  Y <- binned_data[rows_y, , drop = FALSE]
  rownames(X) <- rownames(Y) <- sub_ind$full_name
  stopifnot(identical(rownames(X), rownames(Y)))

  # Match library sizes by RNA sample ID, in exactly the same order as X and Y.
  load("/project2/mstephens/wdenault/GTEX_analysis_Fsusie/sum_count.RData")
  if (is.null(dim(sum_count))) {
    library_ids <- names(sum_count)
    library_values <- sum_count
  } else {
    if (ncol(sum_count) != 1L) stop("Select the library-total column of sum_count.")
    library_ids <- rownames(sum_count)
    library_values <- sum_count[, 1L]
  }
  stopifnot(!is.null(library_ids), is.numeric(library_values))
  library_ids <- chartr("-", ".", library_ids)
  stopifnot(!anyDuplicated(library_ids))
  library_rows <- match(sub_ind$extracted, library_ids)
  if (anyNA(library_rows)) {
    stop("Missing library sizes for: ",
         paste(sub_ind$full_name[is.na(library_rows)], collapse = ", "))
  }
  library_size <- as.numeric(library_values[library_rows])
  stopifnot(all(is.finite(library_size)), all(library_size > 0))
  size_factor <- library_size / mean(library_size)
  Y_cor <- sweep(Y, 1L, size_factor, "/")
  X_cor <- X

  res <- susiF(Y=   log(Y_cor+1)  ,
               X= X_cor   ,
               L_start=5 ,
               L=20,
               pos=pos,
               cal_obj =TRUE,
               nullweight=1,
               verbose=TRUE,
               maxit=40,
               max_SNP_EM = 100,
               cor_small = TRUE   )


  facto= flash(log(Y_cor+1),ebnm_fn = ebnm_point_exponential )
  Y_EBMF=facto$L_pm[, 1: min(20, ncol(  facto$L_pm )), drop=FALSE]
  library(mvsusieR)
X0 <- X_cor  # MAF filtering already removed constant SNPs.
  prior <- create_mixture_prior(R = ncol(Y_EBMF))
  mv_res= mvsusie(Y=Y_EBMF, X=X0,L=10, prior_variance = prior)


  mv_res$sets

  if (length(res$fitted_func)) {
    plot(exp(res$fitted_func[[1]]-1)*apply(log(Y_cor+1), 2, mean))
  }

  out <- list(res=res,
              size_factor_global=size_factor,
              size_factor_local=size_factor,
              Y=Y_cor   ,
              X=X,
              mv_res=mv_res,
              info_SNP=info_SNP ,
              info_ind_gen=sub_ind,
              facto=facto,
              locus=c( min(base_pair), max(base_pair)),
              pos= ( start_bin+0.5*bin_size),
              start_bin = start_bin,
              end_bin = ( start_bin+  bin_size),
              chr=extracted_df$Chromosome[2])
  output_file <- file.path(
    "/project2/mstephens/wdenault/ERC_case_study",
    paste0(ensembl_id, "_all_tissues.RData")
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  save(out, file = output_file)
  message("Saved: ", output_file)

