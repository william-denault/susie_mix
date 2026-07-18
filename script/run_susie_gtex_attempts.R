rm(list=ls())
lf= list.files("/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/")

i = 1
library(data.table)
library(vcfR)
library(fsusieR)
tt = fread(paste0("/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/", lf[i]),
           sep = ",",header = TRUE)


info_ind_gen= data.frame(ind=rep(NA, nrow(tt)),
                         tissue=  rep(NA, nrow(tt)),
                         name= rep(NA, nrow(tt)),
                         full_name= rep(NA, nrow(tt)))
for ( te in 1: nrow(tt)){

  ind_name <-  tt[te,1]
  # Extract everything before the second hyphen using a regular expression
  before_second_hyphen <- sub("^([^\\-]+\\-[^\\-]+).*", "\\1", ind_name)

  # Extract everything after the period
  after_period <- sub(".*\\.", "", ind_name)

  info_ind_gen$ind[te]= before_second_hyphen
  info_ind_gen $tissue[te]=after_period
  info_ind_gen$name[te]=  sub("\\..*", "", ind_name)
  info_ind_gen $ full_name[te] =  unlist(   ind_name)

}


#count= fread("/project2/mstephens/fsusie_gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
count= fread("/project2/mstephens/fsusie_gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

#which(colnames( count)%in% info_ind_gen$name[which(info_ind_gen$tissue=="Brain_Cortex")])
#which(colnames(raw_count)%in% info_ind_gen$name[which(info_ind_gen$tissue=="Brain_Cortex")])
raw_count= (count[,-c(1,2)])

raw_count= as.matrix(raw_count)

brain_count = raw_count [ ,
                          which(colnames(raw_count)%in% info_ind_gen$name[which(info_ind_gen$tissue=="Brain_Cortex")])
]




gen_count_brain =rowSums(brain_count)
#names(gen_count_brain)= count[,1]
hist(log10(gen_count_brain))

gene_names=count[,1]
gene_names_clean = rep(NA, nrow(count))
for ( te in 1: nrow(count )){
  gene_names_clean[te ] =sub("\\..*", "", as.character(gene_names[te]))
}
names(gen_count_brain)= gene_names_clean

lf_clean <- sub("\\..*", "", lf)

subset_gene = which (gene_names_clean %in% lf_clean)
gen_count_brain_sub=  gen_count_brain[subset_gene ]


quantile(gen_count_brain_sub)
ensembl_id =names(gen_count_brain_sub)[order(gen_count_brain_sub, decreasing = TRUE)[1:1000]]




brain_count

for( i in 1:length(ensembl_id) )
{


  tt = fread(paste0("/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/ENSG00000000003.csv.gz"),
             sep = ",",header = TRUE)

  info_ind_gen= data.frame(ind=rep(NA, nrow(tt)),
                           tissue=  rep(NA, nrow(tt)),
                           full_name= rep(NA, nrow(tt)))

  for ( o in 1: nrow(tt)){

    ind_name <-  tt[o,1]
    # Extract everything before the second hyphen using a regular expression
    before_second_hyphen <- sub("^([^\\-]+\\-[^\\-]+).*", "\\1", ind_name)

    # Extract everything after the period
    after_period <- sub(".*\\.", "", ind_name)

    info_ind_gen$ind[o]= before_second_hyphen
    info_ind_gen $tissue[o]=after_period
    info_ind_gen $ full_name[o] =  unlist(   ind_name)

  }

  idx=which(info_ind_gen$tissue=="Brain_Cortex" )

  sub_ind= info_ind_gen[idx,]
  sub_ind=data.frame(sub_ind)
  sub_ind$extracted <- sub("\\..*", "", sub_ind$full_name)
  sub_ind$extracted <- gsub("-", ".", sub_ind$extracted)





  lf= list.files("/project2/mstephens/cfbuenabadn/gtex-stm/code/coverage/counts_filtered/")

  id_gene= grep(ensembl_id[i],lf)
  id_gene


  pheno=  data.frame( y = brain_count[which(gene_names_clean == ensembl_id[i]),],
                      names= colnames(brain_count))

  output_vcf <- paste0( "/project2/mstephens/fsusie_gtex/temp/extracted_snps_",id_gene,".vcf")
  if(file.exists(paste0(output_vcf, ".vcf"))){




    vcf_data <- read.vcfR(paste0(output_vcf, ".vcf"))
    # Extract metadata (positions, chromosomes, etc.)
    vcf_meta <- getFIX(vcf_data)  # Extract fixed data (like CHROM, POS, etc.)
    vcf_genotypes <- extract.gt(vcf_data)  # Extract genotypes


    info_SNP =  data.frame(vcf_meta )


    rm(vcf_data)
    rm(vcf_meta)


    # Function to transform genotypes
    transform_genotype <- function(genotype) {

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
    X<- X[, colSums(is.na(X)) == 0]

    if( ncol(X)>0){

      ind_names= pheno[,2]




      temp_name = rep(NA, nrow(X))
      for (  te in 1:nrow(X)){

        temp_name[te] =   sub("_.*", "", rownames(X)[te])


      }
      # Modify pheno$names to keep only the first 4 characters after "GTEX-"
      pheno$names2 <- sub("^([^\\-]+\\-[^\\-]+).*", "\\1", pheno$names)


      susbset_ind= temp_name[which(temp_name %in%   pheno$names2)]


      X= X[which(temp_name %in% susbset_ind), ]

      temp2= temp_name[which(temp_name %in% susbset_ind)]
      X=X[order(temp2),]

      pheno=  pheno[ which(pheno$names2 %in% temp2),]

      pheno= pheno[order(pheno$names2),]
      maf_filter <-  apply(X , 2, sum)/(2*nrow(X))
      if (any(maf_filter < 0.05)) {
        info_SNP=info_SNP[which(maf_filter >0.05),]

        X <- X[, maf_filter >=  0.05]
      }


      # adjusting for the library size
      load("/project2/mstephens/wdenault/GTEX_analysis_Fsusie/sum_count.RData")
      row_names <- row.names(sum_count)
      modify_row_names <- function(x) {
        # Split each string by "."
        parts <- strsplit(x, "\\.")[[1]]
        # Keep only the first two parts, and replace the first dot with a hyphen
        modified_name <- paste(parts[1], parts[2], sep = "-")
        return(modified_name)
      }

      # Apply the function to all row names
      modified_names <- sapply(row_names, modify_row_names)

      # Print the modified row names
       # Apply the function to all row names
     df_sum_count = data.frame(count= sum_count$sum_count,
                               full_name= row.names(sum_count),
                               modified_names= modified_names)


     df_sum_count =df_sum_count [which(df_sum_count$full_name %in% sub_ind$extracted  ),]

     df_sum_count = df_sum_count[which(df_sum_count$modified_names %in% pheno$names2 ), ]
     df_sum_count = df_sum_count[order(df_sum_count$modified_names),]

     library_size= matrix(df_sum_count$count , ncol=1)
      size_factor= library_size/(mean(library_size))
      ####acounting for the library size

      y_cor =   (  pheno$y /as.vector(size_factor))
      X_cor = as.matrix(X)


      if( var(log1p(y_cor))>1e-4){
        res <- susieR::susie(X=X_cor,y=log1p(y_cor),L=20)
        out <- list(res=res,
                    y=pheno$y,
                    size_factor=size_factor ,
                    y_cor=y_cor  ,
                    X=X,
                    info_SNP=info_SNP ,
                    chr=info_SNP$CHROM[1])
        #save(out, file= paste("/project2/mstephens/fsusie_gtex/wave_results/",
        #                      ensembl_id,
        #                      ".RData", sep=""))
        save(out, file= paste("/project2/mstephens/fsusie_gtex/results_susie/",
                              ensembl_id[i],
                              ".RData", sep=""))

      }

  }

}

}
