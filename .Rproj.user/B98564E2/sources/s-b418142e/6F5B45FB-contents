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
for ( i in 1: nrow(tt)){

  ind_name <-  tt[i,1]
  # Extract everything before the second hyphen using a regular expression
  before_second_hyphen <- sub("^([^\\-]+\\-[^\\-]+).*", "\\1", ind_name)

  # Extract everything after the period
  after_period <- sub(".*\\.", "", ind_name)

  info_ind_gen$ind[i]= before_second_hyphen
  info_ind_gen $tissue[i]=after_period
  info_ind_gen$name[i]=  sub("\\..*", "", ind_name)
  info_ind_gen $ full_name[i] =  unlist(   ind_name)

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


### check overlapp  between preselected gene by Carlos


gen_count_brain =rowSums(brain_count)
#names(gen_count_brain)= count[,1]
hist(log10(gen_count_brain))

gene_names=count[,1]
gene_names_clean = rep(NA, nrow(count))
for ( i in 1: nrow(count )){
  gene_names_clean[i] =sub("\\..*", "", as.character(gene_names[i]))
}
names(gen_count_brain)= gene_names_clean

lf_clean <- sub("\\..*", "", lf)

subset_gene = which (gene_names_clean %in% lf_clean)
gen_count_brain_sub=  gen_count_brain[subset_gene ]

hist(log10(gen_count_brain_sub ))
par(mfrow=c(1,2))
plot( (gen_count_brain_sub[order(gen_count_brain_sub, decreasing = TRUE)]),  ylab = " tpm")
plot(log10(gen_count_brain_sub[order(gen_count_brain_sub, decreasing = TRUE)]),  ylab = "log10 tpm")
par(mfrow=c(1,1))

quantile(gen_count_brain_sub)
ensembl_id =names(gen_count_brain_sub)[order(gen_count_brain_sub, decreasing = TRUE) ]




for( i in 1:length(ensembl_id) )
{


  my_line <- paste('source("/project2/mstephens/wdenault/susie_mix/script/analysis/analysis_region.R");  ',
                   'analyse_region(ensembl_id=\"' , ensembl_id[i],'\", windows=', 200000,')', sep="")
  write(my_line,
        file=paste("/project2/mstephens/wdenault/susie_mix/script/scan_job/job_fsusie_",i,".R", sep=""),
        append=TRUE)

}


