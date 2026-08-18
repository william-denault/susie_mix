ad_gwas_top50 <- c(
  "APOE",
  "BIN1",
  "CR1",
  "CLU",
  "PICALM",
  "ABCA7",
  "SORL1",
  "TREM2",
  "INPP5D",
  "CD33",
  "MS4A4A",
  "MS4A6A",
  "CD2AP",
  "EPHA1",
  "PTK2B",
  "HLA-DQA1",
  "ADAM10",
  "ACE",
  "PLCG2",
  "ABI3",
  "SPI1",
  "FERMT2",
  "CASS4",
  "SLC24A4",
  "RIN3",
  "ADAMTS4",
  "APH1B",
  "ABCA1",
  "ATP8B4",
  "IL34",
  "MME",
  "SORT1",
  "NCK2",
  "WWOX",
  "TSPAN14",
  "SCIMP",
  "SHARPIN",
  "SNX1",
  "BLNK",
  "RBCK1",
  "KAT8",
  "DOC2A",
  "ECHDC3",
  "PILRA",
  "ZCWPW1",
  "CELF1",
  "IQCK",
  "GRN",
  "APP",
  "TMEM106B"
)


o=1
gene_name= ad_gwas_top50[o]
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)

par(mfrow=c(3,2))

k=4
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=5
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))

abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=14
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))

abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


print(gene_name)
o=o+1

#APOE
"PILRA"
