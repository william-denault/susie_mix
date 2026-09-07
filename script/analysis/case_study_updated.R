


## interesting HLA disagreement -----


gene_name= "HLA-DQA2"#"ABO"#"CCZ1"
tissue="Adipose Tissue"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



names(out)
k=which(names(out)=="Adipose Tissue")
par(mfrow=c(4,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=which(names(out)=="Blood")



susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=which(names(out)=="Colon")



susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=which(names(out)=="Nerve")


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))







#### case of most disagreement -----


gene_name= "CNTNAP3B"#"ABO"#"CCZ1"
tissue="Nerve"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



names(out)
k=which(names(out)==tissue)
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

par(mfrow=c(1,1))


gene_name= "NPIPA1"#"ABO"#"CCZ1"
tissue="Thyroid"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



names(out)
k=which(names(out)==tissue)
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

par(mfrow=c(1,1))




gene_name= "PDPR"#"ABO"#"CCZ1"
tissue="Adipose Tissue"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



names(out)
k=which(names(out)==tissue)
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

par(mfrow=c(1,1))


##### ------


gene_name= "GSTM1"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))


par(mfrow=c(1,2))
 k=3




  susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
  susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


  abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
  abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)





## Nice case where SuSiE additive is quite differnet ----

gene_name= "HPR"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=4
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets









### Interesting mutli tissue hetoreogneitu
out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/ZP3.rds")



names(out)
k=1
par(mfrow=c(2,2))

susie_plot(out[[k]]$susie_add, y="PIP", main="Adipose Tissue")
susie_plot(out[[k]]$susie_mix, y="PIP")



abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=14


susie_plot(out[[k]]$susie_add, y="PIP", main="Nerve")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)








### Genes in which susie and susie donot agree over multiple tissue -----






out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/CCZ1.rds")
names(out)
k=1
par(mfrow=c(2,2))

k=4


susie_plot(out[[k]]$susie_add, y="PIP", main="Blood Vessel")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

out[[k]]$susie_add$sets
out[[k]]$susie_mix$sets
1230+2*out[[k]]$n_SNP-out[[k]]$n_rec_rm



out[[k]]$susie_mix$sets$cs[[3]]+2*out[[k]]$n_SNP-out[[k]]$n_rec_rm
out[[k]]$susie_mix$sets$cs[[3]]


k=25

susie_plot(out[[k]]$susie_add, y="PIP", main="Thyroid")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read



1304+2*out[[k]]$n_SNP-out[[k]]$n_rec_rm



1242+2*out[[k]]$n_SNP-out[[k]]$n_rec_rm



### Interesting disagrement + secondary hbit

gene_name= "CCZ1"#"ABO"#
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



par(mfrow=c(2,2))
k=4

susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
k=24

susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)








gene_name= "MUC20"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))


par(mfrow=c(2,2))

k=which(names(out) =="Adipose Tissue" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=which(names(out) =="Blood" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

par(mfrow=c(1,1))
