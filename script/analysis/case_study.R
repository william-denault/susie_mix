####Simple cases study
library(susieR)
## Nice secondary recessive secondary signal that is additive ----
out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/AGAP4.rds")
names(out)
k=7
par(mfrow=c(1,2))

susie_plot(out[[k]]$susie_add, y="PIP", main="Colon")
susie_plot(out[[k]]$susie_mix, y="PIP")


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)

abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

### ABO recoding effect more concentrated posterior-----
gene_name= "ABO"#"ABO"#"CCZ1"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



names(out)
k=20
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v= out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$susie_add$sets
out[[k]]$susie_mix$sets

offset= -out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1

out[[k]]$susie_mix$sets$cs$L1-offset




### Case where SuSiE cannot pick between variant  -----

gene_name= "SCAMP5"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=9
par(mfrow=c(2,2))

k=9
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=20
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets



out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1

k=6
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

## Nice secondary dominat secondary signal ----

gene_name= "GTF2H2"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=6
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets










# Case where SuSiE mix and SuSiE mix agree but SuSiE mix recode variant ----


gene_name= "ZNF232"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=14
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets

out[[k]]$n_rec_SNP+out[[k]]$n_SNP  +out[[k]]$susie_add$sets$cs$L2





## interesting  case ----


gene_name= "GSTM1"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))


par(mfrow=c(4,2))
for ( k in which(names(out) %in% res_idx$tissue[ which (res_idx$gene==gene_name) ])){




  susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
  susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


  abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
 abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

}

k=1 #14 20 25
par(mfrow=c(4,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)



k=14

susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=20

susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=25

susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets

2*out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L2





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



out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1




### Interesting mutli tissue hetoreogneitu
out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/ZP3.rds")



names(out)
k=1
par(mfrow=c(3,2))

susie_plot(out[[k]]$susie_add, y="PIP", main="Adipose Tissue")
susie_plot(out[[k]]$susie_mix, y="PIP")



abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=14


susie_plot(out[[k]]$susie_add, y="PIP", main="Nerve")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=20


susie_plot(out[[k]]$susie_add, y="PIP", main="Skin")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)





### Genes in which susie and susie donot agree over multiple tissue -----






out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/CCZ1.rds")
names(out)
k=1
par(mfrow=c(3,2))

susie_plot(out[[k]]$susie_add, y="PIP", main="Adipose Tissue")
susie_plot(out[[k]]$susie_mix, y="PIP")



abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=4


susie_plot(out[[k]]$susie_add, y="PIP", main="Blood Vessel")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=8


susie_plot(out[[k]]$susie_add, y="PIP", main="Esophagus")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
k=14

susie_plot(out[[k]]$susie_add, y="PIP", main="Nerve")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
k=20

susie_plot(out[[k]]$susie_add, y="PIP", main="Skin")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
k=25

susie_plot(out[[k]]$susie_add, y="PIP", main="Thyroid")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

k=1
out[[k]]$susie_add$sets
out[[k]]$susie_mix$sets

out[[k]]$n_SNP

### OLD------





### OLD------
## Nice secondary recessive secondary signal ----
out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/BIRC7.rds")
names(out)
k=1
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main="Adipose Tissue")
susie_plot(out[[k]]$susie_mix, y="PIP")


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets



out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1






## Weird case in which susie mix ouput less cs but that do overlap with SuSiE ----
"APIP"

"ABO"


## Potential case of pb ----

out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/ABO.rds")
names(out)
par(mfrow=c(1,2))
k=20


susie_plot(out[[k]]$susie_add, y="PIP", main="Skin")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))

out[[k]]$n_SNP


## investigations heavy disagreement ----

out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/HPR.rds")
names(out)
par(mfrow=c(1,2))
k=4

susie_plot(out[[k]]$susie_add, y="PIP", main="Blood Vessel")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))

