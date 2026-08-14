
# Case where SuSiE mix and SuSiE mix agree but SuSiE mix recode variant ----


Nerve
gene_name= "ZNF232"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=14
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets

2*out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L2



## leaninng torward dominant ----

gene_name= "SCAMP5"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=8
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets



out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1





### Case where SuSiE cannot pick between variant  -----

gene_name= "SCAMP5"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=9
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets



out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1


## Nice secondary dominat secondary signal ----

gene_name= "GTF2H2"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=6
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets








## Nice case where SuSiE additive is quite differnet ----

gene_name= "HPR"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=4
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets



out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1






### OLD------
## Nice secondary recessive secondary signal ----
out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/BIRC7.rds")
names(out)
k=1
par(mfrow=c(1,2))


susie_plot(out[[k]]$susie_add, y="PIP", main="Adipose Tissue")
susie_plot(out[[k]]$susie_mix, y="PIP")


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets



out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1


## Nice secondary recessive secondary signal that is additive ----
out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/AGAP4.rds")
names(out)
k=7
par(mfrow=c(1,2))

susie_plot(out[[k]]$susie_add, y="PIP", main="Colon")
susie_plot(out[[k]]$susie_mix, y="PIP")


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read


### Genes in which susie and susie donot agree over multiple tissue






out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/CCZ1.rds")
names(out)
k=1
par(mfrow=c(3,2))

susie_plot(out[[k]]$susie_add, y="PIP", main="Adipose Tissue")
susie_plot(out[[k]]$susie_mix, y="PIP")



abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)

k=4


susie_plot(out[[k]]$susie_add, y="PIP", main="Blood Vessel")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)

k=8


susie_plot(out[[k]]$susie_add, y="PIP", main="Esophagus")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
k=14

susie_plot(out[[k]]$susie_add, y="PIP", main="Nerve")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
k=20

susie_plot(out[[k]]$susie_add, y="PIP", main="Skin")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
k=25

susie_plot(out[[k]]$susie_add, y="PIP", main="Thyroid")
susie_plot(out[[k]]$susie_mix, y="PIP")
abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

k=1
out[[k]]$susie_add$sets
out[[k]]$susie_mix$sets

out[[k]]$n_SNP

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
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
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
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
