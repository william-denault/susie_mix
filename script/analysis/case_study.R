

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
