source("/project2/mstephens/wdenault/susie_mix/script/analysis/descriptive results.R", echo=TRUE)


table( res_idx$ncs_susie, res_idx$ncs_susie_mix)



res_idx[ which(res_idx$ncs_susie==1 & res_idx$ncs_susie_mix==3), 1:20]

tt= res_idx[ which(res_idx$ncs_susie< res_idx$ncs_susie_mix), 1:20]




tt= res_idx[ which(res_idx$ncs_susie==1 &res_idx$ncs_susie_mix==2), 1:20]
## Nice secondary dominat secondary signal ----
l=1
gene_name= tt$gene[l]
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=which(names(out)==tt$tissue[l])
par(mfrow=c(1,2))

susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))


l=l+1
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


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))
out[[k]]$mean_read

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets

2*out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L2
