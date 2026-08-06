
lst_1c_res=c()


for ( k in 1:nrow(res_1cs)){

  # Nice case of susie Mix pointing to dominant coding of an additive version
  out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",
                       res_1cs$gene[k],
                       ".rds"))

  o = which(  names(out) == res_1cs$tissue[k])
  if (  res_1cs$n_rec[k] ==1 ){
    over=  intersect(1*out[[o]]$n_SNP+out[[o]]$susie_add$sets$cs$L1,out[[o]]$susie_mix$sets$cs$L1)
    lst_1c_res= c(lst_1c_res, ifelse(length(over>0), "r","n_rec"))
  }
  if (  res_1cs$n_dom[k] ==1 ){
    over= intersect(2*out[[o]]$n_SNP+out[[o]]$susie_add$sets$cs$L1,out[[o]]$susie_mix$sets$cs$L1)
    lst_1c_res= c(lst_1c_res, ifelse(length(over>0), "d","n_dom"))
  }

print(k)
}

table( lst_1c_res)

names(out)
k=27
par(mfrow=c(1,2))
susie_plot(out[[k]]$susie_add, y="PIP", main="Vagina")
susie_plot(out[[k]]$susie_mix, y="PIP")

abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))

out[[k]]$n_SNP

out[[k]]$susie_add$sets

out[[k]]$susie_mix$sets



2*out[[k]]$n_SNP+out[[k]]$susie_add$sets$cs$L1

