
path_res="/project2/mstephens/wdenault/susie_mix/results/"
lf=list.files(path_res)
source("/project2/mstephens/wdenault/susie_mix/script/scan_tissue_attempt/help.R", echo=FALSE)
res_l=list()
for ( k in 1:length(lf)){
  out= readRDS(paste0(path_res, lf[k]))
  if( length(out)>2){
    print(k)
    dif_elbo=unlist(lapply( 1 : length(out), function(k)
      -max(out[[k]]$susie_add$elbo)+max(out[[k]]$susie_mix$elbo)))

  cs_s=do.call( rbind,lapply( 1 : length(out), function(k) c( length(out[[k]]$susie_add$sets$cs) ,
                                                   length(out[[k]]$susie_mix$sets$cs))))

  overlap=do.call(c,lapply( 1 : length(out),
                          function (k)  sum(susie_cs_overlap(out[[k]]$susie_add,
                                                             out[[k]]$susie_mix)$overlap_matrix)))


    res_l[[k]]= data.frame(dif_elbo = dif_elbo,
                           tissue   = names(out),
                           gene    = rep(sub("\\.rds$", "", lf[k]), length(out)),
                           ncs_susie=cs_s[,1],
                           ncs_susie_mix=cs_s[,2],
                           overlap=overlap)



  }


}


 res_summary= do.call( rbind, res_l)
 hist(res_summary$dif_elbo, nclass=1000,
      main="Difference in ELBO between fit with additive coding vs \n
      additive, recessive and dominant coding")
abline(v=0, col="red", lty=2)

idx= which(res_summary$ncs_susie>0)
hist(res_summary$dif_elbo[idx], nclass=100,
     main="Difference in ELBO between fit with additive coding vs \n
      additive, recessive and dominant coding")
abline(v=0, col="red", lty=2)




sum(res_summary$ncs_susie[idx])
sum(res_summary$ncs_susie_mix[idx])
table( res_summary$ncs_susie, res_summary$ncs_susie_mix)
table( res_summary$ncs_susie[idx], res_summary$ncs_susie_mix[idx])

res_summary[order(res_summary$dif_elbo, decreasing=TRUE)[1:10],]
 out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/APOA4.rds")
k=1
par(mfrow=c(1,2))
plot(out[[k]]$susie_add$pip)
plot(out[[k]]$susie_mix$pip)

out[[k]]$susie_mix$sets$purity
out[[k]]$susie_add$sets$purity
max(out[[k]]$susie_add$elbo)-max(out[[k]]$susie_mix$elbo)
k=k+1

out= ANKRD31
par(mfrow=c(1,1))
dif_elbo=unlist(lapply( 1 : length(out), function(k)
  max(out[[k]]$susie_add$elbo)-max(out[[k]]$susie_mix$elbo)))
plot(dif_elbo)
order(dif_elbo)




