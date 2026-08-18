
load("/project2/mstephens/wdenault/susie_mix/res_summary.RData")
 sum(res_summary$perm_cs_susie)
 sum(res_summary$perm_cs_susie_mix)

 sum(res_summary$ncs_susie)
 sum(res_summary$ncs_susie_mix)

 idx= which(res_summary$min_pv< 5e-8)



 res_idx= res_summary[idx,]
 plot(-log10(res_idx$min_pv))
 abline(col="red", h= -log10(5e-8))
 sum(res_summary$perm_cs_susie[idx])
 sum(res_summary$perm_cs_susie_mix[idx])

 sum(res_summary$ncs_susie[idx])
 sum(res_summary$ncs_susie_mix[idx])




 hist(res_summary$dif_elbo, nclass=1000,
      main="Difference in ELBO between fit with additive coding vs \n
      additive, recessive and dominant coding")
abline(v=0, col="red", lty=2)
hist(-2 *(res_summary$log_lik_add- res_summary$log_lik_mix), nclass=1000,
     main="-2 (log lik additive coding - log lik mix coding")
abline(v=0, col="red", lty=2)

hist(-2 *(res_summary$log_lik_add- res_summary$log_lik_mix), nclass=300,
     xlim=c(-100, 300),
     xlab="-2 (log lik additive coding - log lik mix coding)",
     main="-2 (log lik additive coding - log lik mix coding)")
abline(v=0, col="red", lty=2)



hist(-2 *(res_summary$log_lik_add_perm- res_summary$log_lik_mix_perm), nclass=1000,
     xlim=c(-100, 300),
     main="-2 (log lik diff with additive coding vs \n
      additive, recessive and dominant coding")
abline(v=0, col="red", lty=2)


hist(res_summary$dif_elbo_perm, nclass=1000,
     main="Difference in ELBO perm")
abline(v=0, col="red", lty=2)


idx= which(res_summary$ncs_susie>0)
hist(res_summary$dif_elbo[idx], nclass=100,
     main="Difference in ELBO between fit with additive coding vs \n
      additive, recessive and dominant coding")
abline(v=0, col="red", lty=2)

idx= which(res_summary$ncs_susie>0)
hist(-2 *(res_summary$log_lik_add_perm- res_summary$log_lik_mix_perm)[idx], nclass=300,
     xlab="-2 (log lik additive coding - log lik mix coding)",
     main="-2 (log lik additive coding - log lik mix coding) with more than 0 cs")
abline(v=0, col="red", lty=2)

res_idx= res_idx[which(res_idx$mean_count>100),]




table(res_idx$n_add)
table(res_idx$n_rec)
table(res_idx$n_dom)

table(res_idx$n_add_perm)
table(res_idx$n_rec_perm)

table(res_idx$n_dom_perm)




image( as.matrix(res_idx[,c("n_add", "n_rec" ,"n_dom") ]))
table( res_summary$ncs_susie, res_summary$ncs_susie_mix)
table( res_summary$ncs_susie[idx], res_summary$ncs_susie_mix[idx])

res_idx[order(res_idx$dif_elbo, decreasing=TRUE)[1:50],]

res_idx[order(res_idx$min_pv, decreasing=FALSE)[1:50],]


sum(ifelse(res_idx[,c("n_add")] >0,1,0) * ifelse(res_idx[,c( "n_dom" )]  >0,1,0))

res_idx[which(ifelse(res_idx[,c("n_add")] >0,1,0) * ifelse(res_idx[,c( "n_dom" )]  >0,1,0)==1),]

sum(ifelse(res_idx[,c("n_rec")] >0,1,0) *ifelse(res_idx[,c( "n_dom" )] >0,1,0))

sum(ifelse(res_idx[,c("n_add")]>0,1,0) *ifelse(res_idx[,c("n_rec")] >0,1,0)*ifelse(res_idx[,c( "n_dom" )]>0,1,0))



plot(res_idx$overlap)
table(res_idx$ncs_susie,res_idx$ncs_susie_mix)
which(!(res_idx$ncs_susie ==res_idx$ncs_susie_mix))
#proportion of number of different cs
length(which(!(res_idx$ncs_susie ==res_idx$ncs_susie_mix)) )/ nrow(res_idx)


 res_idx[which(!(res_idx$ncs_susie ==res_idx$ncs_susie_mix)),]


 table(res_idx$tissue[which(!(res_idx$ncs_susie ==res_idx$ncs_susie_mix)) ])


# percentage in which SUSiE and SuSiE mix output 1 CS but disagree on the "variant" coding
table( res_idx$ overlap[which(res_idx$ncs_susie==1 & res_idx$ncs_susie_mix==1  )  ])[1]/ length(which(res_idx$ncs_susie==1 & res_idx$ncs_susie_mix==1  ))



### region in which suise and susie mix give a single CS that are different -----

 res_1cs=  res_idx[which(res_idx$ncs_susie==1 & res_idx$ncs_susie_mix==1 & res_idx$overlap==0),]
res_1cs[ which(res_1cs$n_dom== 1),]

res_1cs$t_log_d=-2*(res_1cs$log_lik_add- res_1cs$log_lik_mix)

table(res_1cs$gene)



res_idx[ order(res_idx$ncs_susie, decreasing = TRUE),]



### example of 1 cs susie and susie mix but different coding choice -----

table(res_1cs$tissue)

# Nice case of susie Mix pointing to dominant coding of an additive version


out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/AATK.rds")
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




# Nice case of susie Mix pointing to dominant coding of an additive version
out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/ABCA3.rds")
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








### example of 0 cs susie mix -----

out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/ABCA3.rds")
names(out)
k=3
par(mfrow=c(1,2))
susie_plot(out[[k]]$susie_add, y="PIP", main="Blood")
susie_plot(out[[k]]$susie_mix, y="PIP")

abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=2*out[[k]]$n_SNP+1, col="red", lty=2)
par(mfrow=c(1,1))



out[[k]]$mean_read




out<- readRDS("/project2/mstephens/wdenault/susie_mix/results/ABCC6.rds")
names(out)
k=9
par(mfrow=c(1,2))
susie_plot(out[[k]]$susie_add, y="PIP", main="Heart")
susie_plot(out[[k]]$susie_mix, y="PIP")
par(mfrow=c(1,1))



out[[k]]$mean_read
out[[k]]$susie_add$sets
out[[k]]$susie_mix$sets

