
path_res="/project2/mstephens/wdenault/susie_mix/results/"
lf=list.files(path_res)
length(lf)
source("/project2/mstephens/wdenault/susie_mix/script/scan_tissue_attempt/help.R", echo=FALSE)
res_l=list()
for ( k in 1:length(lf)){
  out= readRDS(paste0(path_res, lf[k]))
  if( length(out)>2){

    dif_elbo=unlist(lapply( 1 : length(out), function(k)
      -max(out[[k]]$susie_add$elbo)+max(out[[k]]$susie_mix$elbo)))
    dif_elbo_perm =unlist(lapply( 1 : length(out), function(k)
      -max(out[[k]]$susie_add_perm$elbo)+max(out[[k]]$susie_mix_perm$elbo)))
    cs_s=do.call( rbind,lapply( 1 : length(out), function(k) c( length(out[[k]]$susie_add$sets$cs) ,
                                                                length(out[[k]]$susie_mix$sets$cs))))
    min_pv=do.call( rbind,lapply( 1 : length(out), function(k)  out[[k]]$min_pv))
    overlap=do.call(c,
                    lapply( 1 : length(out),
                            function (k)  sum(susie_cs_overlap(out[[k]]$susie_add,
                                                               out[[k]]$susie_mix)$overlap_matrix)
                    )
    )

    log_lik_add=do.call( rbind,lapply( 1 : length(out), function(k) {
      out[[k]]$susie_add$elbo[length(out[[k]]$susie_add$elbo)]+sum(out[[k]]$susie_add$KL)}))

    log_lik_mix=do.call( rbind,lapply( 1 : length(out), function(k) {
      out[[k]]$susie_mix$elbo[length(out[[k]]$susie_mix$elbo)]+sum(out[[k]]$susie_mix$KL)

              }))

    log_lik_add_perm=do.call( rbind,lapply( 1 : length(out), function(k) {
      out[[k]]$susie_add_perm$elbo[length(out[[k]]$susie_add_perm$elbo)]+sum(out[[k]]$susie_add_perm$KL)
    }))


    log_lik_mix_perm=do.call( rbind,lapply( 1 : length(out), function(k) {
      out[[k]]$susie_mix_perm$elbo[length(out[[k]]$susie_mix_perm$elbo)]+sum(out[[k]]$susie_mix_perm$KL)
    }))


    type_cs=do.call(rbind,
                    lapply( 1 : length(out),
                            function(k){
                              n_add=0
                              n_rec=0
                              n_dom=0
                              if (length(out[[k]]$susie_mix$sets$cs)>0){
                                for( l in 1:length(out[[k]]$susie_mix$sets$cs)){
                                  div_num= (out[[k]]$susie_mix$sets$cs[[l]][1]-1) %/% out[[k]]$n_SNP
                                  if(div_num==0){
                                    n_add=n_add+1
                                  }
                                  if(div_num==1){
                                    n_rec=n_rec+1
                                  }
                                  if(div_num==2){
                                    n_dom=n_dom+1
                                  }
                                }
                              }
                              return(c(n_add,
                                       n_rec,
                                       n_dom)
                              )

                            }
                    )
    )

    type_cs=as.data.frame(type_cs)
    colnames(type_cs)=c("n_add",
                        "n_rec",
                        "n_dom")



    type_cs_perm=do.call(rbind,
                         lapply( 1 : length(out),
                                 function(k){
                                   n_add=0
                                   n_rec=0
                                   n_dom=0
                                   if (length(out[[k]]$susie_mix_perm$sets$cs)>0){
                                     for( l in 1:length(out[[k]]$susie_mix_perm$sets$cs)){
                                       div_num= out[[k]]$susie_mix_perm$sets$cs[[l]][1] %/% out[[k]]$n_SNP
                                       if(div_num==0){
                                         n_add=n_add+1
                                       }
                                       if(div_num==1){
                                         n_rec=n_rec+1
                                       }
                                       if(div_num==2){
                                         n_dom=n_dom+1
                                       }
                                     }
                                   }
                                   return(c(n_add,
                                            n_rec,
                                            n_dom)
                                   )

                                 }
                         )
    )

    type_cs_perm=as.data.frame(type_cs_perm)
    colnames(type_cs_perm)=c("n_add_perm",
                             "n_rec_perm",
                             "n_dom_perm")
    perm_cs_susie= do.call(c,
                           lapply( 1 : length(out),
                                   function (k)  length(out[[k]]$susie_add_perm$sets$cs)
                           )
    )
    perm_cs_susie_mix= do.call(c,
                               lapply( 1 : length(out),
                                       function (k)  length(out[[k]]$susie_mix_perm$sets$cs)
                               )
    )

    n_ind            = do.call(c,
                               lapply( 1 : length(out),
                                       function (k)   (out[[k]]$n_ind)
                               )
    )
    median_count= do.call(c,
                        lapply( 1 : length(out),
                                function (k)   (out[[k]]$median_read)
                        )
    )
    mean_count= do.call(c,
                          lapply( 1 : length(out),
                                  function (k)   (out[[k]]$mean_read)
                          )
    )
    temp_r =  data.frame(dif_elbo = dif_elbo,
                         min_pv=min_pv,
                         mean_count=   mean_count,
                         tissue   = names(out),
                         gene    = rep(sub("\\.rds$", "", lf[k]), length(out)),
                         ncs_susie=cs_s[,1],
                         ncs_susie_mix=cs_s[,2],
                         overlap=overlap,
                         perm_cs_susie=perm_cs_susie,
                         perm_cs_susie_mix=  perm_cs_susie_mix,
                         median_count=   median_count,
                         dif_elbo_perm=dif_elbo_perm,
                         n_ind=  n_ind )
    print(k)
    res_l[[k]]= cbind(temp_r, type_cs,type_cs_perm)

  }


}


 res_summary= do.call( rbind, res_l)

save(res_summary, file="res_summary.RData")
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



hist(res_summary$dif_elbo_perm, nclass=1000,
     main="Difference in ELBO perm")
abline(v=0, col="red", lty=2)


idx= which(res_summary$ncs_susie>0)
hist(res_summary$dif_elbo[idx], nclass=100,
     main="Difference in ELBO between fit with additive coding vs \n
      additive, recessive and dominant coding")
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

### region in which suise and susie mix give a single CS that are different


 res_1cs=  res_idx[which(res_idx$ncs_susie==1 & res_idx$ncs_susie_mix==1 & res_idx$overlap==0),]


res_1cs


res_idx[ order(res_idx$ncs_susie_mix, decreasing = TRUE),]
res_idx[ order(res_idx$ncs_susie_mix),]


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
