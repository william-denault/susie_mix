rm(list = ls())
dir_name= "/project2/mstephens/wdenault/susie_mix/results/first_test/"
lf= list.files(dir_name)


out_l= list()
for ( k in 1: length(lf)){

  load(paste0(dir_name, lf[k]))

  out_l[[k]]  =c( length(out$res_add$sets$cs),
                   length(out$res_mix$sets$cs))



}



res_df= do.call(rbind, out_l)
res_df
table(res_df[,1],res_df[,2])
k=159

 lf[k]


 ENSG00000154146
