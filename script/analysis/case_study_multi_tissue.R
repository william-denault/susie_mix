

print( res_idx[ order(res_idx$ncs_susie_mix, decreasing = TRUE)[1:100],1:30] )



table( res_idx$ncs_susie_mix,res_idx$ncs_susie )
print( res_idx[ order(res_idx$ncs_susie_mix, decreasing = TRUE)[1:100],1:30] )




lol = res_idx[ which(res_idx$ncs_susie_mix==2 &res_idx$ncs_susie ==1),1:20]


lol[order(lol$gene),]
gene_name= "HLA-DQA2"#"ABO"#"CCZ1"
tissue="Adipose Tissue"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



names(out)
k=which(names(out)=="Adipose Tissue")
par(mfrow=c(3,2))


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



library(susieR)




gene_name= "CCZ1"#"ABO"#
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))



par(mfrow=c(4,2))

for ( k in 1: length( (out))){


  susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
  susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


  abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
 abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
}


### Reaaly interesting case  -----

gene_name= "SCAMP5"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))
names(out)
k=9
par(mfrow=c(3,2))

k=4
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=6
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=8
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))

k=9
abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=12
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)
k=13

susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=14
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=+out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=16
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


### Not able to select a type of effect -----



gene_name= "NPIPA1"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))


par(mfrow=c(3,2))

k=which(names(out) =="Heart" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=which(names(out) =="Nerve" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=which(names(out) =="Skin" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


par(mfrow=c(1,1))






gene_name= "MUC20"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))


par(mfrow=c(4,2))

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

k=which(names(out) =="Muscle" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

k=which(names(out) =="Thyroid" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

par(mfrow=c(1,1))






res_idx[ order(res_idx$ncs_susie_mix, decreasing = TRUE),]



gene_name= "HPR"
out<- readRDS(paste0("/project2/mstephens/wdenault/susie_mix/results/",gene_name,".rds"))


par(mfrow=c(2,2))

k=which(names(out) =="Blood Vessel" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)


k=which(names(out) =="Adipose Tissue" )
susie_plot(out[[k]]$susie_add, y="PIP", main=paste( gene_name, names(out)[k], "addivite"))
susie_plot(out[[k]]$susie_mix, y="PIP", main=paste( gene_name, names(out)[k], "mix"))


abline(v=out[[k]]$n_SNP+1, col="red", lty=2)
abline(v=-out[[k]]$n_rec_rm+2*out[[k]]$n_SNP+1, col="red", lty=2)

par(mfrow=c(1,1))



