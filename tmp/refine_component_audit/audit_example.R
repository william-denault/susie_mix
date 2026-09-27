
.libPaths(c("C:/Users/willi/AppData/Local/R/win-library/4.6", .libPaths()))
library(susieR)
options(width=160)
cat("susieR version:", as.character(packageVersion("susieR")), "\n")
data(N3finemapping)
X_raw <- N3finemapping$X
n <- nrow(X_raw); p <- ncol(X_raw)
recode_snp <- function(x) {
  sorted_vals <- names(sort(table(x), decreasing=TRUE))
  ad <- match(x, sorted_vals)-1
  list(ad=ad,dom=as.integer(ad>=1),rec=as.integer(ad==2))
}
AD <- DOM <- REC <- matrix(0,n,p)
for(j in seq_len(p)) { rc<-recode_snp(X_raw[,j]); AD[,j]<-rc$ad; DOM[,j]<-rc$dom; REC[,j]<-rc$rec }
cat("Entries recoded to 3 before capping:",sum(AD==3,na.rm=TRUE),"\n")
AD[AD==3] <- 2
cat("Recessive coding inconsistent with capped AD:",sum(REC!=(AD==2),na.rm=TRUE),"\n")
set.seed(1)
causal_snps <- c(12,92,128,165)
y <- rnorm(n,sd=.5)
y <- y+2*scale(AD[,12])+2*scale(DOM[,92])+2*scale(REC[,128])
X <- cbind(AD,DOM,REC)
colnames(X) <- c(paste0("AD_",seq_len(p)),paste0("DOM_",seq_len(p)),paste0("REC_",seq_len(p)))
code_mass <- function(fit) cbind(additive=rowSums(fit$alpha[,seq_len(p),drop=FALSE]),dominant=rowSums(fit$alpha[,p+seq_len(p),drop=FALSE]),recessive=rowSums(fit$alpha[,2*p+seq_len(p),drop=FALSE]))
describe <- function(fit,label) {
  cat("\n",label,"\n")
  cat("Converged",fit$converged,"iterations",fit$niter,"ELBO",tail(fit$elbo,1),"CS",length(fit$sets$cs),"\n")
  mass <- code_mass(fit)
  d <- data.frame(component=seq_len(nrow(fit$alpha)),V=rep_len(fit$V,nrow(fit$alpha)),reported_CS=seq_len(nrow(fit$alpha)) %in% fit$sets$cs_index,max_alpha=apply(fit$alpha,1,max),lbf=fit$lbf,mass)
  print(d,digits=9,row.names=FALSE)
  print(fit$sets$purity)
  keep <- d$V>0
  cat("One-step coding update using V>0 (A,D,R):",colSums(mass[keep,,drop=FALSE])/sum(keep),"\n")
  pure <- susie_get_cs(fit,X=X,min_abs_corr=.5)
  cat("Post-hoc purity>=0.5 components:",pure$cs_index,"\n")
  if(length(pure$cs_index)) cat("CS-only sensitivity (A,D,R):",colMeans(mass[pure$cs_index,,drop=FALSE]),"\n")
  raw <- susie_get_cs(fit,min_abs_corr=0,dedup=FALSE)
  cat("Undeduplicated sets above variance threshold:",names(raw$cs),"\n")
  print(data.frame(component=names(raw$cs),size=lengths(raw$cs),duplicates_previous=duplicated(raw$cs)),row.names=FALSE)
  invisible(d)
}
set.seed(123)
f0 <- susieR::susie(X,y,L=10,verbose=FALSE,min_abs_corr=0)
d0 <- describe(f0,"USER FIT: defaults, purity=0")
flush.console()
set.seed(123)
fr <- susieR::susie(X,y,L=10,verbose=FALSE,min_abs_corr=.5,refine=TRUE,max_iter=1000,tol=1e-5)
dr <- describe(fr,"REFINED: purity=.5, max_iter=1000, tol=1e-5")
flush.console()
set.seed(123)
ft <- susieR::susie(X,y,L=10,verbose=FALSE,min_abs_corr=.5,refine=FALSE,max_iter=1000,tol=1e-5)
dt <- describe(ft,"TIGHT CONVERGENCE WITHOUT REFINEMENT")
cat("Tight vs refined alpha max difference:",max(abs(ft$alpha-fr$alpha)),"\n")
saveRDS(list(user=f0,refined=fr,tight=ft,X=X,y=y,version=packageVersion("susieR")), "tmp/refine_component_audit/example_fits.rds")
write.csv(rbind(transform(d0,fit="user"),transform(dr,fit="refined"),transform(dt,fit="tight")), "tmp/refine_component_audit/components.csv",row.names=FALSE)

