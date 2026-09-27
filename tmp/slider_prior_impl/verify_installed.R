task_lib <- "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/slider_prior_impl/library"
library(susieSlide, lib.loc=task_lib)
set.seed(27)
X <- matrix(rbinom(300*12,2,.4),300,12)
y <- X[,3]-.5*(X[,3]==1)+rnorm(300,sd=.6)
fit <- suppressMessages(susie(X,y,L=2))
stopifnot(fit$converged, identical(fit$delta_grid,seq(-1,1,length.out=17)),
  isTRUE(all.equal(fit$delta_prior,rep(1/17,17))),
  identical(dim(fit$alpha_delta),c(2L,12L,17L)),
  isTRUE(all.equal(predict(fit,newx=X),fit$fitted,tolerance=1e-10)))
for(value in c(.5,NA_real_,Inf,3)) {
  invalid <- X; invalid[1,1] <- value
  rejected <- tryCatch({susie(invalid,y); FALSE},error=function(e) TRUE)
  stopifnot(rejected)
}
w <- rep(0,17); w[c(1,9,17)] <- c(.2,.5,.3)
weighted <- suppressMessages(susie(X,y,L=2,delta_prior=w,model_init=fit))
stopifnot(isTRUE(all.equal(weighted$delta_prior,w)),
          all(weighted$alpha_delta[,,-c(1,9,17)]==0))
cat("Installed package:",as.character(packageVersion("susieSlide")),"\n")
cat("Fixed-prior fit, predictions, warm start and strict 0/1/2 input validation: PASS\n")
