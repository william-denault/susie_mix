# Real slider fits through parallel adoption, two pooled updates and warm starts.
extra <- Sys.getenv("SLIDER_PRIOR_TEST_LIB")
if (nzchar(extra)) .libPaths(c(extra, .libPaths()))
repo <- normalizePath(".", winslash = "/")
dir.create(file.path(repo, "tmp"), showWarnings = FALSE)
for (f in c("em_utils.R", "slide_prior_em_utils.R", "slide_prior_recovery.R",
            "prepare_slide_prior_em_iteration.R", "slide_prior_parallel.R"))
  source(file.path(repo, "script/scan_tissue_attempt", f))
version <- spe_check_packages()
root <- tempfile("parallel_real_",tmpdir=file.path(repo,"tmp")); dir.create(root)
scripts <- file.path(root,"script/scan_tissue_attempt"); dir.create(scripts,recursive=TRUE)
file.copy(file.path(repo,"script/scan_tissue_attempt/slide_prior_recovery.R"),scripts)
dir.create(file.path(root,"results_em/iteration_000"),recursive=TRUE)
write.csv(data.frame(iteration=0L,tissue="Brain"),file.path(root,"results_em/prior_history.csv"),row.names=FALSE)
write.csv(data.frame(chunk=1L,gene=paste0("G",1:4)),file.path(root,"results_em/iteration_000/manifest.csv"),row.names=FALSE)
set.seed(104)
data <- lapply(1:4,function(i) {
  X <- matrix(rbinom(160*6,2,.35),160,6); X[,6] <- 0; X[1:8,6] <- 1
  colnames(X) <- paste0("snp",1:6)
  list(X=X,y=X[,1]+.5*(X[,1]==1)+rnorm(160,sd=.4))
})
calls <- 0L; warm <- 0L
run <- function(target_gene,tissue_priors,previous_result=NULL,...) {
  calls <<- calls+1L; d <- data[[match(target_gene,paste0("G",1:4))]]
  a <- list(X=d$X,y=d$y,L=10,min_obs=5,standardize=FALSE,estimate_prior_method="EM",
            delta_grid=spe_grid(),delta_prior=spe_prior(tissue_priors,"Brain"),max_iter=5000,tol=1e-5)
  if (!is.null(previous_result)) { a$model_init <- previous_result$Brain$fit_slide_prior; warm <<- warm+1L }
  list(Brain=list(fit_slide_prior=do.call(susieRSlidePrior::susie,a)))
}
launch <- spe_prepare_iteration(root,"new",package_version=version)
ctx <- sre_context(root,launch$iteration_dir)
for (g in ctx$manifest$gene) {
  x <- run(g,ctx$priors); attr(x,"em_iteration") <- 0L
  saveRDS(x,sre_result(ctx,g))
}
hashes <- tools::md5sum(sre_result(ctx,ctx$manifest$gene))
q <- file.path(root,"queue.txt"); writeLines(character(),q)
for(iter in 0:2) {
  spp_init(ctx)
  p <- spp_plan(ctx,298,q)
  stopifnot(p$action == if(iter==0) "validate" else "fit")
  w <- spp_batch(ctx,p$batch)
  for(i in seq_len(nrow(w))) spp_worker_gene(ctx,w$gene[i],w$action[i],run)
  serial <- spe_pool(sre_pool_files(ctx),ctx$priors)
  stopifnot(spp_plan(ctx,298,q)$status=="complete")
  actual <- spp_require_complete(root,ctx$dir)
  stopifnot(isTRUE(all.equal(actual,serial,tolerance=1e-12)),
            identical(hashes,tools::md5sum(names(hashes))))
  if(iter<2) {
    launch <- spe_prepare_iteration(root,"new",package_version=version,completion_reader=spp_require_complete)
    ctx <- sre_context(root,launch$iteration_dir)
  }
}
stopifnot(calls==12L,warm==8L,max(abs(spe_prior(ctx$priors,"Brain")-1/17))>1e-4)
cat("PASS: real regional fits; initialization reused; two exact parallel count updates; eight full-posterior warm starts; final results fitted under iteration 2 priors.\n")
