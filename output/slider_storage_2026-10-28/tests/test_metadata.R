.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
bundle <- normalizePath("output/slider_storage_2026-10-28",winslash="/")
source(file.path(bundle,"maintenance.R"));load_helpers(bundle)
root <- tempfile("storage_metadata_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(root)
project <- file.path(root,"project");dir.create(project)
base <- file.path(project,"results_slide_prior_em");dir.create(base)
state <- file.path(root,"state");dir.create(state)
manifest <- data.frame(gene=c("GENEA","GENEB","XGENE"),chunk=1:3)
settings <- spe_settings();settings$package_version <- spe_check_packages()
weights <- as.data.frame(matrix(1/17,1,17,dimnames=list(NULL,spe_columns())))
history <- do.call(rbind,lapply(0:2,function(i) data.frame(iteration=i,source_iteration=i-1L,tissue="Tissue",weights,
  n_fits=0L,n_active=0L,free_mass=0,prior_max_change=0,m_step_q_gain=0,source_elbo=NA_real_,
  update_method="slider_prior_active_counts_v1",created_at="2026-10-28 UTC")))
em_atomic_write(history,file.path(base,"prior_history.csv"),csv=TRUE)
scope <- list(schema=1L,exclude=c("X","Y","MT"),genes=data.frame(gene=manifest$gene,chromosome=c("1","1","X")))
em_atomic_write(scope,file.path(base,"chromosome_scope.rds"))
make_out <- function(i,prior=NULL) {
  alpha <- matrix(1/3,3,3,dimnames=list(NULL,c("s1","s2","s3")))
  fit <- structure(list(alpha=alpha,alpha_delta=array(1/51,c(3L,3L,17L)),
    delta_grid=spe_grid(),delta_prior=rep(1/17,17),delta_prior_counts=matrix(1/17,3,17),
    V=rep(1,3),delta_forced=rep(FALSE,3),input_p=3L,elbo=-12,converged=TRUE,
    mu_grid=array(.1,c(3L,3L,17L)),mu2_grid=array(.2,c(3L,3L,17L)),
    call=as.call(c(list(as.name("f"),X=matrix(seq_len(4000),100)),if(!is.null(prior)) list(model_init=prior$Tissue$fit_slide_prior)))),class=c("susie_slide","susie"))
  structure(list(Tissue=list(fit_slide_prior=fit,sample_ids=c("A","B"))),em_iteration=as.integer(i))
}
previous <- NULL
for(i in 0:2) {
  d <- file.path(base,sprintf("iteration_%03d",i));dir.create(file.path(d,"results"),recursive=TRUE)
  em_atomic_write(history[history$iteration==i,],file.path(d,"priors.csv"),csv=TRUE)
  em_atomic_write(history[history$iteration<=i,],file.path(d,"prior_history.csv"),csv=TRUE)
  em_atomic_write(settings,file.path(d,"settings.rds"));em_atomic_write(manifest,file.path(d,"manifest.csv"),csv=TRUE)
  em_atomic_write(data.frame(weight_column=spe_columns(),delta=spe_grid()),file.path(d,"grid.csv"),csv=TRUE)
  out <- make_out(i,previous)
  for(g in if(i<2) c("GENEA","GENEB") else "GENEA") em_atomic_write(out,file.path(d,"results",paste0(g,".rds")))
  ctx <- sre_context(project,d);spp_init(ctx)
  if(i<2) {
    for(g in c("GENEA","GENEB")) spp_inspect(ctx,g,force=TRUE)
    spp_complete(ctx,spp_index(ctx))
  }
  previous <- out
}
snapshot(bundle,project,state,2L)
s <- readRDS(file.path(state,"snapshot.rds"));preserved <- s$preserved
files <- list.files(base,"\\.rds$",recursive=TRUE,full.names=TRUE)
files <- files[grepl("/results/",files,fixed=TRUE)]
receipts <- list()
for(f in files) {
  before <- md5(f);expected <- strip_calls(readRDS(f))$out
  temp <- file.path(root,"compact.rds");report <- file.path(root,"report.csv")
  compact_file(f,temp,report)
  stopifnot(identical(readRDS(temp),expected))
  stopifnot(file.copy(temp,f,overwrite=TRUE))
  receipts[[length(receipts)+1L]] <- data.frame(relative=substring(f,nchar(project)+2L),old_md5=before,new_md5=md5(f))
}
write.csv(do.call(rbind,receipts),file.path(state,"receipts.csv"),row.names=FALSE)
finalize(bundle,project,state)
stopifnot(identical(preserved,tools::md5sum(names(preserved))))
for(i in 0:1) spp_require_complete(project,file.path(base,sprintf("iteration_%03d",i)))
ctx <- sre_context(project,file.path(base,"iteration_002"))
q <- file.path(root,"queue.txt");writeLines(character(),q)
plan <- spp_plan(ctx,298L,q)
work <- read.csv(file.path(plan$batch,"manifest.csv"))
stopifnot(identical(work$gene,"GENEB"),identical(work$action,"fit"))
# Repeat metadata finalization after interruption: identical priors, same reuse.
finalize(bundle,project,state)
cat("PASS: compaction preserves all scientific fields and history; only missing GENEB is scheduled.\n")
retire_state <- file.path(root,"retire");dir.create(retire_state)
stopifnot(inherits(try(retire_prepare(bundle,project,retire_state,1L,2L),silent=TRUE),"try-error"))
stopifnot(length(list.files(file.path(base,"iteration_001/results")))==2L)
retire_prepare(bundle,project,retire_state,0L,2L)
unlink(list.files(file.path(base,"iteration_000/results"),full.names=TRUE))
retire_finish(bundle,project,retire_state)
spp_require_complete(project,file.path(base,"iteration_001"))
stopifnot(spp_inspect(ctx,"GENEA")$status=="valid",identical(preserved,tools::md5sum(names(preserved))))
cat("PASS: retiring old full results leaves warm-start and current summaries valid.\n")
# Export fixture paths for the Python/R expansion integration check.
writeLines(project,file.path("tmp","slider_storage_fixture_path.txt"))
cat("FIXTURE:",root,"\n")
