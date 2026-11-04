.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
bundle<-normalizePath("output/slider_partial_2026-10-29",winslash="/")
source(file.path(bundle,"maintenance.R"));load_helpers(bundle)
root<-tempfile("partial_metadata_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(root)
project<-file.path(root,"project");dir.create(project)
base<-file.path(project,"results_slide_prior_em");dir.create(base)
original<-file.path(root,"original_state");dir.create(original)
state<-file.path(root,"final_state")
manifest<-data.frame(gene=c("GENEA","GENEB","XGENE"),chunk=1:3)
settings<-spe_settings();settings$package_version<-spe_check_packages()
weights<-as.data.frame(matrix(1/17,1,17,dimnames=list(NULL,spe_columns())))
history<-do.call(rbind,lapply(0:2,function(i) data.frame(iteration=i,source_iteration=i-1L,tissue="Tissue",weights,
  n_fits=0L,n_active=0L,free_mass=0,prior_max_change=0,m_step_q_gain=0,source_elbo=NA_real_,
  update_method="slider_prior_active_counts_v1",created_at="2026-10-29 UTC")))
em_atomic_write(history,file.path(base,"prior_history.csv"),csv=TRUE)
scope<-list(schema=1L,exclude=c("X","Y","MT"),genes=data.frame(gene=manifest$gene,chromosome=c("1","1","X")))
em_atomic_write(scope,file.path(base,"chromosome_scope.rds"))
make_out<-function(i,previous=NULL) {
  fit<-structure(list(alpha=matrix(1/3,3,3,dimnames=list(NULL,c("s1","s2","s3"))),
    alpha_delta=array(1/51,c(3L,3L,17L)),delta_grid=spe_grid(),delta_prior=rep(1/17,17),
    delta_prior_counts=matrix(1/17,3,17),V=rep(1,3),delta_forced=rep(FALSE,3),input_p=3L,
    elbo=-12,converged=TRUE,mu_grid=array(.1,c(3L,3L,17L)),mu2_grid=array(.2,c(3L,3L,17L)),
    call=as.call(c(list(as.name("f"),X=matrix(seq_len(4000),100)),if(!is.null(previous)) list(model_init=previous$Tissue$fit_slide_prior)))),class=c("susie_slide","susie"))
  structure(list(Tissue=list(fit_slide_prior=fit,sample_ids=c("A","B"))),em_iteration=as.integer(i))
}
previous<-NULL
for(i in 0:2) {
  d<-file.path(base,sprintf("iteration_%03d",i));dir.create(file.path(d,"results"),recursive=TRUE)
  em_atomic_write(history[history$iteration==i,],file.path(d,"priors.csv"),csv=TRUE)
  em_atomic_write(history[history$iteration<=i,],file.path(d,"prior_history.csv"),csv=TRUE)
  em_atomic_write(settings,file.path(d,"settings.rds"));em_atomic_write(manifest,file.path(d,"manifest.csv"),csv=TRUE)
  em_atomic_write(data.frame(weight_column=spe_columns(),delta=spe_grid()),file.path(d,"grid.csv"),csv=TRUE)
  out<-make_out(i,previous)
  for(g in if(i<2) c("GENEA","GENEB") else "GENEA") em_atomic_write(out,file.path(d,"results",paste0(g,".rds")))
  ctx<-sre_context(project,d);spp_init(ctx)
  if(i<2) { for(g in c("GENEA","GENEB")) spp_inspect(ctx,g,force=TRUE);spp_complete(ctx,spp_index(ctx)) }
  previous<-out
}
# Match RCC: older full results have retired, but history directories remain.
retired<-file.path(root,"retired");dir.create(retired)
retire_prepare(bundle,project,retired,0L,2L)
unlink(list.files(file.path(base,"iteration_000/results"),full.names=TRUE))
retire_finish(bundle,project,retired)
cert0<-file.path(base,"iteration_000/recovery/parallel/COMPLETE.rds")
stopifnot(file.rename(cert0,paste0(cert0,".retired")))
writeLines('{}',file.path(base,"iteration_000/storage_retired.json"))
snapshot(bundle,project,original,2L)
preserved<-tools::md5sum(history_files(project))
untouched<-file.path(base,c("iteration_001/results/GENEB.rds","iteration_002/results/GENEA.rds"))
before<-tools::md5sum(untouched)
target<-file.path(base,"iteration_001/results/GENEA.rds");expected<-strip_calls(readRDS(target))$out
py<-function(mode) {
  status<-system2(Sys.getenv("PYTHON"),c(shQuote(file.path(bundle,"tests/fixture_helper.py")),shQuote(bundle),shQuote(root),mode))
  if(status!=0L) stop("Python fixture failed: ",mode)
}
py("compact_one");py("finalize");py("finalize")
stopifnot(identical(readRDS(target),expected),identical(before,tools::md5sum(untouched)),
          identical(preserved,tools::md5sum(names(preserved))))
invisible(spp_require_complete(project,file.path(base,"iteration_001")))
report<-read.csv(file.path(state,"validation_report.csv"))
stopifnot(report$deferred_file_validation[report$iteration==2L]==1L)
cat("PASS: mixed compact/original checkpoint validates; untouched files, prior history and counts unchanged.\n")
ctx<-sre_context(project,file.path(base,"iteration_002"))
queue<-file.path(root,"queue.txt");writeLines(character(),queue)
plan<-spp_plan(ctx,298L,queue);work<-read.csv(file.path(plan$batch,"manifest.csv"))
stopifnot(identical(work$gene,"GENEA"),identical(work$action,"validate"))
spp_worker_gene(ctx,"GENEA","validate",run_gene=function(...) stop("A completed gene must not be refitted"))
plan<-spp_plan(ctx,298L,queue);work<-read.csv(file.path(plan$batch,"manifest.csv"))
stopifnot(identical(work$gene,"GENEB"),identical(work$action,"fit"))
stopifnot(identical(sre_previous(ctx,"GENEA"),expected))
stopifnot(md5(file.path(project,"script/scan_tissue_attempt/workhorse_slide_prior_em.R"))==
          md5(file.path(bundle,"patch/script/scan_tissue_attempt/workhorse_slide_prior_em.R")))
cat("PASS: saved current gene is validated without refitting; only missing gene is scheduled; call-stripping worker installed.\n")
# Repeat metadata repair after a hypothetical interruption, using the same snapshot.
source(file.path(bundle,"partial_finalize.R"));finalize_partial(bundle,project,state)
stopifnot(identical(preserved,tools::md5sum(names(preserved))))
cat("PASS: metadata repair is restartable.\n")
# No receipt plus changed result must fail, not silently adopt the new bytes.
changed<-file.path(base,"iteration_001/results/GENEB.rds")
writeBin(charToRaw("external change"),changed)
stopifnot(inherits(try(finalize_partial(bundle,project,state),silent=TRUE),"try-error"))
cat("PASS: unjournaled modifications are rejected.\n")
cat("FIXTURE:",root,"\n")
