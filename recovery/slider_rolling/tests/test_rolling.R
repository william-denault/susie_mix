# Synthetic integration: real validation/M-step, hard links, interrupted cleanup,
# repaired certificates, and three successive iteration transitions.
.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
distribution <- normalizePath("output/slider_rolling_2026-10-28",winslash="/")
root <- tempfile("rolling_integration_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(root)
bundle <- file.path(root,"bundle");dir.create(bundle)
stopifnot(all(file.copy(list.files(distribution,full.names=TRUE),bundle,recursive=TRUE)))
project <- file.path(root,"project");dir.create(project)
base <- file.path(project,"results_slide_prior_em");dir.create(base)
branch <- file.path(project,"slider_expanded_040_2026-10-28")
source(file.path(bundle,"maintenance.R"));load_helpers(bundle)
manifest <- data.frame(gene=c("GENEA","GENEB","XGENE"),chunk=1:3)
settings <- spe_settings();settings$package_version<-spe_check_packages()
weights <- as.data.frame(matrix(1/17,1,17,dimnames=list(NULL,spe_columns())))
history <- do.call(rbind,lapply(0:5,function(i) data.frame(iteration=i,source_iteration=i-1L,tissue="Tissue",weights,
  n_fits=0L,n_active=0L,free_mass=0,prior_max_change=0,m_step_q_gain=0,source_elbo=NA_real_,
  update_method="slider_prior_active_counts_v1",created_at="2026-10-28 UTC")))
em_atomic_write(history,file.path(base,"prior_history.csv"),csv=TRUE)
scope <- list(schema=1L,exclude=c("X","Y","MT"),genes=data.frame(gene=manifest$gene,chromosome=c("1","1","X")))
em_atomic_write(scope,file.path(base,"chromosome_scope.rds"))
make_out <- function(i) {
  fit <- structure(list(alpha=matrix(1/3,3,3,dimnames=list(NULL,c("s1","s2","s3"))),
    alpha_delta=array(1/51,c(3L,3L,17L)),delta_grid=spe_grid(),delta_prior=rep(1/17,17),
    delta_prior_counts=matrix(1/17,3,17),V=rep(1,3),delta_forced=rep(FALSE,3),input_p=3L,
    elbo=-12,converged=TRUE,mu_grid=array(.1,c(3L,3L,17L)),mu2_grid=array(.2,c(3L,3L,17L))),class=c("susie_slide","susie"))
  structure(list(Tissue=list(fit_slide_prior=fit,sample_ids=c("A","B"))),em_iteration=as.integer(i))
}
for(i in 0:5) {
  d<-file.path(base,sprintf("iteration_%03d",i));dir.create(file.path(d,"results"),recursive=TRUE)
  em_atomic_write(history[history$iteration==i,],file.path(d,"priors.csv"),csv=TRUE)
  em_atomic_write(history[history$iteration<=i,],file.path(d,"prior_history.csv"),csv=TRUE)
  em_atomic_write(settings,file.path(d,"settings.rds"));em_atomic_write(manifest,file.path(d,"manifest.csv"),csv=TRUE)
  em_atomic_write(data.frame(weight_column=spe_columns(),delta=spe_grid()),file.path(d,"grid.csv"),csv=TRUE)
  for(g in c("GENEA","GENEB")) em_atomic_write(make_out(i),file.path(d,"results",paste0(g,".rds")))
  ctx<-sre_context(project,d);spp_init(ctx)
  for(g in c("GENEA","GENEB")) spp_inspect(ctx,g,force=TRUE)
  if(i<5) spp_complete(ctx,spp_index(ctx))
}
python<-Sys.getenv("PYTHON");stopifnot(nzchar(python))
helper<-file.path(bundle,"tests/fixture_helper.py")
py<-function(mode,...) {
  status<-system2(python,c(shQuote(helper),shQuote(bundle),mode,shQuote(branch),shQuote(as.character(c(...)))))
  if(status!=0) stop("Fixture Python failed: ",mode)
}
py("code")
dir.create(file.path(bundle,"checks"))
input<-file.path(base,"iteration_000/manifest.csv");info<-file.info(input)[,c("size","mtime"),drop=FALSE]
info$path<-input;rownames(info)<-NULL
map<-read.csv(file.path(bundle,"renamed_gene_map.csv"),stringsAsFactors=FALSE)
saveRDS(list(passed=TRUE,project=normalizePath(project,winslash="/"),map_md5=md5(file.path(bundle,"renamed_gene_map.csv")),
  input_info=info,packages=c(susieRSlidePrior=spe_check_packages()),datadir=root),file.path(bundle,"checks/preflight.rds"))
source(file.path(bundle,"expanded.R"))
stopifnot(inherits(try(seed_expanded(bundle,project,branch),silent=TRUE),"try-error"))
cat("PASS: unfinished original 005 blocks expansion.\n")
parent<-sre_context(project,file.path(base,"iteration_005"));spp_complete(parent,spp_index(parent))
original_hashes<-tools::md5sum(sre_result(parent,c("GENEA","GENEB")))
original_history<-md5(file.path(base,"prior_history.csv"))
seed_expanded(bundle,project,branch)
py("freeze",project)
for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R","slide_prior_parallel.R","prepare_slide_prior_em_iteration.R"))
  source(file.path(branch,"script/scan_tissue_attempt",f))
expanded_root<-file.path(branch,"results_slide_prior_em")
ctx<-sre_context(branch,file.path(expanded_root,"iteration_005"))
stopifnot(identical(original_hashes,tools::md5sum(names(original_hashes))),
  md5(file.path(expanded_root,"prior_history.csv"))==original_history,
  md5(file.path(expanded_root,"original_prior_history.csv"))==original_history)
queue<-file.path(root,"queue.txt");writeLines(character(),queue)
plan<-spp_plan(ctx,298L,queue);work<-read.csv(file.path(plan$batch,"manifest.csv"))
stopifnot(nrow(work)==298L,all(work$gene %in% map$gene),all(work$action=="fit"))
stopifnot(inherits(try(spp_complete(ctx,spp_index(ctx)),silent=TRUE),"try-error"))
unlink(spp_summary_file(ctx,"GENEA"));stopifnot(spp_inspect(ctx,"GENEA",force=TRUE)$status=="valid")
for(j in seq_len(nrow(map))) {
  out<-if(j==1L) make_out(5L) else structure(list(),em_iteration=5L)
  em_atomic_write(out,sre_result(ctx,map$gene[j]));spp_inspect(ctx,map$gene[j],force=TRUE)
}
pooled<-spp_complete(ctx,spp_index(ctx))
stopifnot(sum(pooled$audit$issue=="fit")==3L,max(abs(as.numeric(pooled$counts[1,spe_columns()])-9/17))<1e-12)
cat("PASS: original bytes/priors/history reused; only added IDs fitted at 005; added counts enter the M-step.\n")
launch<-spe_prepare_iteration(branch,"new",completion_reader=spp_require_complete)
stopifnot(basename(launch$iteration_dir)=="iteration_006")
source(file.path(branch,"retention.R"))
state<-file.path(expanded_root,"rolling_cleanup/before_006")
sr_prepare(branch,state,6L)
py("interrupt",6L)
py("step",6L)
stopifnot(!length(list.files(base,"\\.rds$",recursive=TRUE,full.names=TRUE)[grepl("/results/",list.files(base,"\\.rds$",recursive=TRUE,full.names=TRUE),fixed=TRUE)]))
stopifnot(md5(file.path(base,"prior_history.csv"))==original_history)
spp_require_complete(branch,file.path(expanded_root,"iteration_005"))
cat("PASS: source outputs retire after hard-link transfer; interrupted deletion resumes; all source history remains.\n")
for(i in 6:7) {
  ctx<-sre_context(branch,file.path(expanded_root,sprintf("iteration_%03d",i)));spp_init(ctx)
  stopifnot(identical(sre_previous(ctx,"GENEA"),make_out(i-1L)))
  for(g in setdiff(ctx$manifest$gene,ctx$scope$excluded)) {
    out<-sre_previous(ctx,g);attr(out,"em_iteration")<-as.integer(i)
    em_atomic_write(out,sre_result(ctx,g));spp_inspect(ctx,g,force=TRUE)
  }
  # Once completed, the cleanup receipt makes repeats harmless even with current results.
  source(file.path(branch,"script/scan_tissue_attempt/rolling_hook.R"));sr_barrier(ctx,queue)
  spp_complete(ctx,spp_index(ctx))
  launch<-spe_prepare_iteration(branch,"new",completion_reader=spp_require_complete)
  current<-i+1L;state<-file.path(expanded_root,"rolling_cleanup",sprintf("before_%03d",current))
  sr_prepare(branch,state,current)
  py("interrupt",current)
  py("step",current)
  sr_barrier(sre_context(branch,launch$iteration_dir),queue)
  spp_require_complete(branch,ctx$dir)
  stopifnot(spp_inspect(ctx,"GENEA")$status=="valid")
  next_ctx<-sre_context(branch,launch$iteration_dir)
  stopifnot(identical(sre_previous(next_ctx,"GENEA"),make_out(i)))
  old<-file.path(expanded_root,sprintf("iteration_%03d",5:(current-2L)),"results")
  stopifnot(all(vapply(old,function(d) length(list.files(d))==0L,logical(1))))
  stopifnot(file.exists(file.path(expanded_root,sprintf("iteration_%03d",current-2L),"recovery/parallel/pooled.rds")))
  cat("PASS: cleanup before",current,"preserves counts/history, repairs certificate, and keeps",i,"as a valid warm start.\n")
}
# Fit-settings and numerical EM helpers remain the deployed baseline.
for(f in c("em_utils.R","slide_prior_em_utils.R"))
  stopifnot(md5(file.path(branch,"script/scan_tissue_attempt",f))==md5(file.path(bundle,"baseline/script/scan_tissue_attempt",f)))
stopifnot(md5(file.path(base,"prior_history.csv"))==original_history,
          md5(file.path(expanded_root,"original_prior_history.csv"))==original_history)
cat("PASS: unchanged EM numerical helpers and fitting settings; iteration numbering continues 005,006,007,008.\n")
writeLines(root,"tmp/slider_rolling_fixture_path.txt")
cat("FIXTURE:",root,"\n")
