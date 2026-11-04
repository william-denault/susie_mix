.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
distribution<-normalizePath("output/slider_inplace_2026-10-29",winslash="/")
root<-tempfile("inplace_integration_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(root)
bundle<-file.path(root,"bundle");dir.create(bundle)
stopifnot(all(file.copy(list.files(distribution,full.names=TRUE),bundle,recursive=TRUE)))
project<-file.path(root,"project");dir.create(project)
base<-file.path(project,"results_slide_prior_em");dir.create(base)
runtime<-file.path(project,"slider_inplace_runtime_2026-10-29")
meta<-file.path(base,"inplace_expansion_2026-10-29")
source(file.path(bundle,"maintenance.R"));load_helpers(bundle)
manifest<-data.frame(gene=c("GENEA","GENEB","XGENE"),chunk=1:3)
settings<-spe_settings();settings$package_version<-spe_check_packages()
weights<-as.data.frame(matrix(1/17,1,17,dimnames=list(NULL,spe_columns())))
history<-do.call(rbind,lapply(0:5,function(i) data.frame(iteration=i,source_iteration=i-1L,tissue="Tissue",weights,
 n_fits=0L,n_active=0L,free_mass=0,prior_max_change=0,m_step_q_gain=0,source_elbo=NA_real_,
 update_method="slider_prior_active_counts_v1",created_at="2026-10-28 UTC")))
em_atomic_write(history,file.path(base,"prior_history.csv"),csv=TRUE)
scope<-list(schema=1L,exclude=c("X","Y","MT"),genes=data.frame(gene=manifest$gene,chromosome=c("1","1","X")))
em_atomic_write(scope,file.path(base,"chromosome_scope.rds"))
make_out<-function(i,prior=rep(1/17,17),counts=rep(1/17,17)) {
 fit<-structure(list(alpha=matrix(1/3,3,3,dimnames=list(NULL,c("s1","s2","s3"))),
  alpha_delta=array(rep(counts/3,each=9),c(3L,3L,17L)),delta_grid=spe_grid(),delta_prior=prior,
  delta_prior_counts=matrix(rep(counts,each=3),3,17),V=rep(1,3),delta_forced=rep(FALSE,3),input_p=3L,
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
 status<-system2(python,c(shQuote(helper),shQuote(bundle),mode,shQuote(project),shQuote(as.character(c(...)))))
 if(status!=0) stop("Fixture Python failed: ",mode)
}
py("code")
dir.create(file.path(bundle,"checks"))
input<-file.path(base,"iteration_000/manifest.csv");info<-file.info(input)[,c("size","mtime"),drop=FALSE]
info$path<-input;rownames(info)<-NULL
map<-read.csv(file.path(bundle,"renamed_gene_map.csv"),stringsAsFactors=FALSE)
saveRDS(list(passed=TRUE,project=normalizePath(project,winslash="/"),map_md5=md5(file.path(bundle,"renamed_gene_map.csv")),
 input_info=info,packages=c(susieRSlidePrior=spe_check_packages()),datadir=root),file.path(bundle,"checks/preflight.rds"))
source(file.path(bundle,"initialize.R"))
stopifnot(inherits(try(initialize_six(bundle,project),silent=TRUE),"try-error"),!dir.exists(file.path(base,"iteration_006")))
cat("PASS: incomplete original 005 blocks initialization.\n")
parent<-sre_context(project,file.path(base,"iteration_005"));pooled5<-spp_complete(parent,spp_index(parent))
original_hashes<-tools::md5sum(c(sre_result(parent,c("GENEA","GENEB")),file.path(parent$dir,c("manifest.csv","priors.csv","settings.rds","grid.csv","recovery/parallel/pooled.rds"))))
original_history<-md5(file.path(base,"prior_history.csv"));original_scope<-md5(file.path(base,"chromosome_scope.rds"))
initialize_six(bundle,project)
# Owned setup can resume before workers start, without rewriting original data.
initialize_six(bundle,project);py("freeze")
for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R","slide_prior_parallel.R","prepare_slide_prior_em_iteration.R","cohort.R"))
 source(file.path(runtime,"script/scan_tissue_attempt",f))
ctx<-sre_context(project,file.path(base,"iteration_006"));spp_init(ctx)
stopifnot(nrow(ctx$manifest)==683L,length(ctx$scope$excluded)==1L,
 identical(original_hashes,tools::md5sum(names(original_hashes))),md5(file.path(base,"chromosome_scope.rds"))==original_scope,
 md5(file.path(meta,"prior_history.before_expansion.csv"))==original_history,
 max(abs(as.numeric(ctx$priors[1,spe_columns()])-as.numeric(pooled5$priors[1,spe_columns()])))<1e-12,
 identical(sre_previous(ctx,"GENEA"),make_out(5L)),is.null(sre_previous(ctx,map$gene[1])),is.null(spp_previous_cohort(ctx,map$gene[1])))
cat("PASS: 006 keeps the ordinary 005 prior update; 680 added IDs start cold; existing genes use unchanged 005 fits.\n")
# Refuse setup once any worker metadata/results exist.
stopifnot(inherits(try(initialize_six(bundle,project),silent=TRUE),"try-error"))
source(file.path(runtime,"retention.R"))
state<-file.path(base,"rolling_cleanup_inplace/before_006")
sr_prepare(project,state,6L);py("interrupt",6L);py("step",6L)
stopifnot(all(vapply(0:4,function(i) !length(list.files(file.path(base,sprintf("iteration_%03d/results",i)))),logical(1))),
 identical(original_hashes,tools::md5sum(names(original_hashes))))
spp_require_complete(project,parent$dir)
queue<-file.path(root,"queue.txt");writeLines(character(),queue)
source(file.path(runtime,"script/scan_tissue_attempt/rolling_hook.R"));sr_barrier(ctx,queue)
plan<-spp_plan(ctx,298L,queue);work<-read.csv(file.path(plan$batch,"manifest.csv"))
stopifnot(nrow(work)==298L,all(work$action=="fit"))
stopifnot(inherits(try(spp_complete(ctx,spp_index(ctx)),silent=TRUE),"try-error"))
cat("PASS: restartable cleanup retires 000..004 only; 005 certificate/warm starts remain valid; no premature M-step.\n")
for(g in c("GENEA","GENEB")) {em_atomic_write(make_out(6L),sre_result(ctx,g));spp_inspect(ctx,g,force=TRUE)}
skew<-c(.8,rep(.2/16,16))
for(j in seq_len(nrow(map))) {
 out<-if(j==1L) make_out(6L,counts=skew) else structure(list(),em_iteration=6L)
 em_atomic_write(out,sre_result(ctx,map$gene[j]));spp_inspect(ctx,map$gene[j],force=TRUE)
}
pool6<-spp_complete(ctx,spp_index(ctx))
stopifnot(sum(pool6$audit$issue=="fit")==3L,max(abs(as.numeric(pool6$counts[1,spe_columns()])-(6/17+3*skew)))<1e-12)
cat("PASS: added-gene counts enter the expanded E-step and M-step; old gene/tissue cohort is retained.\n")
for(i in 7:8) {
 launch<-spe_prepare_iteration(project,"new",completion_reader=spp_require_complete)
 stopifnot(basename(launch$iteration_dir)==sprintf("iteration_%03d",i))
 state<-file.path(base,"rolling_cleanup_inplace",sprintf("before_%03d",i))
 sr_prepare(project,state,i);py("interrupt",i);py("step",i)
 ctx<-sre_context(project,launch$iteration_dir);spp_init(ctx);sr_barrier(ctx,queue)
 spp_require_complete(project,file.path(base,sprintf("iteration_%03d",i-1L)))
 stopifnot(!is.null(sre_previous(ctx,map$gene[1])),!is.null(sre_previous(ctx,"GENEA")))
 prior<-as.numeric(ctx$priors[1,spe_columns()])
 if(i==7L) stopifnot(prior[1]>1/17,max(abs(prior-as.numeric(pool6$priors[1,spe_columns()])))<1e-12)
 for(g in setdiff(ctx$manifest$gene,ctx$scope$excluded)) {
  out<-sre_previous(ctx,g);attr(out,"em_iteration")<-as.integer(i)
  if(length(out)) out[[1]]$fit_slide_prior$delta_prior<-prior
  em_atomic_write(out,sre_result(ctx,g));spp_inspect(ctx,g,force=TRUE)
 }
 spp_complete(ctx,spp_index(ctx))
 stopifnot(all(vapply(0:(i-2L),function(k) !length(list.files(file.path(base,sprintf("iteration_%03d/results",k)))),logical(1))),
   length(list.files(file.path(base,sprintf("iteration_%03d/results",i-1L))))>0L,
   md5(file.path(meta,"prior_history.before_expansion.csv"))==original_history)
 cat("PASS: before",i,"cleanup preserves all history/counts, keeps previous warm starts, and maintains the expanded cohort.\n")
}
for(f in c("em_utils.R","slide_prior_em_utils.R"))
 stopifnot(md5(file.path(runtime,"script/scan_tissue_attempt",f))==md5(file.path(bundle,"baseline/script/scan_tissue_attempt",f)))
writeLines(root,"tmp/slider_inplace_fixture_path.txt")
cat("PASS: numerical EM helpers unchanged; same results folder through 008. FIXTURE:",root,"\n")
