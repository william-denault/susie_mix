.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
code<-normalizePath("output/slider_online_queue_patch_2026-11-06/runtime",winslash="/")
audit<-normalizePath("recovery/slider_transition_audit_2026-11-06",winslash="/")
source(file.path(code,"engine.R"));mb_load(code)
Sys.setenv(RSCRIPT=file.path(R.home("bin"),"Rscript"),R_LIBS_USER=paste(.libPaths(),collapse=.Platform$path.sep))
set.seed(104);X<-matrix(rbinom(160*6,2,.35),160,6);X[,6]<-0;X[1:8,6]<-1
colnames(X)<-paste0("snp",1:6);y<-X[,1]+.5*(X[,1]==1)+rnorm(160,sd=.4)
priors<-data.frame(tissue="Brain",matrix(1/17,1,17,dimnames=list(NULL,spe_columns())))
a<-list(X=X,y=y,L=10,min_obs=5,standardize=FALSE,estimate_prior_method="EM",
        delta_grid=spe_grid(),delta_prior=spe_prior(priors,"Brain"),max_iter=5000,tol=1e-5)
fit<-do.call(susieRSlidePrior::susie,a)
out<-structure(list(Brain=list(fit_slide_prior=fit,sample_ids=as.character(1:160))),em_iteration=14L)
sre_validate(out,list(iteration=14L,priors=priors,settings=spe_settings()))
genes<-c("T","F",paste0("REAL",3:9))
project<-tempfile("transition_real_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(project)
root<-file.path(project,"results_slide_prior_em");dir.create(root)
branch<-file.path(project,"workflows/slider_online_015_025");dir.create(branch,recursive=TRUE)
seed<-file.path(branch,"seed/snapshot_001");dir.create(file.path(seed,"results"),recursive=TRUE)
for(g in genes) {
  mb_save(out,file.path(seed,"results",paste0(g,".rds")))
  mb_save(out,file.path(root,"iteration_014/results",paste0(g,".rds")))
}
mb_save(priors,file.path(seed,"source_priors.csv"),TRUE)
old<-data.frame(iteration=14L,source_iteration=13L,priors)
mb_save(old,file.path(seed,"source_prior_history.csv"),TRUE);mb_save(old,file.path(root,"prior_history.csv"),TRUE)
sm<-spe_summarize_gene(out,"G",priors)
m<-list(settings=spe_settings(),packages=mb_packages(),genes=genes,datadir="unused",source_project=project,
        map=mb_read(file.path(code,"renamed_gene_map.csv"))[FALSE,],source_iteration=14L,
        initial_fits=9L,cohorts=setNames(rep(list(spp_cohort(out)),9),genes),
        counts=9*sm$counts,initial_priors=mb_weights(9*sm$counts,priors))
mb_save(m,file.path(seed,"metadata.rds"));paths<-file.path(seed,"results",paste0(genes,".rds"));info<-file.info(paths)
mb_save(data.frame(gene=genes,size=info$size,mtime=as.numeric(info$mtime),md5=unname(tools::md5sum(paths))),file.path(seed,"expected.csv"),TRUE)
dir.create(file.path(branch,"plan"));rows<-list();sid<-0L
for(pass in 1:6) {
  chosen<-sample(genes)
  for(batch in 1:9) {
    sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(chosen[batch],file.path(branch,gf))
    group<-sprintf("../../results_slide_prior_em/iteration_%03d",14+pass)
    rows[[sid]]<-data.frame(stage=sid,iteration=14+pass,phase="minibatch",pass=pass,batch=batch,rho=.3,
      group=group,results_dir=paste0(group,"/results"),directory=paste0(group,sprintf("/batches/batch_%02d",batch)),genes_file=gf)
  }
}
for(step in 1:5) {
  sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(genes,file.path(branch,gf))
  group<-sprintf("../../results_slide_prior_em/iteration_%03d",20+step)
  rows[[sid]]<-data.frame(stage=sid,iteration=20+step,phase="full",pass=step,batch=1,rho=1,
    group=group,results_dir=paste0(group,"/results"),directory=paste0(group,"/batches/full_step"),genes_file=gf)
}
plan<-do.call(rbind,rows);mb_save(plan,file.path(branch,"stages.csv"),TRUE);mb_initialize(branch,seed);mb_history(branch)
python<-"C:/Users/willi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
bridge<-function(mode,index) {
  log<-file.path(branch,paste0("audit_",mode,"_",index,".log"))
  status<-system2(python,shQuote(c(file.path(audit,"bridge.py"),code,branch,mode,as.character(index))),stdout=log,stderr=log)
  if(status!=0)stop(paste(readLines(log),collapse="\n"))
}
bridge("init",1L)
calls<-0L;expected_iteration<-setNames(rep(14L,9),genes)
real_gene<-function(target_gene,tissue_priors,previous_result,...) {
  stopifnot(attr(previous_result,"em_iteration")==expected_iteration[[target_gene]])
  calls<<-calls+1L
  args<-a;args$delta_prior<-spe_prior(tissue_priors,"Brain")
  warm<-previous_result$Brain$fit_slide_prior;warm$call<-NULL;args$model_init<-warm
  list(Brain=list(fit_slide_prior=do.call(susieRSlidePrior::susie,args),sample_ids=as.character(1:160)))
}
for(i in 1:59) {
  bridge("own",i);mb_prepare(branch);mb_history(branch)
  d<-file.path(branch,plan$directory[i]);x<-readRDS(file.path(d,"input.rds"));before<-readRDS(file.path(branch,"state.rds"))
  stopifnot(all(x$warm$iteration==x$iteration-1L))
  mb_save(before$latest,file.path(branch,"audit_warm.csv"),TRUE);bridge("cleanup",i)
  for(g in x$genes)mb_worker(branch,d,g,real_gene)
  counts<-Reduce(`+`,lapply(x$genes,function(g)readRDS(file.path(d,"summaries",paste0(g,".rds")))$counts))
  expected<-if(i<=54) .7*before$S+.3*counts/length(x$genes) else counts/length(genes)
  mb_commit(branch,d);mb_history(branch);st<-readRDS(file.path(branch,"state.rds"))
  stopifnot(max(abs(st$S-expected))<1e-12,all(st$latest$iteration %in% c(x$iteration-1L,x$iteration)))
  expected_iteration[x$genes]<-x$iteration
  if(i%%9==0L || i>54L) {
    stopifnot(all(st$latest$iteration==x$iteration),all(file.exists(st$latest$path)),all(file.exists(st$latest$priors)))
    cat("PASS real fits through iteration",x$iteration,"with actual cleanup and per-gene batch priors.\n")
  }
}
stopifnot(calls==99L,all(st$latest$iteration==25L),st$stage==59L)
first_batch_priors<-lapply(plan$directory[1:9],function(d)as.matrix(mb_read(file.path(branch,d,"priors.csv"))[,spe_columns()]))
stopifnot(any(vapply(first_batch_priors[-1],function(p)max(abs(p-first_batch_priors[[1]]))>1e-10,logical(1))))
cat("PASS: 99 actual susieRSlidePrior fits across 015-025, all 59 updates, differing batch priors, T/F identifiers, cleanup, and exact per-gene warm-start provenance.\n")
