.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
code<-normalizePath("recovery/slider_online_extend_024_038_2026-11-09/runtime",winslash="/")
audit<-normalizePath("recovery/slider_online_extend_024_038_2026-11-09/tests",winslash="/")
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
genes<-c("T","F",paste0("REAL",3:4))
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
        initial_fits=4L,cohorts=setNames(rep(list(spp_cohort(out)),4),genes),
        counts=4*sm$counts,initial_priors=mb_weights(4*sm$counts,priors))
mb_save(m,file.path(seed,"metadata.rds"));paths<-file.path(seed,"results",paste0(genes,".rds"));info<-file.info(paths)
mb_save(data.frame(gene=genes,size=info$size,mtime=as.numeric(info$mtime),md5=unname(tools::md5sum(paths))),file.path(seed,"expected.csv"),TRUE)
dir.create(file.path(branch,"plan"));rows<-list();sid<-0L
for(iteration in 15:38) {
  online<-iteration %in% c(15:20,24:33)
  groups<-if(online)2L else 1L
  rho<-if(!online)1 else if(iteration<=16).3 else if(iteration<=20).5 else if(iteration<=25).8 else if(iteration<=28).5 else .3
  chosen<-sample(genes)
  batches<-split(chosen,rep(seq_len(groups),length.out=length(genes)))
  for(batch in seq_len(groups)) {
    sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(batches[[batch]],file.path(branch,gf))
    group<-sprintf("../../results_slide_prior_em/iteration_%03d",iteration)
    rows[[sid]]<-data.frame(stage=sid,iteration=iteration,phase=if(online)"minibatch"else"full",pass=iteration-14,batch=batch,rho=rho,
      group=group,results_dir=paste0(group,"/results"),directory=paste0(group,if(online)sprintf("/batches/batch_%02d",batch)else"/batches/full_step"),genes_file=gf)
  }
}
plan<-do.call(rbind,rows);mb_save(plan,file.path(branch,"stages.csv"),TRUE);mb_initialize(branch,seed);mb_history(branch)
python<-"C:/Users/willi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
bridge<-function(mode,index) {
  log<-file.path(branch,paste0("audit_",mode,"_",index,".log"))
  status<-system2(python,shQuote(c(file.path(audit,"bridge.py"),code,branch,mode,as.character(index))),stdout=log,stderr=log)
  if(status!=0)stop(paste(readLines(log),collapse="\n"))
}
bridge("init",1L)
calls<-0L;expected_iteration<-setNames(rep(14L,4),genes)
real_gene<-function(target_gene,tissue_priors,previous_result,...) {
  stopifnot(attr(previous_result,"em_iteration")==expected_iteration[[target_gene]])
  calls<<-calls+1L
  args<-a;args$delta_prior<-spe_prior(tissue_priors,"Brain")
  warm<-previous_result$Brain$fit_slide_prior;warm$call<-NULL;args$model_init<-warm
  list(Brain=list(fit_slide_prior=do.call(susieRSlidePrior::susie,args),sample_ids=as.character(1:160)))
}
for(i in seq_len(nrow(plan))) {
  bridge("own",i);mb_prepare(branch);mb_history(branch)
  d<-file.path(branch,plan$directory[i]);x<-readRDS(file.path(d,"input.rds"));before<-readRDS(file.path(branch,"state.rds"))
  stopifnot(all(x$warm$iteration==x$iteration-1L))
  mb_save(before$latest,file.path(branch,"audit_warm.csv"),TRUE);bridge("cleanup",i)
  for(g in x$genes)mb_worker(branch,d,g,real_gene)
  counts<-Reduce(`+`,lapply(x$genes,function(g)readRDS(file.path(d,"summaries",paste0(g,".rds")))$counts))
  expected<-(1-x$rho)*before$S+x$rho*counts/length(x$genes)
  mb_commit(branch,d);mb_history(branch);st<-readRDS(file.path(branch,"state.rds"))
  stopifnot(max(abs(st$S-expected))<1e-12,all(st$latest$iteration %in% c(x$iteration-1L,x$iteration)))
  expected_iteration[x$genes]<-x$iteration
  if(i==max(which(plan$iteration==x$iteration))) {
    stopifnot(all(st$latest$iteration==x$iteration),all(file.exists(st$latest$path)),all(file.exists(st$latest$priors)))
    cat("PASS real fits through iteration",x$iteration,"with actual cleanup and per-gene batch priors.\n")
  }
}
stopifnot(calls==96L,all(st$latest$iteration==38L),st$stage==40L)
first_batch_priors<-lapply(plan$directory[1:2],function(d)as.matrix(mb_read(file.path(branch,d,"priors.csv"))[,spe_columns()]))
stopifnot(any(vapply(first_batch_priors[-1],function(p)max(abs(p-first_batch_priors[[1]]))>1e-10,logical(1))))
cat("PASS: 96 actual susieRSlidePrior fits across 015-038, 40 reduced-cohort updates, differing batch priors, T/F identifiers, cleanup, and exact per-gene warm-start provenance.\n")

summary<-mb_read(file.path(root,"iteration_summary.csv"))
stopifnot(all(summary$completed[summary$iteration>=15L]),all(summary$status[summary$iteration>=15L]=="completed"),all(summary$rho[summary$iteration %in% 24:25]==.8),all(summary$rho[summary$iteration %in% 26:28]==.5),all(summary$rho[summary$iteration %in% 29:33]==.3))
