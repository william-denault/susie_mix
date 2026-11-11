.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
code<-normalizePath("recovery/slider_online_extend_024_038_2026-11-09/runtime",winslash="/"); oldcode<-normalizePath("output/slider_online_rho05_18_patch_2026-11-06/runtime",winslash="/")
source(file.path(oldcode,"engine.R"));mb_load(oldcode)
project<-tempfile("transition_extended_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(project)
root<-file.path(project,"results_slide_prior_em");dir.create(root)
branch<-file.path(project,"workflows/slider_online_015_025");dir.create(branch,recursive=TRUE)
seed<-file.path(branch,"seed/snapshot_001");dir.create(file.path(seed,"results"),recursive=TRUE)
fails<-function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
map<-mb_read(file.path(code,"renamed_gene_map.csv"))[1:2,]
genes<-c(map$gene,"T","F",sprintf("GENE%02d",5:61));tissues<-c("Tissue","ZeroMass")
priors<-data.frame(tissue=tissues,matrix(1/17,2,17,dimnames=list(NULL,spe_columns())))
make_out<-function(iter,prior,g) {
  out<-lapply(tissues,function(t) {
    p<-spe_prior(prior,t);q<-p*(.8+seq_len(17)/80+(g%%4)*rev(seq_len(17))/100)
    q<-q/sum(q);joint<-array(0,c(3L,3L,17L))
    for(k in 1:17) joint[,c(1,3),k]<-q[k]/3
    joint[,2,9]<-1/3
    fit<-structure(list(alpha=matrix(1/3,3,3,dimnames=list(NULL,c("s1","s2","s3"))),
      alpha_delta=joint,delta_grid=spe_grid(),delta_prior=p,
      delta_prior_counts=matrix(rep(2*q/3,each=3),3,17),V=if(t=="ZeroMass") rep(0,3) else c(1,0,1),
      delta_forced=c(FALSE,TRUE,FALSE),input_p=3L,elbo=-12,converged=TRUE,
      mu_grid=array(.1,c(3L,3L,17L)),mu2_grid=array(.2,c(3L,3L,17L)),call=quote(dummy(large_X))),class=c("susie_slide","susie"))
    list(fit_slide_prior=fit,sample_ids=c("A","B"))
  })
  structure(setNames(out,tissues),em_iteration=as.integer(iter))
}
outs<-lapply(seq_along(genes),function(g)make_out(14,priors,g));names(outs)<-genes
for(g in genes) { mb_save(outs[[g]],file.path(seed,"results",paste0(g,".rds")));mb_save(outs[[g]],file.path(root,"iteration_014/results",paste0(g,".rds"))) }
sums<-lapply(genes,function(g)spe_summarize_gene(outs[[g]],g,priors))
counts<-Reduce(`+`,lapply(sums,`[[`,"counts"))
stopifnot(abs(sum(counts[1,])-61*4/3)<1e-12,all(counts[2,]==0))
m<-list(settings=spe_settings(),packages=mb_packages(),genes=genes,datadir="unused",map=map,
  source_iteration=14L,source_project=project,initial_fits=122L,cohorts=lapply(outs,spp_cohort),counts=counts,
  initial_priors=mb_weights(counts,priors))
mb_save(m,file.path(seed,"metadata.rds"));mb_save(priors,file.path(seed,"source_priors.csv"),TRUE)
oldhistory<-do.call(rbind,lapply(0:14,function(i)data.frame(iteration=i,source_iteration=i-1L,priors,
 n_fits=122L,n_active=2L,free_mass=24,prior_max_change=0,m_step_q_gain=0,source_elbo=-12,
 update_method="slider_prior_active_counts_v1",created_at="2026-11-04 UTC")))
mb_save(oldhistory,file.path(seed,"source_prior_history.csv"),TRUE)
mb_save(oldhistory,file.path(root,"prior_history.csv"),TRUE)
original_history_md5<-mb_md5(file.path(root,"prior_history.csv"))
paths<-file.path(seed,"results",paste0(genes,".rds"));info<-file.info(paths)
mb_save(data.frame(gene=genes,size=info$size,mtime=as.numeric(info$mtime),md5=unname(tools::md5sum(paths))),file.path(seed,"expected.csv"),TRUE)
source_hashes<-tools::md5sum(file.path(seed,c("metadata.rds","source_priors.csv","source_prior_history.csv","expected.csv")))
rows<-list();sid<-0L;dir.create(file.path(branch,"plan"));set.seed(31)
for(iteration in 15:38) {
  online<-iteration %in% c(15:20,24:33)
  groups<-if(!online)1L else if(iteration<=16L)9L else if(iteration<=20L)18L else if(iteration==24L)61L else 3L
  rho<-if(!online)1 else if(iteration<=16).3 else if(iteration<=20).5 else if(iteration<=25).8 else if(iteration<=28).5 else .3
  shuffled<-sample(genes);batches<-split(shuffled,rep(seq_len(groups),length.out=length(genes)))
  for(batch in seq_len(groups)) {
    sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(batches[[batch]],file.path(branch,gf))
    group<-sprintf("../../results_slide_prior_em/iteration_%03d",iteration)
    rows[[sid]]<-data.frame(stage=sid,phase=if(online)"minibatch"else"full",pass=iteration-14,batch=batch,rho=rho,
      iteration=iteration,group=group,results_dir=paste0(group,"/results"),
      directory=paste0(group,if(online)sprintf("/batches/batch_%02d",batch)else"/batches/full_step"),genes_file=gf)
  }
}
plan<-do.call(rbind,rows);mb_save(plan,file.path(branch,"stages.csv"),TRUE)
mb_initialize(branch,seed);mb_history(branch)
audit<-normalizePath("recovery/slider_online_extend_024_038_2026-11-09/tests",winslash="/")
python<-"C:/Users/willi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
bridge<-function(mode,index) {
  log<-file.path(branch,paste0("audit_",mode,"_",index,".log"))
  status<-system2(python,shQuote(c(file.path(audit,"bridge.py"),code,branch,mode,as.character(index))),stdout=log,stderr=log)
  if(status!=0)stop(paste(readLines(log),collapse="\n"))
}
bridge("init",1L)
expected_iteration<-setNames(rep(14L,length(genes)),genes);calls<-0L
fake<-function(target_gene,tissue_priors,previous_result,target_gene_id=NULL,...) {
  calls<<-calls+1L;g<-if(is.null(target_gene_id))target_gene else map$gene[map$stable_id==target_gene_id]
  stopifnot(length(g)==1L,attr(previous_result,"em_iteration")==expected_iteration[[g]])
  make_out(0,tissue_priors,match(g,genes))
}
for(i in seq_len(nrow(plan))) {
  bridge("own",i);mb_prepare(branch);mb_prepare(branch);mb_history(branch)
  d<-file.path(branch,plan$directory[i]);x<-readRDS(file.path(d,"input.rds"));before<-readRDS(file.path(branch,"state.rds"))
  stopifnot(all(x$warm$iteration==x$iteration-1L),all(file.exists(x$warm$path)),all(file.exists(x$warm$priors)))
  if(x$batch==1L) {
    mb_save(before$latest,file.path(branch,"audit_warm.csv"),TRUE)
    if(x$iteration>15L) {
      obsolete<-file.path(root,sprintf("iteration_%03d/results",x$iteration-2L))
      saveRDS(list(abandoned=TRUE),tempfile(".em_write_",tmpdir=obsolete))
    }
    bridge(if(x$iteration %in% c(24,26,29,34))"partial"else"cleanup",i)
  }
  if(i==1L)fails(mb_commit(branch,d))
  for(g in x$genes) {
    mb_worker(branch,d,g,run_gene=fake)
    if(x$iteration==23L && g==x$genes[1L]) {
      # Upgrade the actual R engine after 023's input and first gene were saved.
      oldinput<-mb_md5(file.path(d,"input.rds"));oldstate<-mb_md5(file.path(branch,"state.rds"))
      saved<-mb_read(file.path(branch,"prior_history.csv"));oldcalls<-calls
      source(file.path(code,"engine.R"));mb_load(code)
      mb_prepare(branch);mb_worker(branch,d,g,run_gene=function(...)stop("Completed gene must be reused"));mb_history(branch)
      stopifnot(calls==oldcalls,mb_md5(file.path(d,"input.rds"))==oldinput,mb_md5(file.path(branch,"state.rds"))==oldstate,
        isTRUE(all.equal(saved,mb_read(file.path(branch,"prior_history.csv")))))
    }
  }
  raw<-Reduce(`+`,lapply(x$genes,function(g)readRDS(file.path(d,"summaries",paste0(g,".rds")))$counts))
  expected<-(1-x$rho)*before$S+x$rho*raw/length(x$genes)
  statefile<-file.path(branch,"state.rds");oldbytes<-readBin(statefile,"raw",n=file.info(statefile)$size)
  mb_commit(branch,d)
  last<-i==max(which(plan$iteration==x$iteration))
  if(last) { writeBin(oldbytes,statefile);mb_commit(branch,d) }
  mb_commit(branch,d);mb_history(branch);after<-readRDS(statefile)
  stopifnot(after$stage==i,max(abs(after$S-expected))<1e-13,
    max(abs(spe_prior(after$priors,"ZeroMass")-rep(1/17,17)))<1e-13)
  expected_iteration[x$genes]<-x$iteration
  if(last) {
    stopifnot(all(after$latest$iteration==x$iteration),all(file.exists(after$latest$path)),all(file.exists(after$latest$priors)))
    if(x$iteration>=23) {
      z<-mb_read(file.path(root,"iteration_summary.csv"))
      stopifnot(z$status[z$iteration==x$iteration]=="completed",z$updates_completed[z$iteration==x$iteration]==sum(plan$iteration==x$iteration),
        all(z$status[z$iteration>x$iteration]=="planned"))
    }
    cat("PASS scan",x$iteration,": exact counts, prior history, warm starts, cleanup and interrupted commit recovery.\n")
  }
}
stopifnot(calls==61L*24L,after$stage==186L,identical(source_hashes,tools::md5sum(names(source_hashes))),
  mb_md5(file.path(branch,"original_prior_history.csv"))==original_history_md5)
wide<-mb_read(file.path(root,"prior_history.csv"));long<-mb_read(file.path(root,"prior_history_long.csv"))
stopifnot(nrow(long)==nrow(wide)*17L,all(wide$batches_in_iteration[wide$iteration==24L]==61L),
  identical(sort(unique(wide$batch[wide$iteration==24L])),0:61),
  all(wide$batch[wide$phase=="full"]==0L),all(wide$rho[wide$iteration %in% 24:25]==.8),
  all(wide$rho[wide$iteration %in% 26:28]==.5),all(wide$rho[wide$iteration %in% 29:33]==.3))
before<-mb_md5(file.path(root,"prior_history.csv"));mb_history(branch);stopifnot(mb_md5(file.path(root,"prior_history.csv"))==before)
stopifnot(all(vapply(14:36,function(i)length(list.files(file.path(root,sprintf("iteration_%03d/results",i)),pattern="\\.rds$"))==0L,logical(1))),
  all(vapply(37:38,function(i)length(list.files(file.path(root,sprintf("iteration_%03d/results",i)),pattern="\\.rds$"))==61L,logical(1))))
cat("PASS: 186 updates through 038, one actual 61-batch pass, 1464 controlled fits, partial-023 engine adoption, all taper transitions and five final full steps.\n")
