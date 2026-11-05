.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
code<-normalizePath("output/slider_online_014_2026-11-05",winslash="/")
source(file.path(code,"engine.R"));mb_load(code)
project<-tempfile("online_project_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(project)
root<-file.path(project,"results_slide_prior_em");dir.create(root)
branch<-file.path(project,"workflows/slider_online_014_024");dir.create(branch,recursive=TRUE)
seed<-file.path(branch,"seed/snapshot_001");dir.create(file.path(seed,"results"),recursive=TRUE)
fails<-function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
map<-mb_read(file.path(code,"renamed_gene_map.csv"))[1:2,]
genes<-c(map$gene,sprintf("GENE%02d",3:18));tissues<-c("Tissue","ZeroMass")
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
outs<-lapply(seq_along(genes),function(g)make_out(13,priors,g));names(outs)<-genes
for(g in genes)mb_save(outs[[g]],file.path(seed,"results",paste0(g,".rds")))
sums<-lapply(genes,function(g)spe_summarize_gene(outs[[g]],g,priors))
counts<-Reduce(`+`,lapply(sums,`[[`,"counts"))
stopifnot(abs(sum(counts[1,])-18*4/3)<1e-12,all(counts[2,]==0))
m<-list(settings=spe_settings(),packages=mb_packages(),genes=genes,datadir="unused",map=map,
  source_iteration=13L,source_project=project,initial_fits=36L,cohorts=lapply(outs,spp_cohort),counts=counts,
  initial_priors=mb_weights(counts,priors))
mb_save(m,file.path(seed,"metadata.rds"));mb_save(priors,file.path(seed,"source_priors.csv"),TRUE)
oldhistory<-do.call(rbind,lapply(0:13,function(i)data.frame(iteration=i,source_iteration=i-1L,priors,
 n_fits=36L,n_active=2L,free_mass=24,prior_max_change=0,m_step_q_gain=0,source_elbo=-12,
 update_method="slider_prior_active_counts_v1",created_at="2026-11-04 UTC")))
mb_save(oldhistory,file.path(seed,"source_prior_history.csv"),TRUE)
mb_save(oldhistory,file.path(root,"prior_history.csv"),TRUE)
original_history_md5<-mb_md5(file.path(root,"prior_history.csv"))
paths<-file.path(seed,"results",paste0(genes,".rds"));info<-file.info(paths)
mb_save(data.frame(gene=genes,size=info$size,mtime=as.numeric(info$mtime),md5=unname(tools::md5sum(paths))),file.path(seed,"expected.csv"),TRUE)
source_hashes<-tools::md5sum(paths)
rows<-list();sid<-0L;dir.create(file.path(branch,"plan"));set.seed(31)
for(pass in 1:6) {
  shuffled<-sample(genes)
  for(batch in 1:9) {
    sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(shuffled[(batch*2-1):(batch*2)],file.path(branch,gf))
    rows[[sid]]<-data.frame(stage=sid,phase="minibatch",pass=pass,batch=batch,rho=.3,
      iteration=13+pass,group=sprintf("../../results_slide_prior_em/iteration_%03d",13+pass),
      results_dir=sprintf("../../results_slide_prior_em/iteration_%03d/results",13+pass),
      directory=sprintf("../../results_slide_prior_em/iteration_%03d/batches/batch_%02d",13+pass,batch),genes_file=gf)
  }
}
for(step in 1:5) {
  sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(genes,file.path(branch,gf))
  rows[[sid]]<-data.frame(stage=sid,phase="full",pass=step,batch=1,rho=1,
    iteration=19+step,group=sprintf("../../results_slide_prior_em/iteration_%03d",19+step),
    results_dir=sprintf("../../results_slide_prior_em/iteration_%03d/results",19+step),
    directory=sprintf("../../results_slide_prior_em/iteration_%03d/batches/full_step",19+step),genes_file=gf)
}
plan<-do.call(rbind,rows);mb_save(plan,file.path(branch,"stages.csv"),TRUE)
mb_initialize(branch,seed);mb_history(branch)
st<-readRDS(file.path(branch,"state.rds"));stopifnot(isTRUE(all.equal(st$S,counts/18)))
expected_iteration<-setNames(rep(13L,18),genes);calls<-0L
fake<-function(target_gene,tissue_priors,previous_result,target_gene_id=NULL,...) {
  calls<<-calls+1L;g<-if(is.null(target_gene_id)) target_gene else map$gene[map$stable_id==target_gene_id]
  stopifnot(length(g)==1L,attr(previous_result,"em_iteration")==expected_iteration[[g]])
  make_out(0,tissue_priors,match(g,genes))
}
for(i in 1:59) {
  group<-file.path(branch,plan$group[i]);dir.create(group,recursive=TRUE,showWarnings=FALSE)
  writeLines("fixture",file.path(group,"ONLINE_OWNER.json"))
  d<-file.path(branch,plan$directory[i]);mb_prepare(branch);mb_prepare(branch);mb_history(branch)
  x<-readRDS(file.path(d,"input.rds"));before<-readRDS(file.path(branch,"state.rds"))
  if(i==1L) fails(mb_commit(branch,d))
  for(g in x$genes)mb_worker(branch,d,g,run_gene=fake)
  if(i==1L) {
    g<-x$genes[1];outpath<-file.path(x$result_dir,paste0(g,".rds"));saved<-readRDS(outpath)
    stopifnot(is.null(saved[[1]]$fit_slide_prior$call))
    unlink(file.path(d,"success",paste0(g,".csv")))
    mb_worker(branch,d,g,run_gene=function(...)stop("must adopt the already-saved result"))
    bad<-saved;bad[[1]]$fit_slide_prior$converged<-FALSE;mb_save(bad,outpath)
    fails(mb_worker(branch,d,g,run_gene=fake));mb_save(saved,outpath)
    bad<-saved;bad[[1]]$fit_slide_prior$delta_prior<-rev(bad[[1]]$fit_slide_prior$delta_prior);mb_save(bad,outpath)
    fails(mb_worker(branch,d,g,run_gene=fake));mb_save(saved,outpath)
    bad<-saved[-1];attr(bad,"em_iteration")<-14L;mb_save(bad,outpath)
    fails(mb_worker(branch,d,g,run_gene=fake));mb_save(saved,outpath)
    mb_worker(branch,d,g,run_gene=fake)
    sf<-file.path(d,"summaries",paste0(g,".rds"));sm<-readRDS(sf);bad<-sm;bad$counts<-bad$counts*2;mb_save(bad,sf)
    fails(mb_commit(branch,d));mb_save(sm,sf)
  }
  stats<-lapply(x$genes,function(g)readRDS(file.path(d,"summaries",paste0(g,".rds"))))
  raw<-Reduce(`+`,lapply(stats,`[[`,"counts"))
  expected<-if(i<=54L) .7*before$S+.3*raw/length(x$genes) else raw/18
  statefile<-file.path(branch,"state.rds");oldbytes<-readBin(statefile,"raw",n=file.info(statefile)$size)
  mb_commit(branch,d)
  if(i==1L) { writeBin(oldbytes,statefile);mb_commit(branch,d) }
  mb_commit(branch,d);mb_history(branch)
  after<-readRDS(statefile)
  stopifnot(after$stage==i,max(abs(after$S-expected))<1e-13,
    max(abs(spe_prior(after$priors,"ZeroMass")-rep(1/17,17)))<1e-13,
    max(abs(spe_prior(after$priors,"Tissue")-expected[1,]/sum(expected[1,])))<1e-13)
  if(i>54) stopifnot(max(abs(spe_prior(after$priors,"Tissue")-spe_prior(spe_pool_summaries(stats,x$priors)$priors,"Tissue")))<1e-13)
  expected_iteration[x$genes]<-x$iteration
}
stopifnot(calls==18*11,identical(source_hashes,tools::md5sum(paths)),
          nrow(mb_read(file.path(branch,"prior_history.csv")))==(14+11+59)*2,
          length(list.files(file.path(root,"iteration_024/results")))==18L,
          mb_md5(file.path(branch,"original_prior_history.csv"))==original_history_md5)
wide<-mb_read(file.path(root,"prior_history.csv"));long<-mb_read(file.path(root,"prior_history_long.csv"))
stopifnot(identical(names(wide)[seq_along(names(oldhistory))],names(oldhistory)),
 isTRUE(all.equal(wide[wide$iteration<=13,names(oldhistory)],oldhistory,check.attributes=FALSE)),
 nrow(long)==nrow(wide)*17,identical(sort(unique(long$slider_value)),spe_grid()))
for(iter in 14:19)stopifnot(identical(sort(unique(wide$batch[wide$iteration==iter])),0:9))
for(iter in 20:24)stopifnot(identical(sort(unique(wide$batch[wide$iteration==iter])),0L))
before_history<-mb_md5(file.path(root,"prior_history.csv"));mb_history(branch)
stopifnot(mb_md5(file.path(root,"prior_history.csv"))==before_history)
tampered<-wide;tampered$w_01[1]<-0;mb_save(tampered,file.path(root,"prior_history.csv"),TRUE)
fails(mb_history(branch));mb_save(wide,file.path(root,"prior_history.csv"),TRUE)
stopifnot(attr(readRDS(file.path(root,"iteration_024/results",paste0(genes[1],".rds"))),"em_iteration")==24L)
cat("PASS: all 54 minibatch updates + five full EM steps, exact sufficient-statistic math, inactive/forced exclusions,\n",
    "warm starts, renamed IDs, zero-mass fallback, invalid fits rejected, interrupted-output adoption,\n",
    "exactly-once commit recovery, all within-scan weights/grid values, 014..024 numbering, legacy history preserved, source files unchanged.\n")
