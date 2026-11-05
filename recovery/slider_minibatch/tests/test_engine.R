.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
code<-normalizePath("output/slider_minibatch_2026-11-04",winslash="/")
source(file.path(code,"engine.R"));mb_load(code)
branch<-tempfile("mb_engine_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(branch)
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
outs<-lapply(seq_along(genes),function(g)make_out(12,priors,g));names(outs)<-genes
for(g in genes)mb_save(outs[[g]],file.path(seed,"results",paste0(g,".rds")))
sums<-lapply(genes,function(g)spe_summarize_gene(outs[[g]],g,priors))
counts<-Reduce(`+`,lapply(sums,`[[`,"counts"))
stopifnot(abs(sum(counts[1,])-18*4/3)<1e-12,all(counts[2,]==0))
m<-list(settings=spe_settings(),packages=mb_packages(),genes=genes,datadir="unused",map=map,
  source_iteration=12L,initial_fits=36L,cohorts=lapply(outs,spp_cohort),counts=counts,
  initial_priors=mb_weights(counts,priors))
mb_save(m,file.path(seed,"metadata.rds"));mb_save(priors,file.path(seed,"source_priors.csv"),TRUE)
paths<-file.path(seed,"results",paste0(genes,".rds"));info<-file.info(paths)
mb_save(data.frame(gene=genes,size=info$size,mtime=as.numeric(info$mtime),md5=unname(tools::md5sum(paths))),file.path(seed,"expected.csv"),TRUE)
source_hashes<-tools::md5sum(paths)
rows<-list();sid<-0L;dir.create(file.path(branch,"plan"));set.seed(31)
for(pass in 1:6) {
  shuffled<-sample(genes)
  for(batch in 1:9) {
    sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(shuffled[(batch*2-1):(batch*2)],file.path(branch,gf))
    rows[[sid]]<-data.frame(stage=sid,phase="minibatch",pass=pass,batch=batch,rho=.3,
      directory=sprintf("passes/pass_%03d/batch_%02d",pass,batch),genes_file=gf)
  }
}
for(step in 1:5) {
  sid<-sid+1L;gf<-sprintf("plan/%03d.txt",sid);writeLines(genes,file.path(branch,gf))
  rows[[sid]]<-data.frame(stage=sid,phase="full",pass=step,batch=0,rho=1,
    directory=sprintf("full_em/iteration_%03d",step),genes_file=gf)
}
plan<-do.call(rbind,rows);mb_save(plan,file.path(branch,"stages.csv"),TRUE)
mb_initialize(branch,seed);mb_history(branch)
st<-readRDS(file.path(branch,"state.rds"));stopifnot(isTRUE(all.equal(st$S,counts/18)))
expected_iteration<-setNames(rep(12L,18),genes);calls<-0L
fake<-function(target_gene,tissue_priors,previous_result,target_gene_id=NULL,...) {
  calls<<-calls+1L;g<-if(is.null(target_gene_id)) target_gene else map$gene[map$stable_id==target_gene_id]
  stopifnot(length(g)==1L,attr(previous_result,"em_iteration")==expected_iteration[[g]])
  make_out(0,tissue_priors,match(g,genes))
}
for(i in 1:59) {
  d<-file.path(branch,plan$directory[i]);mb_prepare(branch);mb_prepare(branch)
  x<-readRDS(file.path(d,"input.rds"));before<-readRDS(file.path(branch,"state.rds"))
  if(i==1L) fails(mb_commit(branch,d))
  for(g in x$genes)mb_worker(branch,d,g,run_gene=fake)
  if(i==1L) {
    g<-x$genes[1];outpath<-file.path(d,"results",paste0(g,".rds"));saved<-readRDS(outpath)
    stopifnot(is.null(saved[[1]]$fit_slide_prior$call))
    unlink(file.path(d,"success",paste0(g,".csv")))
    mb_worker(branch,d,g,run_gene=function(...)stop("must adopt the already-saved result"))
    bad<-saved;bad[[1]]$fit_slide_prior$converged<-FALSE;mb_save(bad,outpath)
    fails(mb_worker(branch,d,g,run_gene=fake));mb_save(saved,outpath)
    bad<-saved;bad[[1]]$fit_slide_prior$delta_prior<-rev(bad[[1]]$fit_slide_prior$delta_prior);mb_save(bad,outpath)
    fails(mb_worker(branch,d,g,run_gene=fake));mb_save(saved,outpath)
    bad<-saved[-1];attr(bad,"em_iteration")<-1L;mb_save(bad,outpath)
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
  expected_iteration[x$genes]<-i
}
stopifnot(calls==18*11,identical(source_hashes,tools::md5sum(paths)),
          nrow(mb_read(file.path(branch,"prior_history.csv")))==60*2,
          length(list.files(file.path(branch,"full_em/iteration_005/results")))==18L)
cat("PASS: all 54 minibatch updates + five full EM steps, exact sufficient-statistic math, inactive/forced exclusions,\n",
    "warm starts, renamed IDs, zero-mass fallback, invalid fits rejected, interrupted-output adoption,\n",
    "exactly-once commit recovery, full prior history, source files unchanged.\n")
