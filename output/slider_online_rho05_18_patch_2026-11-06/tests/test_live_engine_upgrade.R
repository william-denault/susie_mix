.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
args<-commandArgs(TRUE)
inject_orphans<-TRUE
code<-normalizePath(if(length(args))args[1]else"recovery/slider_online_rho05_18_patch_2026-11-06/runtime",winslash="/")
audit<-normalizePath("recovery/slider_online_rho05_18_patch_2026-11-06/tests",winslash="/")
python<-"C:/Users/willi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
Sys.setenv(RSCRIPT=file.path(R.home("bin"),"Rscript"),R_LIBS_USER=paste(.libPaths(),collapse=.Platform$path.sep))
fixture<-paste(readLines(file.path(code,"tests/test_engine.R")),collapse="\n")
replace_once<-function(text,old,new) {
  hits<-gregexpr(old,text,fixed=TRUE)[[1]];stopifnot(length(hits)==1L,hits[1]>0)
  sub(old,new,text,fixed=TRUE)
}
fixture<-replace_once(fixture,'code<-normalizePath("output/slider_online_015_2026-11-05",winslash="/")',
                     paste0('code<-',encodeString(code,quote='"')))
fixture<-replace_once(fixture,'tempfile("online_project_"','tempfile("transition_integrated_"')
fixture<-replace_once(fixture,'genes<-c(map$gene,sprintf("GENE%02d",3:18))',
                     'genes<-c(map$gene,"T","F",sprintf("GENE%02d",5:18))')
fixture<-replace_once(fixture,'for(g in genes)mb_save(outs[[g]],file.path(seed,"results",paste0(g,".rds")))',
 'for(g in genes) { mb_save(outs[[g]],file.path(seed,"results",paste0(g,".rds")));mb_save(outs[[g]],file.path(root,"iteration_014/results",paste0(g,".rds"))) }')
fixture<-replace_once(fixture,'source_hashes<-tools::md5sum(paths)',
 'source_hashes<-tools::md5sum(file.path(seed,c("metadata.rds","source_priors.csv","source_prior_history.csv","expected.csv")))')
fixture<-replace_once(fixture,'identical(source_hashes,tools::md5sum(paths))',
 'identical(source_hashes,tools::md5sum(names(source_hashes)))')
setup<-paste0('
bridge<-function(mode,index) {
  log<-file.path(branch,paste0("audit_",mode,"_",index,".log"))
  status<-system2(',encodeString(python,quote='"'),',shQuote(c(',encodeString(file.path(audit,"bridge.py"),quote='"'),',code,branch,mode,as.character(index))),stdout=log,stderr=log)
  if(status!=0)stop(paste(readLines(log),collapse="\\n"))
  cat(paste(readLines(log),collapse="\\n"),"\\n")
}
bridge("init",1L)
# Seed helper must use the real deployed R engine while queue access is mocked.
dir.create(file.path(branch,"code"));file.copy(file.path(code,"engine.R"),file.path(branch,"code/engine.R"))
dir.create(file.path(branch,"code/script"));file.copy(file.path(code,"script/scan_tissue_attempt"),file.path(branch,"code/script"),recursive=TRUE)
')
fixture<-replace_once(fixture,'expected_iteration<-setNames(rep(14L,18),genes);calls<-0L',
 paste0(setup,'\nexpected_iteration<-setNames(rep(14L,18),genes);calls<-0L'))
fixture<-replace_once(fixture,'  group<-file.path(branch,plan$group[i]);dir.create(group,recursive=TRUE,showWarnings=FALSE)',
 '  group<-file.path(branch,plan$group[i])')
fixture<-replace_once(fixture,'  writeLines("fixture",file.path(group,"ONLINE_OWNER.json"))','  bridge("own",i)')
hook<-'
  # Every warm start must come from the previous scan and its own saved batch prior.
  stopifnot(all(x$warm$iteration==x$iteration-1L),all(file.exists(x$warm$path)),all(file.exists(x$warm$priors)))
  for(j in seq_len(nrow(x$warm))) {
    w<-x$warm[j,];z<-readRDS(w$path)
    stopifnot(identical(mb_md5(w$path),w$md5))
    sre_validate(z,list(iteration=w$iteration,priors=mb_read(w$priors),settings=spe_settings()),x$cohorts[[w$gene]])
  }
  mb_save(before$latest,file.path(branch,"audit_warm.csv"),TRUE)
  if(inject_orphans && x$batch==1L && x$iteration>15L) {
    obsolete<-file.path(root,sprintf("iteration_%03d/results",x$iteration-2L))
    stopifnot(dir.exists(obsolete))
    saveRDS(list(incomplete_atomic_save_fixture=TRUE),tempfile(".em_write_",tmpdir=obsolete))
  }
  # Exercise actual deletion before workers, including interruption at 016 and 017.
  bridge(if(i %in% c(10L,19L))"partial"else"cleanup",i)
  if(i %in% c(10L,19L,55L,59L))bridge("seed",i)
  stopifnot(all(file.exists(x$warm$path)),all(file.exists(x$warm$priors)))
'
fixture<-replace_once(fixture,'  if(i==1L) fails(mb_commit(branch,d))',paste0(hook,'\n  if(i==1L) fails(mb_commit(branch,d))'))
# Controller termination between durable commit and state update at the end of each online pass.
fixture<-replace_once(fixture,'  if(i==1L) { writeBin(oldbytes,statefile);mb_commit(branch,d) }',
 '  if(i %in% c(1L,9L,18L,27L,36L,45L,54L,55L,59L)) { writeBin(oldbytes,statefile);mb_commit(branch,d) }')

fixture<-replace_once(fixture,'  for(batch in 1:9) {','  groups<-if(pass<=2L)9L else 18L; sizes<-rep(length(genes)%/%groups,groups);sizes[seq_len(length(genes)%%groups)]<-sizes[seq_len(length(genes)%%groups)]+1L;ends<-cumsum(sizes)\n  for(batch in seq_len(groups)) {')
fixture<-replace_once(fixture,'shuffled[(batch*2-1):(batch*2)]','shuffled[seq.int(ends[batch]-sizes[batch]+1L,ends[batch])]')
fixture<-replace_once(fixture,'pass=pass,batch=batch,rho=.3,','pass=pass,batch=batch,rho=if(pass<=2L).3 else .5,')
fixture<-replace_once(fixture,'for(i in 1:59) {','for(i in seq_len(nrow(plan))) {')
fixture<-replace_once(fixture,'if(i<=54L) .7*before$S+.3*raw/length(x$genes) else raw/18','(1-x$rho)*before$S+x$rho*raw/length(x$genes)')
fixture<-replace_once(fixture,'if(i>54)','if(x$phase=="full")')
fixture<-replace_once(fixture,'(15+11+59)*2','(15+11+95)*2')
fixture<-replace_once(fixture,'c(1L,9L,18L,27L,36L,45L,54L,55L,59L)','c(1L,9L,18L,36L,54L,71L,72L,90L,91L,95L)')
fixture<-replace_once(fixture,'for(iter in 15:20)stopifnot(identical(sort(unique(wide$batch[wide$iteration==iter])),0:9))',
 'for(iter in 15:20) { n<-if(iter<=16L)9L else 18L; stopifnot(identical(sort(unique(wide$batch[wide$iteration==iter])),0:n),all(wide$batches_in_iteration[wide$iteration==iter]==n)) }')
fixture<-replace_once(fixture,'all 54 minibatch updates + five full EM steps','all 90 mixed-schedule minibatch updates + five full EM steps')


# Start with the deployed engine and swap at a partially fitted 016 batch 7.
baseline<-normalizePath("output/slider_online_rho05_12_patch_2026-11-06/runtime",winslash="/")
fixture<-replace_once(fixture,'source(file.path(code,"engine.R"));mb_load(code)',paste0('source(',encodeString(file.path(baseline,"engine.R"),quote='"'),');mb_load(code)'))
fixture<-replace_once(fixture,'for(i in seq_len(nrow(plan))) {','for(i in 1:19) {')
fixture<-strsplit(fixture,'stopifnot(calls==18*11',fixed=TRUE)[[1]][1]
hook<-'
  if(i==16L) {
    saved_input<-mb_md5(file.path(d,"input.rds"));saved_state<-mb_md5(file.path(branch,"state.rds"))
    g<-readRDS(file.path(d,"input.rds"))$genes[1L]
    mb_worker(branch,d,g,run_gene=fake)
    saved_result<-file.path(readRDS(file.path(d,"input.rds"))$result_dir,paste0(g,".rds"));saved_md5<-mb_md5(saved_result)
    source(file.path(code,"engine.R"));mb_load(code)
    mb_prepare(branch);mb_history(branch)
    stopifnot(mb_md5(file.path(d,"input.rds"))==saved_input,mb_md5(file.path(branch,"state.rds"))==saved_state)
    mb_worker(branch,d,g,run_gene=function(...)stop("Completed gene must not refit after code upgrade"))
    stopifnot(mb_md5(saved_result)==saved_md5)
  }
'
fixture<-replace_once(fixture,'  x<-readRDS(file.path(d,"input.rds"));before<-readRDS(file.path(branch,"state.rds"))',paste0(hook,'\n  x<-readRDS(file.path(d,"input.rds"));before<-readRDS(file.path(branch,"state.rds"))'))
eval(parse(text=fixture),envir=globalenv())
h<-mb_read(file.path(root,"prior_history.csv"))
stopifnot(after$stage==19L,calls==37L,all(h$batches_in_iteration[h$iteration==16]==9L),all(h$batches_in_iteration[h$iteration==17]==18L))
stopifnot(tail(h$rho,1)==.5)
cat("PASS: engine upgrade in partially fitted 016 batch 7 preserves input, state and completed fit; publishes mixed history and starts 017 with rho=.5.\n")
