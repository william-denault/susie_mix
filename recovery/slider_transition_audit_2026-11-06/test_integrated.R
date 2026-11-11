.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
args<-commandArgs(TRUE)
inject_orphans<-length(args)>1L && args[2]=="orphans"
code<-normalizePath(if(length(args))args[1]else"output/slider_online_queue_patch_2026-11-06/runtime",winslash="/")
audit<-normalizePath("recovery/slider_transition_audit_2026-11-06",winslash="/")
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
eval(parse(text=fixture),envir=globalenv())
stopifnot(length(list.files(file.path(root,"iteration_024/results"),pattern="\\.rds$"))==18L,
          length(list.files(file.path(root,"iteration_025/results"),pattern="\\.rds$"))==18L,
          all(vapply(14:23,function(i)length(list.files(file.path(root,sprintf("iteration_%03d/results",i)),pattern="\\.rds$"))==0L,logical(1))))
cat("PASS: integrated 015 through 025; deleted older result files, warm-started every gene, preserved every history/count, recovered interrupted cleanup and pass-boundary commits.\n")
