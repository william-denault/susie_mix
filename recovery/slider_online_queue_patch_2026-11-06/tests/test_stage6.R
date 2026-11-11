.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
args<-commandArgs(TRUE)
patch<-normalizePath(if(length(args))args[1]else"recovery/slider_online_queue_patch_2026-11-06",winslash="/")
code<-file.path(patch,"runtime")
fixture<-paste(readLines(file.path(code,"tests/test_engine.R")),collapse="\n")
replace_once<-function(text,old,new) {
  hits<-gregexpr(old,text,fixed=TRUE)[[1]];stopifnot(length(hits)==1L,hits[1]>0)
  sub(old,new,text,fixed=TRUE)
}
fixture<-replace_once(fixture,'code<-normalizePath("output/slider_online_015_2026-11-05",winslash="/")',
                     paste0('code<-',encodeString(code,quote='"')))
fixture<-replace_once(fixture,'genes<-c(map$gene,sprintf("GENE%02d",3:18))',
                     'genes<-c(map$gene,"T","F",sprintf("GENE%02d",5:18))')
hook<-paste0('
  if(i==6L) {
    protected<-unique(c(file.path(branch,"state.rds"),
      list.files(group,recursive=TRUE,full.names=TRUE)))
    hashes<-tools::md5sum(protected)
    report<-file.path(branch,"queue_repair_preview");dir.create(report)
    log<-file.path(report,"check.log")
    status<-system2(file.path(R.home("bin"),"Rscript"),
      shQuote(c("--vanilla",',encodeString(file.path(patch,"check_commit.R"),quote='"'),',code,branch,report,"2")),stdout=log,stderr=log)
    if(status!=0)stop(paste(readLines(log),collapse="\\n"))
    stopifnot(identical(hashes,tools::md5sum(protected)),
              readRDS(file.path(branch,"state.rds"))$stage==5L,!file.exists(file.path(d,"commit.rds")))
    cat("PASS: production batch-6 dry-run preserved every saved fit, summary, receipt and prior update.\\n")
  }
')
fixture<-replace_once(fixture,'  mb_commit(branch,d)\n',paste0(hook,'  mb_commit(branch,d)\n'))
eval(parse(text=fixture),envir=globalenv())
