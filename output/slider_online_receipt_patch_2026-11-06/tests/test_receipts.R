args<-commandArgs(TRUE)
patch<-normalizePath(if(length(args)) args[1] else "recovery/slider_online_receipt_patch_2026-11-06",winslash="/")
code<-file.path(patch,"runtime")
source(file.path(code,"engine.R"));mb_load(code)
scratch<-tempfile("literal_receipts_",tmpdir=normalizePath("tmp",winslash="/"));dir.create(scratch)
for(g in c("T","F","TRUE","FALSE","NA","0001","9E10","ABCA1")) {
  row<-data.frame(gene=g,stage=3L,input_md5="00000000000000000000000000000001",
    result_md5="12345678901234567890123456789012",result_path=paste0("/saved/",g,".rds"),
    summary_md5="abcdef1234567890abcdef1234567890",n_fits=27L)
  f<-file.path(scratch,paste0(g,".csv"));write.csv(row,f,row.names=FALSE)
  original<-tools::md5sum(f);fixed<-mb_read_success(f)
  stopifnot(identical(fixed,row),identical(original,tools::md5sum(f)))
}
stopifnot(identical(mb_read(file.path(scratch,"T.csv"))$gene,TRUE))
cat("PASS: literal T/F/NA and numeric-looking identifiers remain text; CSV files unchanged.\n")

# Reuse the existing full-schedule fixture, with T first appearing in batch 3.
# Both full EM and online updates are checked against independently pooled counts.
baseline<-normalizePath("output/slider_online_submit_patch_2026-11-05/runtime",winslash="/")
fixture<-paste(readLines(file.path(baseline,"tests/test_engine.R")),collapse="\n")
replace_once<-function(text,old,new) {
  hits<-gregexpr(old,text,fixed=TRUE)[[1]];stopifnot(length(hits)==1L,hits[1]>0)
  sub(old,new,text,fixed=TRUE)
}
fixture<-replace_once(fixture,'code<-normalizePath("output/slider_online_015_2026-11-05",winslash="/")',
                     paste0('code<-',encodeString(code,quote='"')))
fixture<-replace_once(fixture,'genes<-c(map$gene,sprintf("GENE%02d",3:18))',
                     'genes<-c(map$gene,"T","F",sprintf("GENE%02d",5:18))')
fixture<-replace_once(fixture,'  shuffled<-sample(genes)',
  '  shuffled<-sample(genes)\n  if(pass==1L) { other<-setdiff(shuffled,c("T","F"));shuffled<-c(other[1:4],"T",other[5],"F",other[-(1:5)]) }')
hook<-paste0('
  if(i==3L) {
    oldenv<-new.env(parent=globalenv());sys.source(',encodeString(file.path(baseline,"engine.R"),quote='"'),',oldenv)
    olderror<-tryCatch(oldenv$mb_commit(branch,d),error=function(e)e)
    stopifnot(inherits(olderror,"error"),conditionMessage(olderror)=="Missing/changed validated summary: T")
    fits<-list.files(file.path(group,"results"),full.names=TRUE)
    unchanged<-c(fits,file.path(branch,"state.rds"),list.files(file.path(d,"success"),full.names=TRUE))
    saved_hashes<-tools::md5sum(unchanged)
    report<-file.path(branch,"preflight_preview");dir.create(report)
    log<-file.path(report,"check.log")
    status<-system2(file.path(R.home("bin"),"Rscript"),
      shQuote(c("--vanilla",',encodeString(file.path(patch,"check_commit.R"),quote='"'),',code,branch,report,"2")),stdout=log,stderr=log)
    if(status!=0)stop(paste(readLines(log),collapse="\\n"))
    stopifnot(identical(saved_hashes,tools::md5sum(unchanged)),
              readRDS(file.path(branch,"state.rds"))$stage==2L,!file.exists(file.path(d,"commit.rds")))
    cat("PASS: original T commit failure reproduced; production preflight passed without modifying fits, receipts or state.\\n")
  }
')
fixture<-replace_once(fixture,'  mb_commit(branch,d)\n',paste0(hook,'  mb_commit(branch,d)\n'))
eval(parse(text=fixture),envir=globalenv())
