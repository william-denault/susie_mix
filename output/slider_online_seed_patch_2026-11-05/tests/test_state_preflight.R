# Build the existing synthetic 59-stage fixture only up to stage-zero setup.
lines<-readLines("output/slider_online_015_2026-11-05/tests/test_engine.R")
cut<-grep("^for\\(i in 1:59\\)",lines)
stopifnot(length(cut)==1L)
lines<-lines[seq_len(cut-1L)]
lines<-gsub("output/slider_online_015_2026-11-05","recovery/slider_online_seed_patch_2026-11-05/runtime",lines,fixed=TRUE)
eval(parse(text=lines),envir=globalenv())
mb_prepare(branch);mb_history(branch)
check<-normalizePath("recovery/slider_online_seed_patch_2026-11-05/check_state.R",winslash="/")
run_check<-function() {
  output<-system2(Sys.getenv("RSCRIPT"),c("--vanilla",shQuote(check),shQuote(branch),"18"),stdout=TRUE,stderr=TRUE)
  status<-attr(output,"status");if(is.null(status))status<-0L
  list(status=status,output=output)
}
good<-run_check();stopifnot(good$status==0L)
# A different stage must not qualify for a repair that resets failed attempts.
statefile<-file.path(branch,"state.rds");original<-readBin(statefile,"raw",n=file.info(statefile)$size)
changed<-readRDS(statefile);changed$stage<-1L;mb_save(changed,statefile)
bad<-suppressWarnings(run_check());stopifnot(bad$status!=0L)
writeBin(original,statefile)
cat("PASS: real R seed/state/input/prior validation accepts stage 0 and refuses committed updates.\n")
