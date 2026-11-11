args<-commandArgs(TRUE)
code<-normalizePath(args[1],winslash="/",mustWork=TRUE)
branch<-normalizePath(args[2],winslash="/",mustWork=TRUE)
report<-normalizePath(args[3],winslash="/",mustWork=TRUE)
expected_n<-if(length(args)>3L)as.integer(args[4])else 2029L
source(file.path(code,"engine.R"));mb_load(code)
plan<-mb_read(file.path(branch,"stages.csv"))
d<-normalizePath(file.path(branch,plan$directory[plan$stage==6L]),winslash="/",mustWork=TRUE)
statefile<-file.path(branch,"state.rds");inputfile<-file.path(d,"input.rds")
state<-readRDS(statefile);x<-readRDS(inputfile)
before<-tools::md5sum(c(statefile,inputfile))
stopifnot(state$stage==5L,x$stage==6L,x$iteration==15L,x$batch==6L,x$phase=="minibatch",
          x$rho==.3,length(x$genes)==expected_n,
          identical(x$previous_state_md5,mb_md5(statefile)),
          isTRUE(all.equal(x$S,state$S)),isTRUE(all.equal(x$priors,state$priors)),
          !file.exists(file.path(d,"commit.rds")))
# Execute the real commit calculation with every write redirected into memory.
# No fit, existing receipt, history or state file can be changed by this check.
allowed<-normalizePath(c(file.path(d,c("counts.csv","running_counts_per_gene.rds","prior_after.csv",
    "audit.csv","history.csv","commit.rds")),statefile),winslash="/",mustWork=FALSE)
saved<-new.env(parent=emptyenv())
mb_save<-function(value,path,csv=FALSE) {
  target<-normalizePath(path,winslash="/",mustWork=FALSE)
  if(!target %in% allowed || exists(target,envir=saved,inherits=FALSE))stop("Unexpected preflight write: ",path)
  assign(target,value,envir=saved)
}
mb_worker<-function(...)stop("Fitting forbidden in commit preflight")
run_slide_prior_gene<-function(...)stop("Fitting forbidden in commit preflight")
cat("Dry-run commit: validating summaries and computing the update in memory only.\n")
preview<-mb_commit(branch,d)
stopifnot(length(ls(saved))==7L,preview$stage==6L,preview$state$stage==6L,
          identical(sort(preview$genes),sort(x$genes)),
          identical(preview$state$latest$gene,state$latest$gene),
          identical(before,tools::md5sum(c(statefile,inputfile))),
          !file.exists(file.path(d,"commit.rds")))
write.csv(preview$history,file.path(report,"preview_history.csv"),row.names=FALSE)
write.csv(preview$state$priors,file.path(report,"preview_prior_after.csv"),row.names=FALSE)
cat("PASS:",length(x$genes),"summaries validated; batch 6 dry-run commit succeeded.\n")
cat("The on-disk state remains at stage 5; no completed gene was refitted.\n")
