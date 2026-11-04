# Finalize only actual modifications. Untouched files use the same size/mtime
# validity rule as the existing certified parallel controller, not a new scan.
source_helpers <- function(bundle) {
  source(file.path(bundle,"maintenance.R"),local=globalenv());load_helpers(bundle)
}
finalize_partial <- function(bundle,project,state) {
  source_helpers(bundle)
  s<-readRDS(file.path(state,"snapshot.rds"));assert_history(s)
  if(normalizePath(project,winslash="/")!=s$project) stop("Wrong project")
  ledger<-read.csv(file.path(state,"receipts.csv"),stringsAsFactors=FALSE,colClasses="character")
  if(anyDuplicated(ledger$relative)) stop("Duplicate receipt")
  report<-list()
  for(i in 0:s$latest) {
    ctx<-sre_context(project,file.path(project,"results_slide_prior_em",sprintf("iteration_%03d",i)),check_package=TRUE)
    old<-s$iterations[[as.character(i)]]
    if(isTRUE(old$retired)) next
    genes<-ctx$manifest$gene
    # Inspect all small saved records/stamps before modifying this iteration.
    checks<-setNames(vector("list",length(genes)),genes)
    for(j in seq_along(genes)) {
      g<-genes[j];if(g %in% ctx$scope$excluded) next
      x<-old$summaries[[g]]
      before<-old$result_stamps[j,,drop=FALSE];rownames(before)<-NULL
      prev<-if(i) old$previous_stamps[j,,drop=FALSE] else spp_stamp(character());rownames(prev)<-NULL
      good<-spp_valid_summary(x,ctx,g,before,prev)
      rel<-paste0("results_slide_prior_em/",basename(ctx$dir),"/results/",g,".rds")
      receipt<-ledger[ledger$relative==rel,,drop=FALSE]
      exists<-file.exists(sre_result(ctx,g))
      if(!nrow(receipt)) {
        if(!identical(before,spp_stamp(sre_result(ctx,g)))) stop("Unjournaled result changed: ",rel)
      } else {
        if(!exists) stop("Committed result disappeared: ",rel)
        if(good && x$row$result_md5!=receipt$old_md5) stop("Original cached hash differs from transaction: ",rel)
        if(good && x$row$status=="reviewed_exclusion" && receipt$old_md5!=receipt$new_md5 && g %in% sre_exclusions(ctx)$gene)
          stop("Reviewed exclusion needs separate migration: ",g)
      }
      if(!is.null(old$cert) && (!good || !exists || !sre_ready(x$row$status))) stop("Completed checkpoint has invalid saved summary: ",g)
      checks[[g]]<-list(good=good,exists=exists,receipt=receipt)
    }
    deferred<-0L
    for(g in setdiff(genes,ctx$scope$excluded)) {
      check<-checks[[g]];x<-old$summaries[[g]]
      if(check$good && check$exists) {
        x$result_stamp<-spp_stamp(sre_result(ctx,g));x$previous_stamp<-spp_stamp(spp_previous_file(ctx,g))
        if(nrow(check$receipt)) x$row$result_md5<-check$receipt$new_md5
        if(!identical(sre_read(spp_summary_file(ctx,g)),x)) save_checked(x,spp_summary_file(ctx,g))
      } else {
        if(check$exists) deferred<-deferred+1L
        # Uncached 005 results are validated by the existing parallel workers.
        # Do not deserialize hundreds of huge fits serially on this coordinator.
        stale<-spp_summary_file(ctx,g)
        if(file.exists(stale)) {
          archive<-file.path(state,"deferred_summaries",basename(ctx$dir),paste0(g,".rds"))
          dir.create(dirname(archive),recursive=TRUE,showWarnings=FALSE)
          if(!file.exists(archive) && !file.copy(stale,archive)) stop("Cannot preserve old summary")
          if(unlink(stale)!=0L) stop("Cannot invalidate stale summary")
        }
      }
    }
    index<-spp_index(ctx)
    ready<-vapply(index$entries,function(x) !is.null(x) && sre_ready(x$row$status),logical(1))
    if(!is.null(old$cert)) {
      if(!all(ready)) stop("Completed checkpoint lost validation")
      if(old$pooled_md5!=md5(file.path(spp_root(ctx),"pooled.rds"))) stop("Pooled counts changed")
      cert<-old$cert
      cert$result_stamps<-spp_stamp(sre_result(ctx,genes))
      cert$previous_stamps<-if(i) spp_stamp(vapply(genes,function(g) spp_previous_file(ctx,g),character(1))) else NULL
      cert$summary_stamps<-spp_stamp(vapply(genes,function(g) spp_summary_file(ctx,g),character(1)))
      save_checked(cert,file.path(spp_root(ctx),"COMPLETE.rds"));spp_require_complete(project,ctx$dir)
    }
    audit_file<-file.path(ctx$root,"audit.csv")
    if(file.exists(audit_file)) {
      a<-read.csv(audit_file,stringsAsFactors=FALSE)
      for(g in intersect(a$gene,genes)) {
        x<-index$entries[[g]]
        if(!is.null(x)) a$result_md5[a$gene==g]<-x$row$result_md5
      }
      em_atomic_write(a,audit_file,csv=TRUE)
    }
    if(i==s$latest) {
      archive<-file.path(state,"quota_attempt_receipts.rds")
      reviewed<-if(file.exists(archive)) readRDS(archive) else list()
      for(g in genes[!ready]) {
        a<-sre_read(sre_attempt_file(ctx,g))
        if(!is.null(a$error) && grepl("error writing|Disk quota|No space",a$error,ignore.case=TRUE)) {
          reviewed[[g]]<-a;save_checked(reviewed,archive);unlink(sre_attempt_file(ctx,g))
        }
      }
      batches<-list.files(file.path(spp_root(ctx),"batches"),"^array_job_id\\.txt$",recursive=TRUE,full.names=TRUE)
      for(f in batches) { spp_batch(ctx,dirname(f));save_checked(TRUE,file.path(dirname(f),"closed.rds")) }
      save_checked(list(reviewed_batches=dirname(batches),reason="Partial storage finalization; no live slider jobs"),file.path(state,"closed_batches.rds"))
    }
    report[[as.character(i)]]<-data.frame(iteration=i,existing_results=sum(file.exists(sre_result(ctx,genes))),
      validated_genes=sum(ready & !genes %in% ctx$scope$excluded),scope_excluded=length(ctx$scope$excluded),
      pending=sum(!ready),deferred_file_validation=deferred,complete=!is.null(old$cert))
    cat("Repaired iteration",i,"; validated genes:",sum(ready & !genes %in% ctx$scope$excluded),
        "; saved files deferred to parallel validation:",deferred,"\n")
  }
  assert_history(s)
  write.csv(do.call(rbind,report),file.path(state,"validation_report.csv"),row.names=FALSE)
  cat("Partial compaction finalized. All priors/history/pooled counts preserved; no fitting and no bulk recompression.\n")
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE)
  if(length(a)!=3L) stop("Usage: partial_finalize.R BUNDLE PROJECT STATE")
  finalize_partial(a[1],a[2],a[3])
}
