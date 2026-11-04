# No model fitting: remove only fit$call, then repair content hashes/stamps.
load_helpers <- function(bundle) {
  for (f in c("em_utils.R", "slide_prior_em_utils.R", "slide_prior_recovery.R", "slide_prior_parallel.R"))
    sys.source(file.path(bundle, "baseline/script/scan_tissue_attempt", f), globalenv())
}
md5 <- function(p) unname(tools::md5sum(p))
save_checked <- function(x, p) {
  dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)
  temp <- paste0(p, ".storage-tmp")
  saveRDS(x, temp)
  stopifnot(identical(readRDS(temp), x))
  if (.Platform$OS.type == "windows" && file.exists(p)) {
    if (!file.copy(temp,p,overwrite=TRUE)) stop("Cannot replace ",p)
    unlink(temp)
  } else if (!file.rename(temp,p)) stop("Cannot replace ",p)
}
strip_calls <- function(out) {
  if (!is.list(out)) stop("Gene result is not a list")
  n <- 0L
  if (is.null(out$error)) for (t in names(out)) {
    if (!is.list(out[[t]]) || !inherits(out[[t]]$fit_slide_prior, "susie_slide"))
      stop("Unrecognized tissue result: ", t)
    if (!is.null(out[[t]]$fit_slide_prior$call)) {
      out[[t]]$fit_slide_prior$call <- NULL
      n <- n + 1L
    }
  }
  list(out=out, removed=n)
}
compact_file <- function(input, output, report) {
  x <- strip_calls(readRDS(input))
  saveRDS(x$out, output, compress="gzip")
  if (!identical(readRDS(output), x$out)) stop("Round-trip changed scientific fields")
  write.csv(data.frame(removed_calls=x$removed), report, row.names=FALSE)
}
history_files <- function(project) {
  root <- file.path(project,"results_slide_prior_em")
  dirs <- list.dirs(root,recursive=FALSE,full.names=TRUE)
  dirs <- dirs[grepl("^iteration_[0-9]+$",basename(dirs))]
  files <- c(file.path(root,"prior_history.csv"), unlist(lapply(dirs, function(d)
    file.path(d,c("prior_history.csv","priors.csv","grid.csv","settings.rds","manifest.csv",
                  "source_audit.csv","pooled_counts.csv","slider_counts.csv","proposed_priors.csv")))))
  files[file.exists(files)]
}
snapshot <- function(bundle,project,state,expected) {
  load_helpers(bundle)
  root <- file.path(project,"results_slide_prior_em")
  history <- read.csv(file.path(root,"prior_history.csv"))
  latest <- max(history$iteration)
  if (latest != expected) stop("Latest iteration is ",latest,"; expected ",expected)
  result <- list(project=normalizePath(project,winslash="/"), latest=latest,
                 preserved=tools::md5sum(history_files(project)), iterations=list())
  for (i in 0:latest) {
    ctx <- sre_context(project,file.path(root,sprintf("iteration_%03d",i)),check_package=TRUE)
    if (file.exists(file.path(ctx$dir,"storage_retired.json"))) {
      latest_cert <- file.path(root,sprintf("iteration_%03d",latest),"recovery/parallel/COMPLETE.rds")
      if(i >= latest || (i == latest-1L && !file.exists(latest_cert))) stop("Required warm-start checkpoint has been retired")
      result$iterations[[as.character(i)]] <- list(retired=TRUE)
      next
    }
    cert <- sre_read(file.path(spp_root(ctx),"COMPLETE.rds"))
    if (i < latest || !is.null(cert)) spp_require_complete(project,ctx$dir)
    genes <- ctx$manifest$gene
    entries <- setNames(lapply(genes,function(g) sre_read(spp_summary_file(ctx,g))),genes)
    result$iterations[[as.character(i)]] <- list(cert=cert, summaries=entries,
      result_stamps=spp_stamp(sre_result(ctx,genes)),
      previous_stamps=if(i) spp_stamp(vapply(genes,function(g) spp_previous_file(ctx,g),character(1))) else NULL,
      pooled_md5=if(!is.null(cert)) md5(file.path(spp_root(ctx),"pooled.rds")) else NULL)
    cat("Snapshot iteration",i,"; complete:",!is.null(cert),"\n")
  }
  save_checked(result,file.path(state,"snapshot.rds"))
  file.copy(file.path(root,"prior_history.csv"),file.path(state,"prior_history.preserved.csv"),overwrite=TRUE)
}
assert_history <- function(s) {
  if (!identical(s$preserved,tools::md5sum(names(s$preserved))))
    stop("A frozen input or prior-history file changed; migration stopped")
}
finalize <- function(bundle,project,state) {
  load_helpers(bundle)
  s <- readRDS(file.path(state,"snapshot.rds")); assert_history(s)
  if (!identical(normalizePath(project,winslash="/"),s$project)) stop("Wrong project")
  ledger <- read.csv(file.path(state,"receipts.csv"),stringsAsFactors=FALSE,colClasses="character")
  if (anyDuplicated(ledger$relative)) stop("Duplicate receipt")
  report <- list()
  for (i in 0:s$latest) {
    ctx <- sre_context(project,file.path(project,"results_slide_prior_em",sprintf("iteration_%03d",i)))
    old <- s$iterations[[as.character(i)]]
    if (isTRUE(old$retired)) next
    genes <- ctx$manifest$gene
    for (j in seq_along(genes)) {
      g <- genes[j]
      if (g %in% ctx$scope$excluded) next
      x <- old$summaries[[g]]
      before <- old$result_stamps[j,,drop=FALSE]; rownames(before) <- NULL
      prev_before <- if (i) old$previous_stamps[j,,drop=FALSE] else spp_stamp(character())
      rownames(prev_before) <- NULL
      good <- spp_valid_summary(x,ctx,g,before,prev_before)
      rel <- paste0("results_slide_prior_em/",basename(ctx$dir),"/results/",g,".rds")
      receipt <- ledger[ledger$relative == rel,,drop=FALSE]
      if (file.exists(sre_result(ctx,g)) && nrow(receipt) != 1L) stop("Missing verified receipt: ",g)
      if (good && file.exists(sre_result(ctx,g))) {
        if (x$row$result_md5 != receipt$old_md5) stop("Cached content hash differs from original file: ",i," / ",g)
        if (x$row$status == "reviewed_exclusion" && receipt$old_md5 != receipt$new_md5 &&
            g %in% sre_exclusions(ctx)$gene) stop("Reviewed exclusion requires a separate hash migration: ",g)
        x$result_stamp <- spp_stamp(sre_result(ctx,g))
        x$previous_stamp <- spp_stamp(spp_previous_file(ctx,g))
        x$row$result_md5 <- receipt$new_md5
        save_checked(x,spp_summary_file(ctx,g))
      } else if (!is.null(old$cert)) stop("Completed iteration has invalid original summary: ",i," / ",g)
      else if (file.exists(sre_result(ctx,g))) spp_inspect(ctx,g,force=TRUE)
    }
    index <- spp_index(ctx)
    ready <- vapply(index$entries,function(x) !is.null(x) && sre_ready(x$row$status),logical(1))
    if (!is.null(old$cert)) {
      if (!all(ready)) stop("Completed iteration lost validation")
      if (!identical(old$pooled_md5,md5(file.path(spp_root(ctx),"pooled.rds")))) stop("Pooled counts changed")
      cert <- old$cert
      cert$result_stamps <- spp_stamp(sre_result(ctx,genes))
      cert$previous_stamps <- if(i) spp_stamp(vapply(genes,function(g) spp_previous_file(ctx,g),character(1))) else NULL
      cert$summary_stamps <- spp_stamp(vapply(genes,function(g) spp_summary_file(ctx,g),character(1)))
      save_checked(cert,file.path(spp_root(ctx),"COMPLETE.rds"))
      spp_require_complete(project,ctx$dir)
    }
    audit_file <- file.path(ctx$root,"audit.csv")
    if (file.exists(audit_file)) {
      a <- read.csv(audit_file,stringsAsFactors=FALSE)
      for (g in intersect(a$gene,genes)) {
        x <- index$entries[[g]]
        if (!is.null(x)) a$result_md5[a$gene==g] <- x$row$result_md5
      }
      em_atomic_write(a,audit_file,csv=TRUE)
    }
    if (i == s$latest) {
      reviewed <- list()
      archive <- file.path(state,"quota_attempt_receipts.rds")
      if(file.exists(archive)) reviewed <- readRDS(archive)
      for(g in genes[!ready]) {
        a <- sre_read(sre_attempt_file(ctx,g))
        if (!is.null(a$error) && grepl("error writing|Disk quota|No space",a$error,ignore.case=TRUE)) {
          reviewed[[g]] <- a
          save_checked(reviewed,archive)
          unlink(sre_attempt_file(ctx,g))
        }
      }
      # All workers were proven terminal by the Python maintenance lock/queue check.
      # A quota-full filesystem can prevent the initial started receipt itself.
      batches <- list.files(file.path(spp_root(ctx),"batches"),"^array_job_id\\.txt$",recursive=TRUE,full.names=TRUE)
      for (f in batches) {
        spp_batch(ctx,dirname(f))
        save_checked(TRUE,file.path(dirname(f),"closed.rds"))
      }
      save_checked(list(reviewed_batches=dirname(batches),reason="Disk-quota recovery; no live slider jobs"),
                   file.path(state,"closed_batches.rds"))
    }
    report[[as.character(i)]] <- data.frame(iteration=i,existing_results=sum(file.exists(sre_result(ctx,genes))),
       validated_genes=sum(ready & !genes %in% ctx$scope$excluded),
       scope_excluded=length(ctx$scope$excluded),pending=sum(!ready),complete=!is.null(old$cert))
    cat("Repaired metadata for iteration",i,"; validated in-scope genes:",sum(ready & !genes %in% ctx$scope$excluded),"\n")
  }
  assert_history(s)
  write.csv(do.call(rbind,report),file.path(state,"validation_report.csv"),row.names=FALSE)
  cat("Prior histories, frozen priors and pooled counts preserved. No gene was fitted.\n")
}
retire_prepare <- function(bundle,project,state,through,expected) {
  snapshot(bundle,project,state,expected)
  s <- readRDS(file.path(state,"snapshot.rds"))
  last_complete <- !is.null(s$iterations[[as.character(s$latest)]]$cert)
  limit <- s$latest - if(last_complete) 1L else 2L
  if(through < 0L || through > limit)
    stop("Retirement limit is iteration ",limit,"; retain the warm-start checkpoint")
  save_checked(list(through=through,latest=s$latest),file.path(state,"retire_plan.rds"))
}
retire_finish <- function(bundle,project,state) {
  load_helpers(bundle)
  s <- readRDS(file.path(state,"snapshot.rds")); assert_history(s)
  plan <- readRDS(file.path(state,"retire_plan.rds"))
  i <- plan$through+1L
  ctx <- sre_context(project,file.path(project,"results_slide_prior_em",sprintf("iteration_%03d",i)))
  old <- s$iterations[[as.character(i)]]
  genes <- ctx$manifest$gene
  if(!identical(old$result_stamps,spp_stamp(sre_result(ctx,genes)))) stop("Retained results changed during retirement")
  for(j in seq_along(genes)) {
    g <- genes[j]
    if(g %in% ctx$scope$excluded) next
    x <- old$summaries[[g]]
    before <- old$result_stamps[j,,drop=FALSE]; rownames(before)<-NULL
    prev <- old$previous_stamps[j,,drop=FALSE]; rownames(prev)<-NULL
    if(!spp_valid_summary(x,ctx,g,before,prev)) stop("Cannot retire dependency of an invalid summary: ",g)
    x$previous_stamp <- spp_stamp(spp_previous_file(ctx,g))
    save_checked(x,spp_summary_file(ctx,g))
  }
  spp_index(ctx)
  if(!is.null(old$cert)) {
    cert <- old$cert
    cert$previous_stamps <- spp_stamp(vapply(genes,function(g) spp_previous_file(ctx,g),character(1)))
    cert$summary_stamps <- spp_stamp(vapply(genes,function(g) spp_summary_file(ctx,g),character(1)))
    save_checked(cert,file.path(spp_root(ctx),"COMPLETE.rds"))
    spp_require_complete(project,ctx$dir)
  }
  assert_history(s)
  cat("Retained iteration",i,"still validates; prior histories are unchanged.\n")
}
if (sys.nframe()==0L) {
  a <- commandArgs(TRUE)
  if(a[1]=="file" && length(a)==4L) compact_file(a[2],a[3],a[4])
  else if(a[1]=="snapshot" && length(a)==5L) snapshot(a[2],a[3],a[4],as.integer(a[5]))
  else if(a[1]=="finalize" && length(a)==4L) finalize(a[2],a[3],a[4])
  else if(a[1]=="retire-prepare" && length(a)==6L) retire_prepare(a[2],a[3],a[4],as.integer(a[5]),as.integer(a[6]))
  else if(a[1]=="retire-finish" && length(a)==4L) retire_finish(a[2],a[3],a[4])
  else stop("Invalid maintenance arguments")
}
