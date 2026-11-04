# Small metadata only. No fitting and no full-posterior decompression.
sr_load <- function(project) {
  runtime<-file.path(project,"slider_inplace_runtime_2026-10-29")
  source(file.path(runtime,"maintenance.R"))
  for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R","slide_prior_parallel.R","cohort.R"))
    sys.source(file.path(runtime,"script/scan_tissue_attempt",f),globalenv())
}
sr_history <- function(project) {
  files <- history_files(project)
  root <- file.path(project,"results_slide_prior_em")
  other <- list.files(root,"^(pooled\\.rds|audit\\.csv|original_prior_history\\.csv|seed_provenance\\.rds)$",recursive=TRUE,full.names=TRUE)
  unique(c(files,other))
}
sr_verify <- function(project,state) {
  sr_load(project); s <- readRDS(file.path(state,"snapshot.rds"))
  if(s$project != normalizePath(project,winslash="/")) stop("Wrong retention project")
  assert_history(s)
  ctx <- sre_context(project,s$keep_dir,check_package=TRUE)
  if(!identical(s$result_stamps,spp_stamp(sre_result(ctx,ctx$manifest$gene)))) stop("Retained checkpoint changed")
  if(md5(file.path(spp_root(ctx),"pooled.rds")) != s$cert$pooled_md5) stop("Retained pooled counts changed")
  current <- sre_context(project,s$current_dir)
  if(!identical(current$frozen,s$current_frozen) || max(read.csv(file.path(project,"results_slide_prior_em/prior_history.csv"))$iteration)!=s$current)
    stop("Current iteration changed")
  if(length(list.files(file.path(current$dir,"results"),all.files=TRUE,no..=TRUE))) stop("Current workers have already written results")
  invisible(s)
}
sr_prepare <- function(project,state,current) {
  sr_load(project)
  if(current < 6L || current > 40L) stop("Cleanup only runs before 006..040")
  root <- file.path(project,"results_slide_prior_em")
  if(max(read.csv(file.path(root,"prior_history.csv"))$iteration)!=current) stop("Not the latest iteration")
  ctx <- sre_context(project,file.path(root,sprintf("iteration_%03d",current-1L)),check_package=TRUE)
  spp_require_complete(project,ctx$dir)
  next_ctx <- sre_context(project,file.path(root,sprintf("iteration_%03d",current)))
  if(length(list.files(file.path(next_ctx$dir,"results"),all.files=TRUE,no..=TRUE))) stop("Current results already exist")
  genes <- ctx$manifest$gene
  entries <- setNames(lapply(genes,function(g) sre_read(spp_summary_file(ctx,g))),genes)
  cert <- readRDS(file.path(spp_root(ctx),"COMPLETE.rds"))
  for(g in setdiff(genes,ctx$scope$excluded))
    if(!spp_valid_summary(entries[[g]],ctx,g,spp_stamp(sre_result(ctx,g)),spp_stamp(spp_previous_file(ctx,g))) || !sre_ready(entries[[g]]$row$status))
      stop("Invalid retained summary: ",g)
  dirs <- file.path(root,sprintf("iteration_%03d",seq.int(0L,current-2L)))
  if(any(!dir.exists(dirs))) stop("Missing historical iteration directory")
  preserved <- tools::md5sum(sr_history(project))
  if(anyNA(preserved)) stop("Unreadable preserved history")
  s <- list(project=normalizePath(project,winslash="/"),current=current,current_dir=next_ctx$dir,
    current_frozen=next_ctx$frozen,keep_dir=ctx$dir,cert=cert,summaries=entries,
    result_stamps=spp_stamp(sre_result(ctx,genes)),
    previous_stamps=if(ctx$iteration>0L) spp_stamp(vapply(genes,function(g) spp_previous_file(ctx,g),character(1))) else NULL,
    preserved=preserved,retire_dirs=dirs)
  dir.create(state,recursive=TRUE,showWarnings=FALSE)
  save_checked(s,file.path(state,"snapshot.rds"))
  write.csv(data.frame(directory=dirs),file.path(state,"directories.csv"),row.names=FALSE)
  sr_verify(project,state)
  cat("Retirement prepared before",current,"; keep",current-1L,"and",current,"full outputs.\n")
}
sr_finish <- function(project,state) {
  s <- sr_verify(project,state)
  for(d in s$retire_dirs)
    if(length(list.files(file.path(d,"results"),all.files=TRUE,no..=TRUE))) stop("Old results remain: ",d)
  ctx <- sre_context(project,s$keep_dir); genes <- ctx$manifest$gene
  if(ctx$iteration>0L) for(j in seq_along(genes)) {
    g <- genes[j]; if(g %in% ctx$scope$excluded) next
    x <- s$summaries[[g]]
    before <- s$result_stamps[j,,drop=FALSE]; rownames(before)<-NULL
    prev <- s$previous_stamps[j,,drop=FALSE]; rownames(prev)<-NULL
    if(!spp_valid_summary(x,ctx,g,before,prev)) stop("Invalid saved summary: ",g)
    x$previous_stamp <- spp_stamp(spp_previous_file(ctx,g))
    save_checked(x,spp_summary_file(ctx,g))
  }
  spp_index(ctx)
  cert <- s$cert
  cert$previous_stamps <- if(ctx$iteration>0L) spp_stamp(vapply(genes,function(g) spp_previous_file(ctx,g),character(1))) else NULL
  cert$summary_stamps <- spp_stamp(vapply(genes,function(g) spp_summary_file(ctx,g),character(1)))
  save_checked(cert,file.path(spp_root(ctx),"COMPLETE.rds"))
  spp_require_complete(project,ctx$dir)
  assert_history(s)
  cat("Retained checkpoint validates; every prior/history/count file is unchanged.\n")
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE)
  if(a[1]=="prepare" && length(a)==4L) sr_prepare(a[2],a[3],as.integer(a[4]))
  else if(a[1]=="verify" && length(a)==3L) sr_verify(a[2],a[3])
  else if(a[1]=="finish" && length(a)==3L) sr_finish(a[2],a[3])
  else stop("Invalid retention command")
}
