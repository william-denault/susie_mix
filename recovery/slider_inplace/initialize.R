initialize_six <- function(bundle,project) {
  source(file.path(bundle,"maintenance.R"));load_helpers(bundle)
  meta<-file.path(project,"results_slide_prior_em/inplace_expansion_2026-10-29")
  runtime<-file.path(project,"slider_inplace_runtime_2026-10-29")
  root<-file.path(project,"results_slide_prior_em")
  checks<-readRDS(file.path(bundle,"checks/preflight.rds"))
  if(!isTRUE(checks$passed) || checks$project!=normalizePath(project,winslash="/")) stop("Preflight has not passed for this project")
  if(checks$map_md5!=md5(file.path(bundle,"renamed_gene_map.csv"))) stop("Gene map changed")
  now<-file.info(checks$input_info$path)
  if(anyNA(now$size) || any(now$size!=checks$input_info$size) || any(as.numeric(now$mtime)!=as.numeric(checks$input_info$mtime))) stop("Preflight inputs changed")
  if(checks$packages[["susieRSlidePrior"]]!=spe_check_packages()) stop("Slider package changed")
  history<-read.csv(file.path(root,"prior_history.csv"),stringsAsFactors=FALSE)
  if(!max(history$iteration) %in% c(5L,6L)) stop("Setup requires completed 005 or its owned partial 006 preparation")
  parent<-sre_context(project,file.path(root,"iteration_005"),check_package=TRUE)
  pooled<-spp_require_complete(project,parent$dir)
  map<-read.csv(file.path(bundle,"renamed_gene_map.csv"),stringsAsFactors=FALSE)
  if(nrow(map)!=680L || anyDuplicated(map$gene) || anyDuplicated(map$stable_id) || length(intersect(map$gene,parent$manifest$gene))) stop("Invalid/colliding expansion map")
  # Archive only small provenance/metadata. No posterior files are copied.
  original_history<-file.path(meta,"prior_history.before_expansion.csv")
  if(!file.exists(original_history)) {
    if(max(history$iteration)!=5L) stop("Missing original history snapshot")
    if(!file.copy(file.path(root,"prior_history.csv"),original_history)) stop("Cannot preserve prior history")
  }
  before_history<-read.csv(original_history,stringsAsFactors=FALSE)
  if(!sre_same_priors(before_history,history[history$iteration<=5L,,drop=FALSE])) stop("Original history changed")
  snapshot_files<-c("priors.csv","manifest.csv","settings.rds","grid.csv","prior_history.csv",
    "source_audit.csv","pooled_counts.csv","slider_counts.csv","proposed_priors.csv",
    "recovery/parallel/pooled.rds","recovery/parallel/COMPLETE.rds")
  for(rel in snapshot_files) {
    src<-file.path(parent$dir,rel);dst<-file.path(meta,"original_005_metadata",rel)
    if(!file.exists(src)) next
    dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)
    if(file.exists(dst)) { if(md5(src)!=md5(dst)) stop("Original checkpoint metadata changed: ",rel) }
    else if(!file.copy(src,dst)) stop("Cannot preserve original checkpoint metadata")
  }
  manifest<-em_rechunk_manifest(rbind(parent$manifest[c("gene","chunk")],data.frame(gene=map$gene,chunk=1L)),298L)
  scope<-readRDS(file.path(root,"chromosome_scope.rds"))
  chromosomes<-c(setNames(scope$genes$chromosome,scope$genes$gene),setNames(as.character(map$chromosome),map$gene))
  scope$genes<-data.frame(gene=manifest$gene,chromosome=unname(chromosomes[manifest$gene]))
  if(anyNA(scope$genes$chromosome)) stop("Missing chromosome mapping")
  provenance_file<-file.path(meta,"provenance.rds")
  if(file.exists(provenance_file)) {
    provenance<-readRDS(provenance_file)
    if(!identical(provenance$original_frozen,parent$frozen) || provenance$map_md5!=checks$map_md5 ||
       !identical(provenance$original_genes,parent$manifest$gene) || !identical(provenance$added_genes,map$gene)) stop("Expansion provenance differs")
  } else {
    provenance<-list(original_iteration=5L,first_expanded_iteration=6L,first_expanded_prior_update=7L,
      original_frozen=parent$frozen,original_genes=parent$manifest$gene,added_genes=map$gene,
      map_md5=checks$map_md5,datadir=checks$datadir,packages=checks$packages,
      history_md5=md5(original_history),created_at=format(Sys.time(),tz="UTC",usetz=TRUE))
    save_checked(provenance,provenance_file)
  }
  save_expected<-function(value,path,csv=FALSE) {
    tmp<-paste0(path,".inplace_tmp")
    if(csv) write.csv(value,tmp,row.names=FALSE) else saveRDS(value,tmp)
    if(file.exists(path)) {
      if(md5(path)!=md5(tmp)) stop("Owned setup output differs: ",path)
      unlink(tmp)
    } else if(!file.rename(tmp,path)) stop("Cannot install setup output: ",path)
  }
  save_expected(manifest,file.path(meta,"manifest.csv"),TRUE)
  save_expected(scope,file.path(meta,"chromosome_scope.rds"))
  # Exactly the ordinary 005 -> 006 prior update, used as the INITIAL prior
  # for the expanded population. Added genes enter the 006 E-step, hence 007's
  # prior update. Old posteriors are warm starts, never counted as new fits.
  priors<-data.frame(iteration=6L,source_iteration=5L,pooled$priors,
    update_method="slider_prior_active_counts_v1",created_at=provenance$created_at)
  spe_validate_priors(priors)
  new_history<-rbind(before_history,priors)
  d<-file.path(root,"iteration_006")
  outputs<-c("manifest.csv","priors.csv","settings.rds","grid.csv","source_audit.csv","pooled_counts.csv","prior_history.csv")
  if(dir.exists(d)) {
    entries<-list.files(d,all.files=TRUE,no..=TRUE)
    if(any(!entries %in% c(outputs,paste0(outputs,".inplace_tmp"),"results","completed","logs"))) stop("Unexpected files in partial iteration 006")
    for(sub in c("results","completed","logs"))
      if(length(list.files(file.path(d,sub),all.files=TRUE,no..=TRUE))) stop("Iteration 006 already has work; setup will not overwrite it")
  }
  for(sub in c("results","completed","logs")) dir.create(file.path(d,sub),recursive=TRUE,showWarnings=FALSE)
  save_expected(manifest,file.path(d,"manifest.csv"),TRUE)
  save_expected(priors,file.path(d,"priors.csv"),TRUE)
  save_expected(readRDS(file.path(parent$dir,"settings.rds")),file.path(d,"settings.rds"))
  save_expected(data.frame(weight_column=spe_columns(),delta=spe_grid()),file.path(d,"grid.csv"),TRUE)
  save_expected(pooled$audit,file.path(d,"source_audit.csv"),TRUE)
  save_expected(pooled$counts,file.path(d,"pooled_counts.csv"),TRUE)
  save_expected(new_history,file.path(d,"prior_history.csv"),TRUE)
  if(max(history$iteration)==6L && !sre_same_priors(history,new_history)) stop("Iteration 006 history differs")
  em_atomic_write(new_history,file.path(root,"prior_history.csv"),csv=TRUE)
  for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R","slide_prior_parallel.R","cohort.R"))
    sys.source(file.path(runtime,"script/scan_tissue_attempt",f),globalenv())
  next_ctx<-sre_context(project,d,check_package=TRUE)
  if(nrow(next_ctx$manifest)!=nrow(parent$manifest)+680L) stop("Expanded manifest size differs")
  spp_require_complete(project,parent$dir)
  write.csv(data.frame(first_expanded_iteration=6L,existing_in_scope=nrow(parent$manifest)-length(parent$scope$excluded),
      added_expression_ids=680L,expanded_in_scope=nrow(next_ctx$manifest)-length(next_ctx$scope$excluded),
      first_prior_update_including_added_genes=7L),file.path(meta,"setup_report.csv"),row.names=FALSE)
  cat("Prepared iteration_006 in the same results folder. Original 005 fits, priors and counts preserved.\n")
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE);if(length(a)!=2L) stop("Usage: initialize.R BUNDLE PROJECT")
  initialize_six(a[1],a[2])
}
