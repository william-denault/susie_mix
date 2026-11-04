seed_expanded <- function(bundle, source_project, branch, source_iteration) {
  source(file.path(bundle,"maintenance.R")); load_helpers(bundle)
  checks <- readRDS(file.path(bundle,"checks/preflight.rds"))
  if(!isTRUE(checks$passed) || checks$project != normalizePath(source_project,winslash="/")) stop("Preflight has not passed for this project")
  if(checks$map_md5 != md5(file.path(bundle,"renamed_gene_map.csv"))) stop("Gene map changed")
  now <- file.info(checks$input_info$path)
  if(anyNA(now$size) || any(now$size != checks$input_info$size) ||
     any(as.numeric(now$mtime) != as.numeric(checks$input_info$mtime))) stop("Preflight inputs changed")
  if(checks$packages[["susieRSlidePrior"]] != spe_check_packages()) stop("Slider package changed")
  map <- read.csv(file.path(bundle,"renamed_gene_map.csv"),stringsAsFactors=FALSE)
  if(nrow(map)!=680L || anyDuplicated(map$gene) || anyDuplicated(map$stable_id)) stop("Invalid candidate map")
  parent <- sre_context(source_project,file.path(source_project,"results_slide_prior_em",sprintf("iteration_%03d",source_iteration)),check_package=TRUE)
  spp_require_complete(source_project,parent$dir)
  if(length(intersect(map$gene,parent$manifest$gene))) stop("New output key collides with original cohort")
  root <- file.path(branch,"results_slide_prior_em"); d <- file.path(root,"iteration_000")
  dir.create(file.path(d,"results"),recursive=TRUE,showWarnings=FALSE)
  priors <- parent$priors
  priors$iteration <- 0L; priors$source_iteration <- -1L
  priors$update_method <- "same_prior_seed_from_original_slider"
  priors$created_at <- format(Sys.time(),tz="UTC",usetz=TRUE)
  manifest <- rbind(parent$manifest[c("gene","chunk")],data.frame(gene=map$gene,chunk=1L))
  manifest <- em_rechunk_manifest(manifest,298L)
  scope <- readRDS(file.path(source_project,"results_slide_prior_em/chromosome_scope.rds"))
  chromosomes <- c(setNames(scope$genes$chromosome,scope$genes$gene),setNames(as.character(map$chromosome),map$gene))
  scope$genes <- data.frame(gene=manifest$gene,chromosome=unname(chromosomes[manifest$gene]))
  if(anyNA(scope$genes$chromosome)) stop("Missing chromosome mapping")
  em_atomic_write(scope,file.path(root,"chromosome_scope.rds"))
  em_atomic_write(manifest,file.path(d,"manifest.csv"),csv=TRUE)
  em_atomic_write(priors,file.path(d,"priors.csv"),csv=TRUE)
  em_atomic_write(priors,file.path(d,"prior_history.csv"),csv=TRUE)
  em_atomic_write(priors,file.path(root,"prior_history.csv"),csv=TRUE)
  for(f in c("settings.rds","grid.csv")) {
    if(!file.copy(file.path(parent$dir,f),file.path(d,f),overwrite=TRUE)) stop("Cannot copy frozen settings")
  }
  file.copy(file.path(source_project,"results_slide_prior_em/prior_history.csv"),
            file.path(root,"original_prior_history.csv"),overwrite=TRUE)
  ctx <- sre_context(branch,d,check_package=TRUE); spp_init(ctx)
  seed_hashes <- character()
  in_scope <- setdiff(parent$manifest$gene,parent$scope$excluded)
  for(g in in_scope) {
    original <- sre_result(parent,g); target <- sre_result(ctx,g)
    x <- sre_read(spp_summary_file(parent,g))
    if(!spp_valid_summary(x,parent,g,spp_stamp(original),spp_stamp(spp_previous_file(parent,g))) || !sre_ready(x$row$status))
      stop("Source gene is not validated: ",g)
    link <- Sys.readlink(target)
    if(!is.na(link) && nzchar(link)) {
      if(normalizePath(target,winslash="/") != normalizePath(original,winslash="/")) stop("Wrong seed link: ",g)
    } else if(file.exists(target)) stop("Existing seed output is not a source link: ",g)
    else if(!file.symlink(original,target)) stop("Cannot link source result: ",g)
    seed_hashes[g] <- x$row$result_md5
    x$key <- spp_key(ctx); x$result_stamp <- spp_stamp(target); x$previous_stamp <- spp_stamp(character())
    x$validated_at <- format(Sys.time(),tz="UTC",usetz=TRUE)
    em_atomic_write(x,spp_summary_file(ctx,g))
  }
  provenance <- list(source_project=normalizePath(source_project,winslash="/"),source_iteration=source_iteration,
       source_priors_md5=md5(file.path(parent$dir,"priors.csv")),source_frozen=parent$frozen,
       result_md5=seed_hashes,map_md5=checks$map_md5,datadir=checks$datadir,packages=checks$packages,
       reused_genes=length(in_scope),added_expression_ids=nrow(map),added_symbols=length(unique(map$expression_symbol)))
  em_atomic_write(provenance,file.path(root,"seed_provenance.rds"))
  # Keys on current summaries refer to the same frozen weights, not the prior
  # update proposed by the original completed checkpoint.
  index <- spp_index(ctx)
  if(sum(vapply(index$entries,function(x) !is.null(x) && sre_ready(x$row$status),logical(1))) != nrow(parent$manifest))
    stop("Seed reuse count differs from the original cohort")
  spp_require_complete(source_project,parent$dir)
  write.csv(data.frame(source_iteration=source_iteration,reused_genes=length(in_scope),
      added_expression_ids=nrow(map),in_scope_entries=length(in_scope)+nrow(map),
      initial_prior_md5=md5(file.path(d,"priors.csv"))),file.path(branch,"seed_report.csv"),row.names=FALSE)
  cat("Expanded seed prepared:",length(in_scope),"existing genes reused;",nrow(map),"new expression IDs await fitting.\n")
}
if(sys.nframe()==0L) {
  a <- commandArgs(TRUE)
  if(length(a)!=4L) stop("Usage: expanded.R BUNDLE SOURCE_PROJECT BRANCH SOURCE_ITERATION")
  seed_expanded(a[1],a[2],a[3],as.integer(a[4]))
}
