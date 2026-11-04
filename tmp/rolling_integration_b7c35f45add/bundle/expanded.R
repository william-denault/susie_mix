seed_expanded <- function(bundle, source_project, branch, source_iteration=5L) {
  if(source_iteration != 5L) stop("This branch must start at original 005")
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
  parent_root <- file.path(source_project,"results_slide_prior_em")
  if(max(read.csv(file.path(parent_root,"prior_history.csv"))$iteration)!=5L) stop("Original run must stop at 005")
  parent <- sre_context(source_project,file.path(parent_root,"iteration_005"),check_package=TRUE)
  spp_require_complete(source_project,parent$dir)
  if(length(intersect(map$gene,parent$manifest$gene))) stop("New output key collides with original cohort")
  in_scope <- setdiff(parent$manifest$gene,parent$scope$excluded)
  entries <- setNames(lapply(in_scope,function(g) {
    x <- sre_read(spp_summary_file(parent,g))
    if(!spp_valid_summary(x,parent,g,spp_stamp(sre_result(parent,g)),spp_stamp(spp_previous_file(parent,g))) || !sre_ready(x$row$status))
      stop("Source gene is not validated: ",g)
    x
  }),in_scope)
  before <- spp_stamp(sre_result(parent,in_scope))
  root <- file.path(branch,"results_slide_prior_em"); d <- file.path(root,"iteration_005")
  dir.create(file.path(d,"results"),recursive=TRUE,showWarnings=FALSE)
  manifest <- em_rechunk_manifest(rbind(parent$manifest[c("gene","chunk")],data.frame(gene=map$gene,chunk=1L)),298L)
  scope <- readRDS(file.path(parent_root,"chromosome_scope.rds"))
  chromosomes <- c(setNames(scope$genes$chromosome,scope$genes$gene),setNames(as.character(map$chromosome),map$gene))
  scope$genes <- data.frame(gene=manifest$gene,chromosome=unname(chromosomes[manifest$gene]))
  if(anyNA(scope$genes$chromosome)) stop("Missing chromosome mapping")
  em_atomic_write(scope,file.path(root,"chromosome_scope.rds"))
  em_atomic_write(manifest,file.path(d,"manifest.csv"),csv=TRUE)
  for(f in c("priors.csv","settings.rds","grid.csv","prior_history.csv"))
    if(!file.copy(file.path(parent$dir,f),file.path(d,f),overwrite=TRUE)) stop("Cannot copy frozen settings/history")
  for(f in c("prior_history.csv","original_prior_history.csv"))
    if(!file.copy(file.path(parent_root,"prior_history.csv"),file.path(root,f),overwrite=TRUE)) stop("Cannot copy prior history")
  # The source is certified with the ORIGINAL helpers. Only the new branch
  # treats 005 as its cohort boundary; original outputs are not rewritten.
  for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R","slide_prior_parallel.R"))
    sys.source(file.path(branch,"script/scan_tissue_attempt",f),globalenv())
  ctx <- sre_context(branch,d,check_package=TRUE); spp_init(ctx)
  seed_hashes <- character()
  for(g in in_scope) {
    original <- sre_result(parent,g); target <- sre_result(ctx,g)
    if(!file.exists(target)) {
      if(!file.link(original,target)) stop("Cannot hard-link seed. Branch must use the same filesystem: ",g)
    } else {
      # A partial init may be retried. Python checks same-inode ownership before
      # source retirement; checksum/stamps also protect this reuse here.
      if(md5(target)!=entries[[g]]$row$result_md5) stop("Existing seed output differs: ",g)
    }
    x <- entries[[g]]; seed_hashes[g] <- x$row$result_md5
    x$key <- spp_key(ctx); x$result_stamp <- spp_stamp(target); x$previous_stamp <- spp_stamp(character())
    em_atomic_write(x,spp_summary_file(ctx,g))
  }
  if(!identical(before,spp_stamp(sre_result(parent,in_scope)))) stop("Original results changed during initialization")
  provenance <- list(source_project=normalizePath(source_project,winslash="/"),source_iteration=5L,
    source_priors_md5=md5(file.path(parent$dir,"priors.csv")),source_frozen=parent$frozen,
    result_md5=seed_hashes,map_md5=checks$map_md5,datadir=checks$datadir,packages=checks$packages,
    reused_genes=length(in_scope),added_expression_ids=nrow(map),added_symbols=length(unique(map$expression_symbol)),
    link_type="hardlink",cohort_boundary_iteration=5L)
  em_atomic_write(provenance,file.path(root,"seed_provenance.rds"))
  index <- spp_index(ctx)
  if(sum(vapply(index$entries,function(x) !is.null(x) && sre_ready(x$row$status),logical(1))) != nrow(parent$manifest))
    stop("Seed reuse count differs from original cohort")
  write.csv(data.frame(source_iteration=5L,reused_genes=length(in_scope),added_expression_ids=nrow(map),
    in_scope_entries=length(in_scope)+nrow(map),initial_prior_md5=md5(file.path(d,"priors.csv"))),
    file.path(branch,"seed_report.csv"),row.names=FALSE)
  cat("Expanded iteration 005:",length(in_scope),"existing fits reused;",nrow(map),"new IDs await fitting under the same priors.\n")
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE)
  if(length(a)!=4L) stop("Usage: expanded.R BUNDLE SOURCE BRANCH 5")
  seed_expanded(a[1],a[2],a[3],as.integer(a[4]))
}
