.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
root <- readLines("tmp/slider_expansion_fixture_path.txt")
bundle <- file.path(root,"test_bundle"); project <- file.path(root,"project");branch <- file.path(root,"expanded")
source(file.path(bundle,"maintenance.R"));load_helpers(bundle)
parent <- sre_context(project,file.path(project,"results_slide_prior_em/iteration_002"))
file.copy(sre_result(parent,"GENEA"),sre_result(parent,"GENEB"),overwrite=TRUE)
spp_inspect(parent,"GENEB",force=TRUE);spp_complete(parent,spp_index(parent))
dir.create(file.path(bundle,"checks"),showWarnings=FALSE)
input <- file.path(project,"results_slide_prior_em/iteration_000/manifest.csv")
info <- file.info(input)[,c("size","mtime"),drop=FALSE];info$path<-input;rownames(info)<-NULL
map <- read.csv(file.path(bundle,"renamed_gene_map.csv"),stringsAsFactors=FALSE)
saveRDS(list(passed=TRUE,project=normalizePath(project,winslash="/"),map_md5=md5(file.path(bundle,"renamed_gene_map.csv")),
             input_info=info,packages=c(susieRSlidePrior=spe_check_packages()),datadir=root),file.path(bundle,"checks/preflight.rds"))
source(file.path(bundle,"expanded.R"))
if(.Platform$OS.type=="windows") {
  file.symlink <- function(from,to) file.copy(from,to)
  cat("Windows fixture: copy stand-in for Linux symlinks; identical source bytes.\n")
}
before <- tools::md5sum(c(sre_result(parent,"GENEA"),sre_result(parent,"GENEB")))
seed_expanded(bundle,project,branch,2L)
for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R","slide_prior_parallel.R","prepare_slide_prior_em_iteration.R"))
  source(file.path(branch,"script/scan_tissue_attempt",f))
ctx <- sre_context(branch,file.path(branch,"results_slide_prior_em/iteration_000"))
stopifnot(nrow(ctx$manifest)==683L)
# Force a saved original-iteration result through the private seed reader.
unlink(spp_summary_file(ctx,"GENEA"))
stopifnot(spp_inspect(ctx,"GENEA",force=TRUE)$status=="valid")
q <- file.path(root,"expanded_queue.txt");writeLines(character(),q)
plan <- spp_plan(ctx,298L,q);work<-read.csv(file.path(plan$batch,"manifest.csv"))
stopifnot(nrow(work)==298L,all(work$gene %in% map$gene),all(work$action=="fit"))
stopifnot(inherits(try(spp_complete(ctx,spp_index(ctx)),silent=TRUE),"try-error"))
cat("PASS: original fits reused, only new IDs scheduled; incomplete expansion blocks the M-step.\n")
# One new informative fit + 679 explicit no-eligible-tissue outputs.
added <- readRDS(sre_result(parent,"GENEA"));attr(added,"em_iteration")<-0L
for(j in seq_len(nrow(map))) {
  g<-map$gene[j]
  out<-if(j==1L) added else structure(list(),em_iteration=0L)
  em_atomic_write(out,sre_result(ctx,g));spp_inspect(ctx,g,force=TRUE)
}
pooled<-spp_complete(ctx,spp_index(ctx))
stopifnot(sum(pooled$audit$issue=="fit")==3L)
stopifnot(max(abs(as.numeric(pooled$counts[1,spe_columns()])-9/17))<1e-12)
launch <- spe_prepare_iteration(branch,"new",completion_reader=spp_require_complete)
stopifnot(basename(launch$iteration_dir)=="iteration_001",identical(before,tools::md5sum(names(before))))
next_ctx<-sre_context(branch,launch$iteration_dir)
stopifnot(identical(sre_previous(next_ctx,"GENEA"),readRDS(sre_result(parent,"GENEA"))))
cat("PASS: added fit contributes to the unchanged M-step; next iteration warm-starts existing fits; originals unchanged.\n")
# Test ID dispatch independently of PLINK and fitting.
captured<-NULL
run_slide_prior_gene<-function(target_gene,target_gene_id=NULL,...) {captured<<-list(symbol=target_gene,id=target_gene_id);list()}
source(file.path(bundle,"expanded_worker_wrapper.R"))
for(g in map$gene[map$expression_symbol=="C2orf61"]) {
  run_slide_prior_gene(g,ctx$priors,branch)
  row<-map[map$gene==g,,drop=FALSE]
  stopifnot(captured$symbol=="C2orf61",captured$id==row$stable_id)
}
run_slide_prior_gene("GENEA",ctx$priors,branch)
stopifnot(captured$symbol=="GENEA",is.null(captured$id))
cat("PASS: duplicate historical symbols dispatch to their separate stable IDs; original-gene selection unchanged.\n")
