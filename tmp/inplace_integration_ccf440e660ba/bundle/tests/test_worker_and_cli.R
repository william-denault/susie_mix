.libPaths(c("tmp/slider_prior_impl/library","C:/Users/willi/AppData/Local/R/win-library/4.6",.libPaths()))
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
root<-trimws(readLines("tmp/slider_646_fixture_path.txt"));project<-file.path(root,"project")
runtime<-file.path(project,"slider_inplace_runtime_2026-10-29");code<-file.path(runtime,"script/scan_tissue_attempt")
for(f in c("em_utils.R","slide_prior_em_utils.R","slide_prior_recovery.R","slide_prior_parallel.R","cohort.R","workhorse_slide_prior_em.R")) source(file.path(code,f))
map<-read.csv(file.path(runtime,"renamed_gene_map.csv"));ctx<-sre_context(project,file.path(project,"results_slide_prior_em/iteration_008"))
.inplace_run_gene<-function(...) list(arguments=list(...))
old<-run_slide_prior_gene("GENEA",ctx$priors,project,previous_result=list())
stopifnot(old$arguments$target_gene=="GENEA",is.null(old$arguments$target_gene_id))
for(g in unique(c(map$gene[1],map$gene[map$expression_symbol %in% c("C2orf61","CTAGE5")]))) {
 row<-map[map$gene==g,];new<-run_slide_prior_gene(g,ctx$priors,project,previous_result=NULL)
 stopifnot(new$arguments$target_gene==row$expression_symbol,new$arguments$target_gene_id==row$stable_id,
  identical(attr(new,"recovery_identity"),as.list(row)),identical(new$arguments$tissue_priors,ctx$priors))
}
cat("PASS: existing symbols retain selection; renamed/duplicate symbols use distinct frozen stable IDs.\n")
for(g in c("GENEA",map$gene[1])) {
 file<-sre_result(ctx,g);saved<-readRDS(file)
 em_atomic_write(structure(list(),em_iteration=8L),file)
 status<-spp_inspect(ctx,g,force=TRUE)
 stopifnot(!sre_ready(status$status))
 em_atomic_write(saved,file);stopifnot(spp_inspect(ctx,g,force=TRUE)$status=="valid")
}
invisible(spp_complete(ctx,spp_index(ctx)))
cat("PASS: dropping an original or added gene's established tissue fits is rejected.\n")
Sys.setenv(SUSIE_SLIDE_TARGET_ITERATION="8")
queue<-file.path(root,"queue.txt")
out<-system2(Sys.getenv("RSCRIPT"),c("--vanilla",shQuote(file.path(code,"run_slide_prior_parallel.R")),"plan",shQuote(project),"latest","0",shQuote(queue)),stdout=TRUE)
stopifnot(is.null(attr(out,"status")),identical(out,c("complete","8")))
cat("PASS: actual private CLI loads correctly, recognizes validated results, and stops at target.\n")
