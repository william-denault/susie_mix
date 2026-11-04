# Original genes keep their original symbol-based input selection. Added genes
# use the verified stable IDs; all fit settings and EM math are unchanged.
.inplace_run_gene <- run_slide_prior_gene
run_slide_prior_gene <- function(target_gene,tissue_priors,project_dir,...) {
  provenance<-si_info(project_dir)
  path<-file.path(project_dir,"slider_inplace_runtime_2026-10-29/renamed_gene_map.csv")
  if(unname(tools::md5sum(path))!=provenance$map_md5) stop("Expanded map changed")
  map<-read.csv(path,stringsAsFactors=FALSE);row<-map[map$gene==target_gene,,drop=FALSE]
  args<-c(list(target_gene=if(nrow(row)) row$expression_symbol else target_gene,
              tissue_priors=tissue_priors,project_dir=project_dir,datadir=provenance$datadir),list(...))
  if(nrow(row)) args$target_gene_id<-row$stable_id
  out<-do.call(.inplace_run_gene,args)
  if(nrow(row)) attr(out,"recovery_identity")<-as.list(row)
  out
}
