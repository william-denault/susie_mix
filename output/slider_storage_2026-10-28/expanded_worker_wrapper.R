# Only the added genes use stable-ID matching. Original genes retain their
# original expression/annotation selection and numerical fitting settings.
.expanded_run_gene <- run_slide_prior_gene
run_slide_prior_gene <- function(target_gene, tissue_priors, project_dir, ...) {
  root <- file.path(project_dir,"results_slide_prior_em")
  provenance <- readRDS(file.path(root,"seed_provenance.rds"))
  path <- file.path(project_dir,"renamed_gene_map.csv")
  if(unname(tools::md5sum(path)) != provenance$map_md5) stop("Expanded gene map changed")
  map <- read.csv(path,stringsAsFactors=FALSE)
  row <- map[map$gene==target_gene,,drop=FALSE]
  args <- c(list(target_gene=if(nrow(row)) row$expression_symbol else target_gene,
                 tissue_priors=tissue_priors,project_dir=project_dir,datadir=provenance$datadir),list(...))
  if(nrow(row)) args$target_gene_id <- row$stable_id
  out <- do.call(.expanded_run_gene,args)
  if(nrow(row)) attr(out,"recovery_identity") <- as.list(row)
  out
}
