si_metadata <- function(project) file.path(project,"results_slide_prior_em/inplace_expansion_2026-10-29")
si_info <- function(project) readRDS(file.path(si_metadata(project),"provenance.rds"))
si_new_at_six <- function(ctx,gene) ctx$iteration==6L && gene %in% si_info(ctx$project)$added_genes
