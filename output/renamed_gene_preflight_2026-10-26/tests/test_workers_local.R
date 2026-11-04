# Local regression: same expression ID produces the same fits after a rename.
# Uses the repository's existing synthetic GTEx/PLINK fixture, not RCC data.
.libPaths(c("tmp/slider_prior_impl/library", "C:/Users/willi/AppData/Local/R/win-library/4.6", .libPaths()))
fixture_source <- readLines("script/scan_tissue_attempt/tests/test_slide_prior_em_worker.R")
# Use its real-fitter fixture, excluding the separate recovery-manager tests.
fixture_end <- grep("^# Integrate the recovery wrapper",fixture_source)[1]-1L
eval(parse(text=fixture_source[seq_len(fixture_end)]),envir=globalenv())
bundle <- normalizePath("output/renamed_gene_preflight_2026-10-26",winslash="/")
repo <- normalizePath(".",winslash="/")
expression2 <- expression[2:1,,drop=FALSE]
expression2$Description <- "GENE" # Duplicate old symbol; the correct ID is now second.
mapped_expr <- file.path(root,"mapped_expression.gct")
writeLines(c("#1.2",paste(2,n,sep="\t")),mapped_expr)
suppressWarnings(write.table(expression2,mapped_expr,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,append=TRUE))
mapped_annot <- file.path(root,"mapped_annotations.R")
writeLines('get_gene_annotations <- function(path) data.frame(gene_name=c("NEW_NAME","OTHER"),ensembl=c("E1","E2"),chromosome="1",start=1000,end=2000,strand="+")',mapped_annot)
mock_plink <- worker_env$system
uniform <- data.frame(tissue=c("Brain","Liver"),matrix(1/17,2,17,dimnames=list(NULL,spe_columns())))
common <- list(target_gene="GENE",datadir=root,subject_pheno_file=file.path(root,"subjects.tsv"),
               sample_attr_file=file.path(root,"samples.tsv"),L=1L,max_iter=2000)
for(model in c("mix","slider")) {
  fname <- if(model=="mix") "workhorse_em.R" else "workhorse_slide_prior_em.R"
  fun <- if(model=="mix") "run_susie_gene" else "run_slide_prior_gene"
  a <- new.env(parent=globalenv());b <- new.env(parent=globalenv())
  a$system <- b$system <- mock_plink
  a$susie <- susieR::susie
  sys.source(file.path(repo,"script/scan_tissue_attempt",fname),a)
  sys.source(file.path(bundle,"script/scan_tissue_attempt",fname),b)
  priors <- if(model=="mix") data.frame(tissue=c("Brain","Liver"),pi_add=.5,pi_rec=.25,pi_dom=.25,source_fit="weighted_fit_mix") else uniform
  base_args <- c(common,list(project_dir=repo,tissue_priors=priors,expr_file=file.path(root,"expression.gct"),
      gene_annot_fun=file.path(root,"annotations.R"),temp_dir=file.path(root,paste0(model,"_original"))))
  expected <- do.call(a[[fun]],base_args)
  mapped_args <- base_args
  mapped_args$project_dir <- bundle; mapped_args$expr_file <- mapped_expr
  mapped_args$gene_annot_fun <- mapped_annot; mapped_args$target_gene_id <- "E1"
  mapped_args$temp_dir <- file.path(root,paste0(model,"_mapped"))
  observed <- do.call(b[[fun]],mapped_args)
  stopifnot(length(expected)==2L,identical(names(expected),names(observed)))
  fitname <- if(model=="mix") "weighted_fit_mix" else "fit_slide_prior"
  for(t in names(expected)) stopifnot(isTRUE(observed[[t]][[fitname]]$converged),
      isTRUE(all.equal(expected[[t]],observed[[t]],tolerance=1e-10)))
  cat("PASS:",model,"stable-ID selection matches original numerical fits despite duplicate old symbols and renamed annotation.\n")
}
