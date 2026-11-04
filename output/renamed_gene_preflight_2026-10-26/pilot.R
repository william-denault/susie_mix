# Six independent gene-ID pilot jobs, each fitting both models without any EM update.
pilot_genes <- function(map) {
  selected <- map[map$expression_symbol %in% c("AAED1", "WARS", "CTAGE5", "C2orf61", "HIST1H3D"), , drop=FALSE]
  selected[order(selected$gene), , drop=FALSE]
}
pilot <- function(bundle, task) {
  bundle <- normalizePath(bundle, winslash="/", mustWork=TRUE)
  source(file.path(bundle, "preflight.R"), local=TRUE)
  for (f in c("em_utils.R", "slide_prior_em_utils.R"))
    source(file.path(bundle, "script/scan_tissue_attempt", f), local=TRUE)
  state <- readRDS(file.path(bundle,"checks/preflight.rds"))
  if (!isTRUE(state$passed)) stop("Metadata preflight has not passed.")
  if (!identical(state$map_md5, unname(tools::md5sum(file.path(bundle,"renamed_gene_map.csv"))))) stop("Recovery map changed.")
  current <- file.info(state$input_info$path)
  if (anyNA(current$size) || any(current$size != state$input_info$size) ||
      any(as.numeric(current$mtime) != as.numeric(state$input_info$mtime))) stop("Preflight source metadata changed; recheck inputs.")
  genes <- pilot_genes(read_rename_map(bundle))
  if (length(task)!=1 || !is.finite(task) || task!=floor(task) || task<1 || task>nrow(genes)) stop("Invalid pilot task.")
  gene <- genes[task,,drop=FALSE]
  checkpoints <- read.csv(file.path(bundle,"checks/pilot_checkpoints.csv"), stringsAsFactors=FALSE)
  if (any(!nzchar(checkpoints$path)) || any(!file.exists(file.path(checkpoints$path,"priors.csv"))))
    stop("A completed model checkpoint is unavailable; inspect pilot_checkpoints.csv.")
  if (!file.exists(file.path(checkpoints$path[checkpoints$model=="slider"],"recovery/parallel/COMPLETE.rds")))
    stop("Slider checkpoint is not marked complete.")
  output <- file.path(bundle,"pilot",gene$gene); dir.create(output, recursive=TRUE, showWarnings=FALSE)
  report <- list(); failures <- character()
  for (model in c("mix","slider")) {
    source_file <- file.path(checkpoints$path[checkpoints$model==model],"priors.csv")
    prior_md5 <- unname(tools::md5sum(source_file))
    priors <- read.csv(source_file,stringsAsFactors=FALSE)
    tryCatch({
      if (model=="mix") em_validate_priors(priors) else spe_validate_priors(priors)
      package <- if (model=="mix") "susieR" else "susieRSlidePrior"
      if (!requireNamespace(package,quietly=TRUE) || as.character(utils::packageVersion(package))!=state$packages[[package]])
        stop("Fitting package changed or is unavailable.")
      env <- new.env(parent=environment())
      file <- if (model=="mix") "workhorse_em.R" else "workhorse_slide_prior_em.R"
      sys.source(file.path(bundle,"script/scan_tissue_attempt",file),env)
      fit_gene <- if (model=="mix") env$run_susie_gene else env$run_slide_prior_gene
      file <- file.path(output,paste0(model,".rds"))
      args <- list(target_gene=gene$expression_symbol,target_gene_id=gene$stable_id,
          tissue_priors=priors,project_dir=bundle,datadir=state$datadir,
          temp_dir=file.path(output,paste0("temp_",model)),max_iter=5000)
      reused <- file.exists(file)
      if (reused) {
        out <- readRDS(file)
        if (!identical(attr(out,"recovery_identity")$stable_id,gene$stable_id) ||
            !identical(attr(out,"source_priors")$md5,prior_md5)) stop("Existing pilot identity/prior differs; retained without overwriting.")
      } else out <- do.call(fit_gene,args)
      if (!is.null(out$error) || length(attr(out,"tissue_errors",exact=TRUE))) stop("Gene/tissue errors remain.")
      if (model=="slider") summary <- spe_summarize_gene(out,gene$gene,priors) else {
        for (t in names(out)) {
          fit <- out[[t]]$weighted_fit_mix
          expected <- em_predictor_weights(out[[t]]$mix_coding,em_prior_for_tissue(priors,t))
          if (!isTRUE(fit$converged) || !isTRUE(all.equal(unname(fit$pi),unname(expected),tolerance=1e-10)))
            stop("Unconverged mix fit or incorrect prior.")
        }
      }
      if (!identical(prior_md5,unname(tools::md5sum(source_file)))) stop("Frozen source prior changed during pilot.")
      attr(out,"recovery_identity") <- as.list(gene)
      attr(out,"source_priors") <- list(path=source_file,md5=prior_md5,iteration=unique(priors$iteration),package_version=state$packages[[package]])
      if (!reused) em_atomic_write(out,file)
      report[[model]] <- data.frame(gene=gene$gene,stable_id=gene$stable_id,model=model,
          status=if(length(out)) "PASS" else "NO_ELIGIBLE_TISSUES",n_fits=length(out),detail="")
    },error=function(e) {
      failures <<- c(failures,paste(model,conditionMessage(e)))
      report[[model]] <<- data.frame(gene=gene$gene,stable_id=gene$stable_id,model=model,status="FAIL",n_fits=NA_integer_,detail=conditionMessage(e))
    })
  }
  write.csv(do.call(rbind,report),file.path(output,"summary.csv"),row.names=FALSE)
  print(do.call(rbind,report),row.names=FALSE)
  if(length(failures)) stop(paste(failures,collapse="\n"))
}
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(length(args)!=2L) stop("Usage: Rscript pilot.R BUNDLE TASK(1..6)")
  pilot(args[1],as.numeric(args[2]))
}
