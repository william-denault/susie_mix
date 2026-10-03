# A chunk saves full discrete fits, including the joint posterior for warm starts.
spe_run_chunk <- function(project_dir, iteration_dir, chunk, run_gene = run_slide_prior_gene) {
  priors <- read.csv(file.path(iteration_dir, "priors.csv"), stringsAsFactors = FALSE)
  spe_validate_priors(priors)
  iteration <- unique(priors$iteration)
  source_iteration <- unique(priors$source_iteration)
  if (length(iteration) != 1L || length(source_iteration) != 1L ||
      source_iteration != iteration - 1L) stop("Invalid slider iteration metadata.")
  settings <- readRDS(file.path(iteration_dir, "settings.rds"))
  settings$max_iter <- spe_fit_budget(settings$max_iter)
  manifest <- read.csv(file.path(iteration_dir, "manifest.csv"), stringsAsFactors = FALSE)
  genes <- manifest$gene[manifest$chunk == chunk]
  if (!length(genes)) stop("No genes for chunk ", chunk)
  marker <- file.path(iteration_dir, "completed", sprintf("chunk_%03d.done", chunk))
  if (file.exists(marker)) return(invisible(NULL))
  source_dir <- file.path(project_dir, "results_slide_prior_em",
                          sprintf("iteration_%03d", source_iteration), "results")
  started_at <- Sys.time()
  validate <- function(out, previous, converged = TRUE) {
    if (!is.list(out)) stop("Invalid slider gene result.")
    if (!is.null(previous)) {
      if (!is.null(previous$error)) {
        if (!identical(out$error, previous$error)) stop("Baseline gene error record changed.")
        return(invisible(TRUE))
      }
      if (!is.null(out$error) || !identical(names(out), names(previous)) ||
          !identical(attr(out, "tissue_errors", exact = TRUE),
                     attr(previous, "tissue_errors", exact = TRUE)))
        stop("Slider gene/tissue cohort changed after initialization.")
    }
    if (!is.null(out$error)) return(invisible(TRUE))
    if (length(out) && (is.null(names(out)) || anyDuplicated(names(out))))
      stop("Invalid tissue names in slider result.")
    for (tissue in names(out))
      spe_fit_counts(out[[tissue]]$fit_slide_prior, spe_prior(priors, tissue), converged)
    invisible(TRUE)
  }
  summary <- data.frame(gene = genes, status = "", n_tissue_errors = 0L, n_nonconverged = 0L)
  for (i in seq_along(genes)) {
    gene <- genes[i]
    previous <- if (iteration == 0L) NULL else readRDS(file.path(source_dir, paste0(gene, ".rds")))
    file <- file.path(iteration_dir, "results", paste0(gene, ".rds"))
    out <- if (file.exists(file)) tryCatch(readRDS(file), error = function(e) NULL) else NULL
    reusable <- is.list(out) && identical(attr(out, "em_iteration"), iteration) &&
      is.null(out$error) && !length(attr(out, "tissue_errors", exact = TRUE)) &&
      isTRUE(tryCatch(validate(out, previous), error = function(e) FALSE))
    if (!reusable) {
      cat(sprintf("[%s] slider iteration %d, chunk %03d: %s\n", Sys.time(), iteration, chunk, gene))
      if (!is.null(previous$error) || (!is.null(previous) && !length(previous))) {
        out <- previous
      } else {
        args <- c(list(target_gene = gene, tissue_priors = priors, project_dir = project_dir,
                       previous_result = previous,
                       temp_dir = file.path(iteration_dir, "temp_plink", sprintf("chunk_%03d", chunk))),
                  settings[setdiff(names(settings), c("schema", "package_version", "cohort_source"))])
        out <- tryCatch(do.call(run_gene, args), error = function(e) {
          if (inherits(e, "em_fit_error") || iteration > 0L) stop(e)
          list(gene = gene, error = conditionMessage(e))
        })
      }
      validate(out, previous, converged = FALSE)
      attr(out, "em_iteration") <- iteration
      attr(out, "fit_max_iter") <- settings$max_iter
      em_atomic_write(out, file)
    }
    summary$status[i] <- if (is.null(out$error)) "finished" else "gene_error"
    summary$n_tissue_errors[i] <- length(attr(out, "tissue_errors", exact = TRUE))
    if (is.null(out$error)) summary$n_nonconverged[i] <- sum(vapply(out,
      function(x) !isTRUE(x$fit_slide_prior$converged), logical(1)))
  }
  em_atomic_write(summary, file.path(iteration_dir, "completed", sprintf("chunk_%03d.csv", chunk)), csv = TRUE)
  if (any(summary$n_nonconverged > 0))
    stop("Unconverged slider fits; no completion marker written. Inspect and resume this iteration.")
  em_atomic_write(list(iteration = iteration, chunk = chunk, n_genes = length(genes),
                       elapsed_seconds = as.numeric(difftime(Sys.time(), started_at, units = "secs"))), marker)
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: run_slide_prior_em_chunk.R PROJECT_DIR ITERATION_DIR CHUNK")
  for (file in c("em_utils.R", "slide_prior_em_utils.R", "workhorse_slide_prior_em.R"))
    source(file.path(args[1], "script/scan_tissue_attempt", file))
  version <- spe_check_packages()
  saved <- readRDS(file.path(args[2], "settings.rds"))
  if (!identical(version, saved$package_version) ||
      !identical(spe_settings(), saved[names(spe_settings())]))
    stop("Slider settings/package version changed since preparation.")
  spe_run_chunk(args[1], args[2], as.integer(args[3]))
}
