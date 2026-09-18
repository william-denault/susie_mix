# One array task runs the existing gene chunk with the frozen tissue priors.
em_run_chunk <- function(project_dir, iteration_dir, chunk, run_gene = run_susie_gene) {
  priors <- read.csv(file.path(iteration_dir, "priors.csv"), stringsAsFactors = FALSE)
  em_validate_priors(priors)
  iteration <- unique(priors$iteration)
  if (length(iteration) != 1L) stop("Expected exactly one iteration in priors.csv.")
  manifest <- read.csv(file.path(iteration_dir, "manifest.csv"), stringsAsFactors = FALSE)
  genes <- manifest$gene[manifest$chunk == chunk]
  if (!length(genes)) stop("No genes for chunk ", chunk)
  marker <- file.path(iteration_dir, "completed", sprintf("chunk_%03d.done", chunk))
  if (file.exists(marker)) {
    message("Chunk ", chunk, " already completed.")
    return(invisible(NULL))
  }
  started_at <- Sys.time()
  source_iteration <- unique(priors$source_iteration)
  if (length(source_iteration) != 1L) stop("Expected exactly one source iteration.")
  source_dir <- if (source_iteration == 0L) file.path(project_dir, "results") else
    file.path(project_dir, "results_em", sprintf("iteration_%03d", source_iteration), "results")
  count_nonconverged <- function(out) {
    if (!is.list(out) || !is.null(out[["error"]])) return(0L)
    sum(vapply(out, function(tissue) !identical(tissue$weighted_fit_mix$converged, TRUE), logical(1)))
  }
  summary <- data.frame(gene = genes, status = "", n_tissue_errors = 0L, n_nonconverged = 0L)
  for (i in seq_along(genes)) {
    gene <- genes[i]
    file <- file.path(iteration_dir, "results", paste0(gene, ".rds"))
    out <- if (file.exists(file)) tryCatch(readRDS(file), error = function(e) NULL) else NULL
    reusable <- is.list(out) && identical(attr(out, "em_iteration"), iteration) &&
      is.null(out[["error"]]) && !length(attr(out, "tissue_errors", exact = TRUE)) &&
      count_nonconverged(out) == 0L
    if (!reusable) {
      cat(sprintf("[%s] iteration %d, chunk %03d: %s\n", Sys.time(), iteration, chunk, gene))
      # Preparation already checked every source file. A later missing/corrupt
      # source is fatal to this chunk, not a gene-level error to silently omit.
      previous_result <- readRDS(file.path(source_dir, paste0(gene, ".rds")))
      out <- tryCatch(
        run_gene(target_gene = gene, tissue_priors = priors, project_dir = project_dir,
                 previous_result = previous_result,
                 temp_dir = file.path(iteration_dir, "temp_plink", sprintf("chunk_%03d", chunk))),
        error = function(e) {
          if (inherits(e, "em_fit_error")) stop(e)
          list(gene = gene, error = conditionMessage(e))
        }
      )
      attr(out, "em_iteration") <- iteration
      em_atomic_write(out, file)
    }
    summary$status[i] <- if (!is.null(out[["error"]])) "gene_error" else "finished"
    summary$n_tissue_errors[i] <- length(attr(out, "tissue_errors", exact = TRUE))
    summary$n_nonconverged[i] <- count_nonconverged(out)
    if (summary$status[i] == "gene_error") message("ERROR [", gene, "]: ", out[["error"]])
  }
  # A marker means every gene was attempted and saved, including explicit
  # error records (the same convention as the original scan drivers).
  em_atomic_write(summary, file.path(iteration_dir, "completed",
                                    sprintf("chunk_%03d.csv", chunk)), csv = TRUE)
  if (any(summary$n_nonconverged > 0)) {
    stop("Chunk ", chunk, " has unconverged SuSiE fits. No completion marker written; ",
         "increase the fitting budget if necessary and resume before advancing EM.")
  }
  finished_at <- Sys.time()
  em_atomic_write(list(iteration = iteration, chunk = chunk, n_genes = length(genes),
                       started_at = format(started_at, tz = "UTC", usetz = TRUE),
                       finished_at = format(finished_at, tz = "UTC", usetz = TRUE),
                       elapsed_seconds = as.numeric(difftime(finished_at, started_at, units = "secs"))), marker)
  message("Chunk ", chunk, " finished: ", sum(summary$status == "gene_error"),
          " gene errors; ", sum(summary$n_tissue_errors), " tissue errors.")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: run_em_chunk.R PROJECT_DIR ITERATION_DIR CHUNK")
  for (package in c("susieR", "data.table", "matrixStats")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package is missing: ", package)
  }
  source(file.path(args[1], "script/scan_tissue_attempt/em_utils.R"))
  source(file.path(args[1], "script/scan_tissue_attempt/workhorse_em.R"))
  em_run_chunk(args[1], args[2], as.integer(args[3]))
}
