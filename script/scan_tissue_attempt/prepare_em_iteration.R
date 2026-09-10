# Estimate priors before any fine mapping, or resume an unfinished iteration.
# CLI: Rscript --vanilla prepare_em_iteration.R PROJECT_DIR [new|resume]
em_prepare_iteration <- function(project_dir, mode = "new") {
  if (!mode %in% c("new", "resume")) stop("Mode must be new or resume.")
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  em_dir <- file.path(project_dir, "results_em")
  history_file <- file.path(em_dir, "prior_history.csv")
  history <- previous <- NULL
  latest <- 0L
  if (file.exists(history_file)) {
    history <- read.csv(history_file, stringsAsFactors = FALSE)
    if (!"iteration" %in% names(history) || !nrow(history) ||
        !is.numeric(history$iteration) || anyNA(history$iteration) ||
        any(history$iteration < 1 | history$iteration != floor(history$iteration)) ||
        !identical(sort(unique(as.integer(history$iteration))), seq_len(max(history$iteration)))) {
      stop("Invalid iteration numbers in ", history_file)
    }
    for (i in unique(history$iteration)) em_validate_priors(history[history$iteration == i, ])
    latest <- max(history$iteration)
    previous <- history[history$iteration == latest, , drop = FALSE]
  }
  if (mode == "resume" && latest == 0L) stop("No EM iteration exists to resume.")

  if (latest > 0L) {
    previous_dir <- file.path(em_dir, sprintf("iteration_%03d", latest))
    manifest <- read.csv(file.path(previous_dir, "manifest.csv"), stringsAsFactors = FALSE)
    snapshot <- read.csv(file.path(previous_dir, "priors.csv"), stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(previous, snapshot, check.attributes = FALSE))) {
      stop("Latest prior history rows differ from the iteration's priors.csv.")
    }
    pending <- em_pending_chunks(previous_dir, manifest)
    if (mode == "resume") {
      if (!length(pending)) stop("Latest iteration is complete; submit without 'resume' for the next one.")
      return(list(iteration_dir = previous_dir, chunks = pending))
    }
    if (length(pending)) {
      stop("Iteration ", latest, " is unfinished (", length(pending),
           " chunks). Wait for the jobs, or use: sbatch em_susie_mix resume")
    }
    source_dir <- file.path(previous_dir, "results")
    fit_name <- "weighted_fit_mix"
  } else {
    manifest <- em_read_manifest(file.path(project_dir, "data/temp_index"))
    source_dir <- file.path(project_dir, "results")
    fit_name <- "susie_mix"
  }

  # This checks every expected gene, including explicit failed-fit records.
  files <- em_check_result_files(source_dir, manifest)
  message("Reading all ", length(files), " gene results from ", source_dir)
  pooled <- em_sum_pips(files, fit_name, previous)
  iteration <- latest + 1L
  priors <- data.frame(iteration = iteration, source_iteration = latest,
                       source_fit = fit_name, pooled$priors,
                       n_source_files = length(files),
                       n_source_gene_errors = sum(pooled$audit$issue == "gene_error"),
                       n_source_tissue_errors = sum(pooled$audit$issue == "tissue_error"),
                       created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       stringsAsFactors = FALSE)
  iteration_dir <- file.path(em_dir, sprintf("iteration_%03d", iteration))
  if (dir.exists(iteration_dir)) {
    stop("Unrecorded iteration directory already exists: ", iteration_dir,
         ". Inspect it before preparing this iteration again.")
  }
  for (subdir in c("results", "completed", "logs")) {
    dir.create(file.path(iteration_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  em_atomic_write(manifest, file.path(iteration_dir, "manifest.csv"), csv = TRUE)
  em_atomic_write(priors, file.path(iteration_dir, "priors.csv"), csv = TRUE)
  em_atomic_write(pooled$audit, file.path(iteration_dir, "source_audit.csv"), csv = TRUE)
  em_atomic_write(rbind(history, priors), history_file, csv = TRUE)
  message("Prepared iteration ", iteration, " for ", nrow(priors), " tissues; ",
          nrow(pooled$audit), " source audit entries. Priors saved before job submission.")
  list(iteration_dir = iteration_dir, chunks = sort(unique(manifest$chunk)))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || length(args) > 2L) stop("Usage: prepare_em_iteration.R PROJECT_DIR [new|resume]")
  source(file.path(args[1], "script/scan_tissue_attempt/em_utils.R"))
  launch <- em_prepare_iteration(args[1], if (length(args) == 2L) args[2] else "new")
  # Exactly two stdout lines, consumed by job/em_susie_mix.
  cat(launch$iteration_dir, "\n", paste(launch$chunks, collapse = ","), "\n", sep = "")
}
