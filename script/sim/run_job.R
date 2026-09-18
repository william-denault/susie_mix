# One manifest row per SLURM task. Can also be sourced for a small local run.
.run_sources <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
.run_source <- if (any(nzchar(.run_sources))) tail(.run_sources[nzchar(.run_sources)], 1L) else
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.run_source <- normalizePath(.run_source, winslash = "/", mustWork = TRUE)
source(file.path(dirname(.run_source), "sim_workhorse.R"), local = TRUE)
sim_project_dir <- dirname(dirname(dirname(.run_source)))
rm(.run_source, .run_sources)

run_simulation_job <- function(job_id, project_dir = sim_project_dir,
                               temp_dir = file.path(project_dir, "temp_plink"),
                               manifest_file = file.path(project_dir, "script/sim/jobs_slide/manifest.csv")) {
  manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
  if (length(job_id) != 1L || !is.finite(job_id) || job_id != floor(job_id))
    stop("job_id must be one integer from the manifest.")
  row <- manifest[manifest$job_id == job_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Job id is missing or repeated in the manifest: ", job_id)
  counts <- as.numeric(row[sim_count_columns])
  if (row$schema_version != sim_schema_version || row$delta_prec != -.5 || row$delta_pdom != .5 ||
      row$K != sum(counts) || row$name != sim_scenario_name(counts))
    stop("Unsupported or inconsistent simulation design; regenerate the manifest.")
  if (!dir.exists(temp_dir) || !length(list.files(temp_dir, pattern = "\\.raw$")))
    stop("No genotype .raw files found in ", temp_dir)
  # Fail once on an installation problem instead of saving 400 error records.
  packages <- c("susieR", "susieSlide", "data.table", "matrixStats")
  for (package in packages) if (!requireNamespace(package, quietly = TRUE))
    stop("Install required package before running jobs: ", package)
  versions <- setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), ""), packages)
  if (!"min_obs" %in% names(formals(susieSlide::susie)))
    stop("The installed susieSlide package does not provide the genotype slider.")

  argument_names <- c("pve", "n", "L", sim_count_columns, "all_additive",
                       "min_maf", "hwe_thresh", "min_n_rec", "slide_min_obs")
  args <- as.list(row[argument_names])
  args$temp_dir <- normalizePath(temp_dir, winslash = "/", mustWork = TRUE)
  checkpoint_settings <- list(design = row, arguments = args, package_versions = versions)
  output <- file.path(project_dir, row$output_file)
  results <- list()
  if (file.exists(output)) {
    saved <- new.env(parent = emptyenv())
    load(output, envir = saved)
    if (!identical(saved$checkpoint_settings, checkpoint_settings))
      stop("Checkpoint design, input path or package versions differ; use a new output directory.")
    results <- saved$results
    if (!is.list(results) || length(results) > row$reps_per_chunk)
      stop("Invalid checkpoint results.")
    for (i in seq_along(results)) {
      expected_rep <- (row$chunk - 1L) * row$reps_per_chunk + i
      x <- results[[i]]
      if (!identical(as.numeric(x$seed), as.numeric(row$seed_base + expected_rep)) ||
          !identical(as.numeric(x$replication), as.numeric(expected_rep)) ||
          (is.null(x$error) && (!identical(x$settings$schema_version, sim_schema_version) ||
            !identical(x$metrics$method, c("SuSiE", "SuSiE-mix", "SuSiE-slide")))))
        stop("Incomplete or incompatible result in checkpoint at replicate ", i)
    }
  }
  if (length(results) == row$reps_per_chunk) {
    message("Job ", job_id, " is already complete: ", output)
    return(invisible(output))
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  for (o in seq.int(length(results) + 1L, row$reps_per_chunk)) {
    replication <- (row$chunk - 1L) * row$reps_per_chunk + o
    args$seed <- row$seed_base + replication
    result <- tryCatch(do.call(sim_mix, args),
                       error = function(e) list(error = conditionMessage(e)))
    result$seed <- args$seed
    result$replication <- replication
    results[[o]] <- result
    temporary <- paste0(output, ".tmp")
    save(results, checkpoint_settings, file = temporary)
    if (!file.rename(temporary, output)) stop("Could not replace checkpoint: ", output)
    message("Job ", job_id, ": ", o, "/", row$reps_per_chunk,
            if (!is.null(result$error)) paste0(" ERROR: ", result$error) else "")
  }
  invisible(output)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (!length(arguments)) stop("Usage: Rscript run_job.R JOB_ID [PROJECT_DIR] [GENOTYPE_DIR]")
  root <- if (length(arguments) >= 2L) arguments[2L] else sim_project_dir
  genotypes <- if (length(arguments) >= 3L) arguments[3L] else file.path(root, "temp_plink")
  run_simulation_job(as.numeric(arguments[1L]), root, genotypes)
}
