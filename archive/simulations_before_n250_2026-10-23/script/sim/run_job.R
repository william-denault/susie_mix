# Shared runner for generated R jobs or a manifest row; owns save/resume logic.
.run_sources <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
.run_source <- if (any(nzchar(.run_sources))) tail(.run_sources[nzchar(.run_sources)], 1L) else
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.run_source <- normalizePath(.run_source, winslash = "/", mustWork = TRUE)
source(file.path(dirname(.run_source), "sim_workhorse.R"), local = TRUE)
source(file.path(dirname(.run_source), "additive_init_genotypes.R"), local = TRUE)
sim_project_dir <- dirname(dirname(dirname(.run_source)))
rm(.run_source, .run_sources)

# Reuse the same real-gene loader as the initialization experiment. Existing
# .raw pools keep their original sampling; an absent pool selects GTEx genes.
prepare_simulation_genotypes <- function(project_dir = sim_project_dir,
    temp_dir = Sys.getenv("SUSIE_MIX_GENOTYPE_DIR", file.path(project_dir, "temp_plink"))) {
  if (!nzchar(temp_dir)) temp_dir <- file.path(project_dir, "temp_plink")
  prepare_additive_genotypes(genotype_dir = temp_dir, project_dir = project_dir)
}

sim_comparable_checkpoint <- function(settings) {
  # PLINK workspace size changes execution, not the simulated data.
  if (identical(settings$genotype_inputs$mode, "gtex_fallback"))
    settings$genotype_inputs$plink_memory <- NULL
  settings
}

simulate_slide_replicate <- function(arguments, genotype_source) {
  region <- NULL
  if (identical(genotype_source$mode, "raw")) {
    current <- list.files(genotype_source$directory, "\\.raw$", full.names = TRUE)
    if (!identical(normalizePath(current, winslash = "/", mustWork = TRUE), genotype_source$files))
      stop("The simulation genotype file list changed during this run.")
    arguments$temp_dir <- genotype_source$directory
  } else if (identical(genotype_source$mode, "gtex")) {
    # A seed chooses the same gene across scenarios, PVEs and resumed runs.
    set.seed(arguments$seed)
    index <- sample.int(nrow(genotype_source$loci), 1L)
    region <- extract_additive_region(genotype_source, index)
    on.exit(region$cleanup(), add = TRUE)
    arguments$temp_dir <- region$directory
  } else stop("Unknown simulation genotype source mode.")

  # Extraction errors above stop the job, leaving the last successful
  # checkpoint resumable. Preserve the existing handling of simulation errors.
  result <- tryCatch(do.call(sim_mix, arguments),
                     error = function(e) list(error = conditionMessage(e)))
  if (!is.null(region)) {
    result$genotype_mode <- genotype_source$settings$mode
    result$gene <- region$gene
    result$region <- region$region
  }
  result
}

run_simulation_job <- function(job_id, project_dir = sim_project_dir,
                               temp_dir = Sys.getenv("SUSIE_MIX_GENOTYPE_DIR", file.path(project_dir, "temp_plink")),
                               manifest_file = file.path(project_dir, "script/sim/jobs_slide/manifest.csv"),
                               job = NULL) {
  if (length(job_id) != 1L || !is.finite(job_id) || job_id != floor(job_id))
    stop("job_id must be one integer from the manifest.")
  if (is.null(job)) {
    manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
    row <- manifest[manifest$job_id == job_id, , drop = FALSE]
    if (nrow(row) != 1L) stop("Job id is missing or repeated in the manifest: ", job_id)
  } else {
    if (!is.data.frame(job) || nrow(job) != 1L || !identical(as.numeric(job$job_id), as.numeric(job_id)))
      stop("The embedded job settings must describe exactly job ", job_id)
    row <- job
  }
  counts <- as.numeric(row[sim_count_columns])
  if (row$schema_version != sim_schema_version || row$delta_prec != -.5 || row$delta_pdom != .5 ||
      row$K != sum(counts) || row$name != sim_scenario_name(counts))
    stop("Unsupported or inconsistent simulation design; regenerate the manifest.")
  # Fail once on an installation problem instead of saving 400 error records.
  packages <- c("susieR", "susieSlide", "susieRSlidePrior", "data.table", "matrixStats")
  for (package in packages) if (!requireNamespace(package, quietly = TRUE))
    stop("Install required package before running jobs: ", package)
  versions <- setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), ""), packages)
  if (!"min_obs" %in% names(formals(susieSlide::susie)))
    stop("The installed susieSlide package does not provide the genotype slider.")
  if (!all(c("delta_grid", "delta_prior") %in% names(formals(susieRSlidePrior::susie))))
    stop("Update susieRSlidePrior to the finite-prior implementation.")
  init_argument <- sim_initialization_argument()
  genotype_source <- prepare_simulation_genotypes(project_dir, temp_dir)

  argument_names <- c("pve", "n", "L", sim_count_columns, "all_additive",
                       "min_maf", "hwe_thresh", "min_n_rec", "slide_min_obs")
  args <- as.list(row[argument_names])
  if (identical(genotype_source$mode, "raw")) args$temp_dir <- genotype_source$directory
  checkpoint_settings <- list(design = row, arguments = args, package_versions = versions,
    methods = sim_methods, delta_grid = sim_delta_grid, delta_prior = sim_delta_prior,
    init_argument = init_argument)
  # Keep cached-pool checkpoints compatible. GTEx runs additionally record
  # their fixed source files and eligible gene pool, never a temporary path.
  if (identical(genotype_source$mode, "gtex"))
    checkpoint_settings$genotype_inputs <- genotype_source$settings
  output <- file.path(project_dir, row$output_file)
  results <- list()
  if (file.exists(output)) {
    saved <- new.env(parent = emptyenv())
    load(output, envir = saved)
    if (!identical(sim_comparable_checkpoint(saved$checkpoint_settings),
                   sim_comparable_checkpoint(checkpoint_settings)))
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
            !identical(x$metrics$method, sim_methods) ||
            !identical(x$settings$slide_prior_grid, sim_delta_grid) ||
            !identical(x$settings$slide_prior_weights, sim_delta_prior) ||
            !identical(x$settings$init_argument, init_argument) ||
            any(vapply(sim_pip_fields, function(field)
              length(x[[field]]) != length(x$susie_pip) || !length(x[[field]]) ||
                any(!is.finite(x[[field]])), logical(1))))))
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
    result <- simulate_slide_replicate(args, genotype_source)
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
  genotypes <- if (length(arguments) >= 3L) arguments[3L] else
    Sys.getenv("SUSIE_MIX_GENOTYPE_DIR", file.path(root, "temp_plink"))
  run_simulation_job(as.numeric(arguments[1L]), root, genotypes)
}
