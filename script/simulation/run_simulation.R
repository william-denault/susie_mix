#!/usr/bin/env Rscript
# Run from any directory. See README.md for local and Slurm commands.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "../.."),
                      winslash = "/", mustWork = TRUE)
source(file.path(root, "script/simulation/config.R"))
source(file.path(root, "script/simulation/simulation_utils.R"))

main <- function(args = commandArgs(TRUE)) {
  config <- simulation_config(root)
  flags <- c("pilot", "dry-run", "retry-failed", "save-fits")
  values <- c("config", "output", "replications", "scenario", "replicate", "shard-id",
              "n-shards", "raw-file", "seed")
  opts <- list()
  while (length(args)) {
    key <- sub("^--", "", args[1])
    if (!startsWith(args[1], "--") || !key %in% c(flags, values))
      stop("Unknown argument: ", args[1], ". See script/simulation/README.md.")
    if (key %in% flags) {
      opts[[key]] <- TRUE
      args <- args[-1]
    } else {
      if (length(args) < 2L || startsWith(args[2], "--")) stop("Missing value for --", key)
      opts[[key]] <- args[2]
      args <- args[-c(1, 2)]
    }
  }
  if (!is.null(opts$config)) {
    env <- new.env(parent = globalenv())
    env$config <- config
    sys.source(opts$config, env)
    config <- env$config
  }
  if (isTRUE(opts$pilot)) config$pilot <- TRUE
  if (isTRUE(opts[["save-fits"]])) config$save_fits <- TRUE
  if (!is.null(opts[["raw-file"]])) config$raw_file <- normalizePath(opts[["raw-file"]],
                                                               winslash = "/", mustWork = TRUE)
  if (!is.null(opts$replications)) config$replications <- as.integer(opts$replications)
  if (!is.null(opts$seed)) config$seed <- as.integer(opts$seed)
  if (config$pilot && is.null(opts$output))
    config$output_dir <- file.path(root, "simulation results", "pilot")
  if (!is.null(opts$output)) config$output_dir <- opts$output
  simulation_validate(config)
  scenarios <- simulation_scenarios(config$pve)
  selected <- scenarios
  parse_indices <- function(text, maximum) {
    x <- suppressWarnings(as.numeric(strsplit(text, ",", fixed = TRUE)[[1]]))
    if (!length(x) || any(!is.finite(x) | x != trunc(x) | x < 1 | x > maximum))
      stop("Indices must be integers between 1 and ", maximum)
    unique(as.integer(x))
  }
  if (!is.null(opts$scenario))
    selected <- scenarios[parse_indices(opts$scenario, nrow(scenarios)), , drop = FALSE]
  reps <- seq_len(config$replications)
  if (!is.null(opts$replicate)) reps <- parse_indices(opts$replicate, config$replications)
  shard <- if (is.null(opts[["shard-id"]])) 1L else as.integer(opts[["shard-id"]])
  n_shards <- if (is.null(opts[["n-shards"]])) 1L else as.integer(opts[["n-shards"]])
  if (anyNA(c(shard, n_shards)) || shard < 1 || n_shards < 1 || shard > n_shards)
    stop("Require 1 <= shard-id <= n-shards.")
  jobs <- expand.grid(scenario_id = selected$scenario_id, replication = reps)
  jobs$trial_id <- (jobs$replication - 1L) * nrow(scenarios) + jobs$scenario_id
  jobs <- jobs[(jobs$trial_id - 1L) %% n_shards == shard - 1L, , drop = FALSE]
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  config$output_dir <- normalizePath(config$output_dir, winslash = "/", mustWork = TRUE)
  message(sprintf("%s: %d scenarios x %d replications = %d datasets / %d fits; N=200, L=10.",
    if (config$pilot) "PILOT" else "PRODUCTION", nrow(scenarios), config$replications,
    nrow(scenarios) * config$replications, 2L * nrow(scenarios) * config$replications))
  message(nrow(jobs), " trials selected for shard ", shard, "/", n_shards)
  plan_dir <- file.path(config$output_dir, "plans")
  dir.create(plan_dir, showWarnings = FALSE)
  write.csv(transform(scenarios, replications = config$replications, n = config$n, L = config$L),
    file.path(plan_dir, sprintf("scenarios_shard_%04d.csv", shard)), row.names = FALSE)
  if (isTRUE(opts[["dry-run"]])) {
    message("Plan written. No genotypes extracted or models fitted.")
    return(invisible(NULL))
  }
  if (!nrow(jobs)) stop("No trials selected.")
  for (p in c("data.table", "matrixStats", "susieR", "R.utils")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Install the R package ", p)
  }
  helpers <- simulation_helpers(root)
  catalogue <- fixed_qc <- NULL
  if (is.null(config$raw_file)) {
    catalogue <- simulation_catalogue(config, helpers)
  } else {
    fixed_qc <- simulation_global_qc(simulation_read_raw(config$raw_file), config,
                                      helpers, "fixed_raw_pilot")
    if (is.null(fixed_qc)) stop("No variants survive QC in the pilot raw file.")
  }
  inputs <- if (is.null(config$raw_file)) c(config$gene_list, config$gtf_file,
    paste0(config$genotype_prefix, c(".bed", ".bim", ".fam"))) else config$raw_file
  scientific_config <- config
  scientific_config$output_dir <- NULL
  source_paths <- file.path(root, c("script/simulation/config.R",
    "script/simulation/simulation_utils.R", "script/simulation/run_simulation.R",
    "script/scan_tissue_attempt/workhorse.R",
    "script/scan_tissue_attempt/workhorse_utils.R",
    "script/scan_tissue_attempt/get_gene_annotations.R"))
  manifest <- list(config = scientific_config, scenarios = scenarios,
    source_md5 = tools::md5sum(source_paths),
    input_metadata = file.info(inputs)[, c("size", "mtime"), drop = FALSE],
    input_small_md5 = tools::md5sum(inputs[file.info(inputs)$size <= 1e7]),
    R = R.version.string,
    package_versions = vapply(c("susieR", "data.table", "matrixStats", "R.utils"),
      function(p) as.character(utils::packageVersion(p)), character(1)))
  stamp <- tempfile()
  saveRDS(manifest, stamp)
  signature <- unname(tools::md5sum(stamp))
  unlink(stamp)
  # Separate shards have independent manifests; checkpoints enforce compatibility.
  manifest$signature <- signature
  manifest$session_info <- sessionInfo()
  manifest_file <- file.path(config$output_dir, sprintf("manifest_shard_%04d.rds", shard))
  if (file.exists(manifest_file) && readRDS(manifest_file)$signature != signature)
    stop("Configuration/source/data changed. Use a new output directory.")
  if (!file.exists(manifest_file)) simulation_save(manifest, manifest_file)
  failed <- 0L
  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, ]
    scenario <- scenarios[job$scenario_id, , drop = FALSE]
    path <- file.path(config$output_dir, "replicates", sprintf("scenario_%03d", job$scenario_id),
                      sprintf("rep_%04d.rds", job$replication))
    if (file.exists(path)) {
      previous <- readRDS(path)
      if (!identical(previous$signature, signature))
        stop("Incompatible checkpoint: ", path, ". Use a new output directory.")
      if (previous$status == "ok") next
      if (!isTRUE(opts[["retry-failed"]])) {
        failed <- failed + 1L
        next
      }
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    lock <- paste0(path, ".lock")
    if (!dir.create(lock, showWarnings = FALSE))
      stop("Trial locked by another worker: ", lock, ". Remove stale lock only after worker has stopped.")
    message(sprintf("[%d/%d] scenario %d (%s, K=%d, PVE=%g), replication %d",
      i, nrow(jobs), job$scenario_id, scenario$architecture, scenario$k, scenario$pve, job$replication))
    seed <- simulation_seed(config$seed, job$scenario_id, job$replication)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    record <- list(signature = signature, scenario = scenario,
                   replication = job$replication, seed = seed, pilot = config$pilot)
    tryCatch({
      result <- tryCatch({
        dataset <- simulation_dataset(scenario, config, helpers, catalogue, fixed_qc)
        fits <- simulation_fit(dataset, config)
        list(status = if (all(vapply(fits, function(x) isTRUE(x$converged), logical(1))))
                        "ok" else "fit_failure",
          gene = dataset$gene, sample_ids = dataset$design$ids,
          qc = dataset$qc, rejected_loci = dataset$rejected,
          phenotype = dataset$phenotype, methods = fits)
      }, error = function(e) list(status = "error", error = conditionMessage(e),
                                  rejected_loci = e$rejected_loci))
      record <- c(record, result)
      simulation_save(record, path)
    }, finally = unlink(lock, recursive = TRUE))
    if (record$status != "ok") {
      failed <- failed + 1L
      message("Trial incomplete: ", record$status, " ", record$error)
    }
  }
  if (failed) stop(failed, " selected trials failed or did not converge. Check checkpoints; ",
                   "--retry-failed reruns the same seeds/data, without replacing difficult loci.")
  message("All selected trials complete. Run summarize_simulation.R to aggregate and plot.")
}

main()
