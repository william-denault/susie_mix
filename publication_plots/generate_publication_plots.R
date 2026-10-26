# Open this file in RStudio and click Source (or run source() on its full path).
# Required for Figure 1: install.packages(c("ggplot2", "patchwork"))
# Rebuild summaries from the synced n=250 checkpoints when they change, then
# generate every publication panel and Figure 1. No models are fit.

# Locate this script, so Source works even when RStudio's working directory differs.
.publication_script_file <- local({
  files <- Filter(function(x) !is.null(x) && length(x) == 1L,
                  lapply(sys.frames(), function(frame) frame$ofile))
  path <- if (length(files)) tail(files, 1L)[[1L]] else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg)) sub("^--file=", "", arg[1]) else NA_character_
  }
  # source(..., chdir=TRUE) may retain a relative ofile after changing directory.
  if (!is.na(path) && !file.exists(path) && file.exists(basename(path)))
    path <- basename(path)
  if (!is.na(path)) normalizePath(path, winslash = "/", mustWork = TRUE) else path
})

refresh_publication_summaries <- function(project_dir, results_dir, input_dir,
                                           force = FALSE) {
  repository_dir <- dirname(project_dir)
  chunk_dir <- file.path(results_dir, "chunks")
  files <- sort(list.files(chunk_dir, pattern = "\\.RData$", full.names = TRUE))
  if (!length(files)) stop("No synced checkpoints found in: ", chunk_dir)
  info <- file.info(files)
  signature <- data.frame(file = basename(files), bytes = info$size,
                          modified = as.numeric(info$mtime))
  source_paths <- file.path(repository_dir, "script/sim", c(
    "plot_simulations.R", "simulation_design.R", "simulation_plot_helpers.R",
    "simulation_metric_helpers.R"))
  if (any(!file.exists(source_paths))) stop("Cannot find the simulation summary scripts in: ", repository_dir)
  source_hashes <- tools::md5sum(source_paths)
  summary_paths <- file.path(input_dir, c("metric_summary.csv", "roc_counts.rds",
    "pip_calibration_seed_counts.rds", "file_audit.csv", "configuration_counts.csv"))
  manifest_path <- file.path(input_dir, "publication_summary_manifest.rds")
  manifest <- if (file.exists(manifest_path))
    tryCatch(readRDS(manifest_path), error = function(e) NULL) else NULL
  current <- !force && !is.null(manifest) && all(file.exists(summary_paths)) &&
    identical(manifest$checkpoints, signature) &&
    identical(manifest$sources, source_hashes) &&
    identical(manifest$summaries, tools::md5sum(summary_paths))
  if (current) {
    message("Using verified summaries for ", length(files), " unchanged checkpoints.")
  } else {
    message("Rebuilding summaries from ", length(files), " synced checkpoints; no model fits.")
    previous_root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
    previous_options <- options(
      susie.sim.chunk_dir = chunk_dir, susie.sim.summary_dir = input_dir,
      susie.sim.methods = c("SuSiE", "SuSiE-mix", "SuSiE-slide-prior"),
      susie.sim.reuse_saved_summaries = FALSE, susie.sim.summaries_only = TRUE)
    on.exit({
      options(previous_options)
      if (is.na(previous_root)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else
        Sys.setenv(SUSIE_MIX_PROJECT_DIR = previous_root)
    }, add = TRUE)
    Sys.setenv(SUSIE_MIX_PROJECT_DIR = repository_dir)
    source(source_paths[1], local = new.env(parent = environment()))
    # Do not cache summaries if a sync changed the input files during the read.
    after_files <- sort(list.files(chunk_dir, pattern = "\\.RData$", full.names = TRUE))
    after <- file.info(after_files)
    after_signature <- data.frame(file = basename(after_files), bytes = after$size,
                                   modified = as.numeric(after$mtime))
    if (!identical(signature, after_signature))
      stop("Checkpoints changed during summarization. Finish syncing and run again.")
    saveRDS(list(checkpoints = signature, sources = source_hashes,
                 summaries = tools::md5sum(summary_paths)), manifest_path)
  }
  audit <- read.csv(file.path(input_dir, "file_audit.csv"), stringsAsFactors = FALSE)
  message("Checkpoint audit: ", nrow(audit), " files; ", sum(audit$included),
          " included replicates; ", sum(audit$errors), " error records; ",
          sum(audit$saved < 300L), " files with fewer than 300 recorded attempts.")
  if (any(audit$errors > 0L) || any(audit$saved < 300L))
    message("Review file_audit.csv: error records are excluded consistently from all three methods.")
  invisible(audit)
}

generate_publication_figures <- function(
    project_dir = if (!is.na(.publication_script_file))
      dirname(.publication_script_file) else getwd(),
    input_dir = file.path(results_dir, "figures"),
    output_dir = file.path(project_dir, "figures"),
    results_dir = file.path(dirname(project_dir), "simulation results/n250_r300_v1"),
    max_causal = 5L,
    rebuild_summaries = FALSE) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  # Do not inherit stale manuscript autorun settings from another RStudio session.
  previous_options <- options(susie.manuscript.autorun = FALSE)
  on.exit(options(previous_options), add = TRUE)
  source(file.path(project_dir, "script/sim/plot_publication_simulations.R"), local = TRUE)
  source(file.path(project_dir, "script/sim/main_manuscript_fig1.R"), local = TRUE)
  missing <- c("ggplot2", "patchwork")[!vapply(c("ggplot2", "patchwork"),
    requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the missing packages: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))")
  audit <- refresh_publication_summaries(project_dir, results_dir, input_dir,
                                         force = rebuild_summaries)
  panels <- make_publication_simulation_figures(project_dir, input_dir, output_dir,
                                                max_causal = max_causal)
  figure1 <- make_manuscript_simulation_figure(project_dir, input_dir,
    output_dir = file.path(output_dir, "figure1"), max_causal = max_causal)
  invisible(list(panels = panels, figure1 = figure1, audit = audit))
}

# Set this option to FALSE only when you want to call the function yourself.
if (isTRUE(getOption("susie.publication.autorun", TRUE))) {
  publication_figures <- generate_publication_figures()
}
