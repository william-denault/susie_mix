# Run from the repository root. Uses a tiny fake summary builder, never real fits.
local({
  old_options <- options(susie.publication.autorun = FALSE)
  on.exit(options(old_options))
  entry <- normalizePath("publication_plots/generate_publication_plots.R", winslash = "/")
  source("publication_plots/generate_publication_plots.R", local = TRUE, chdir = TRUE)
  stopifnot(identical(.publication_script_file, entry))
  root <- tempfile("publication-refresh-")
  pub <- file.path(root, "publication_plots")
  results <- file.path(root, "results")
  input <- file.path(results, "figures")
  scripts <- file.path(root, "script/sim")
  dir.create(pub, recursive = TRUE)
  dir.create(file.path(results, "chunks"), recursive = TRUE)
  dir.create(scripts, recursive = TRUE)
  checkpoint <- file.path(results, "chunks", "condition.RData")
  writeLines("checkpoint version one", checkpoint)
  writeLines("unfinished sync", paste0(checkpoint, ".filepart"))
  writeLines("unfinished checkpoint", paste0(checkpoint, ".tmp"))
  for (name in c("simulation_design.R", "simulation_plot_helpers.R", "simulation_metric_helpers.R"))
    writeLines("# fixture", file.path(scripts, name))
  builder <- c(
    'out <- getOption("susie.sim.summary_dir")',
    'stopifnot(isTRUE(getOption("susie.sim.summaries_only")), !getOption("susie.sim.reuse_saved_summaries"))',
    'stopifnot(length(list.files(getOption("susie.sim.chunk_dir"), pattern="[.]RData$")) == 1L)',
    'dir.create(out, recursive=TRUE, showWarnings=FALSE)',
    'counter <- file.path(out, "build_count.rds")',
    'saveRDS(if (file.exists(counter)) readRDS(counter)+1L else 1L, counter)',
    'for (f in c("metric_summary.csv", "configuration_counts.csv")) write.csv(data.frame(value=1), file.path(out,f), row.names=FALSE)',
    'write.csv(data.frame(saved=300L,included=299L,errors=1L), file.path(out,"file_audit.csv"), row.names=FALSE)',
    'for (f in c("roc_counts.rds", "pip_calibration_seed_counts.rds")) saveRDS(list(value=1), file.path(out,f))')
  writeLines(builder, file.path(scripts, "plot_simulations.R"))
  run <- function(force = FALSE) refresh_publication_summaries(pub, results, input, force)
  count <- function() readRDS(file.path(input, "build_count.rds"))
  previous_env <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
  previous_summary_option <- getOption("susie.sim.summary_dir")
  audit <- run()
  stopifnot(count() == 1L, audit$included == 299L, audit$errors == 1L,
            identical(Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset=NA_character_), previous_env),
            identical(getOption("susie.sim.summary_dir"), previous_summary_option))
  run()
  stopifnot(count() == 1L)
  writeLines("checkpoint version two with additional results", checkpoint)
  run()
  stopifnot(count() == 2L)
  writeLines(c(builder, "# changed summary calculation"), file.path(scripts, "plot_simulations.R"))
  run()
  stopifnot(count() == 3L)
  write.csv(data.frame(value=999), file.path(input, "metric_summary.csv"), row.names=FALSE)
  run()
  stopifnot(count() == 4L)
  run(force=TRUE)
  stopifnot(count() == 5L)
  writeLines(c(builder, 'stop("fixture failure")'), file.path(scripts, "plot_simulations.R"))
  stopifnot(inherits(try(run(), silent=TRUE), "try-error"),
            identical(Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset=NA_character_), previous_env),
            identical(getOption("susie.sim.summary_dir"), previous_summary_option))
  cat("PASS: cache reuse, changed checkpoints/scripts/summaries, forced rebuild, temporary-file exclusion, and session restoration.\n")
})
