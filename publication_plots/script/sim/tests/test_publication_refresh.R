# Run from the repository root. Check the default refresh uses old saved inputs.
local({
  old_options <- options(susie.publication.autorun = FALSE, susie.manuscript.autorun = TRUE)
  on.exit(options(old_options))
  entry <- normalizePath("publication_plots/generate_publication_plots.R", winslash = "/")
  source("publication_plots/generate_publication_plots.R", local = TRUE, chdir = TRUE)
  stopifnot(identical(.publication_script_file, entry))
  pub <- tempfile("publication-entry-")
  scripts <- file.path(pub, "script/sim")
  dir.create(scripts, recursive = TRUE)
  fixture <- function(name) sprintf(
    '%s <- function(project_dir, input_dir, output_dir, max_causal) list(input = input_dir, output = output_dir, K = max_causal)', name)
  writeLines(fixture("make_publication_simulation_figures"),
             file.path(scripts, "plot_publication_simulations.R"))
  writeLines(fixture("make_manuscript_simulation_figure"),
             file.path(scripts, "main_manuscript_fig1.R"))
  result <- generate_publication_figures(project_dir = pub)
  expected <- file.path(normalizePath(pub, winslash = "/"), "input_data")
  stopifnot(identical(result$panels$input, expected),
            identical(result$figure1$input, expected),
            identical(result$panels$K, 4L), identical(result$figure1$K, 4L),
            isTRUE(getOption("susie.manuscript.autorun")))
  cat("PASS: RStudio Source resolves the entry path; both generators default to snapshots and K=1:4 without reading raw checkpoints.\n")
})
