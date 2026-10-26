# Open this file in RStudio and click Source (or run source() on its full path).
# Required for Figure 1: install.packages(c("ggplot2", "patchwork"))
# Use the original n=500 publication snapshots, including the finite slide prior.
# Generate all panels and Figure 1 for causal counts 1:4; no models are fit.

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

generate_publication_figures <- function(
    project_dir = if (!is.na(.publication_script_file))
      dirname(.publication_script_file) else getwd(),
    input_dir = file.path(project_dir, "input_data"),
    output_dir = file.path(project_dir, "figures"),
    max_causal = 4L) {
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
  panels <- make_publication_simulation_figures(project_dir, input_dir, output_dir,
                                                max_causal = max_causal)
  figure1 <- make_manuscript_simulation_figure(project_dir, input_dir,
    output_dir = file.path(output_dir, "figure1"), max_causal = max_causal)
  invisible(list(panels = panels, figure1 = figure1))
}

# Set this option to FALSE only when you want to call the function yourself.
if (isTRUE(getOption("susie.publication.autorun", TRUE))) {
  publication_figures <- generate_publication_figures()
}
