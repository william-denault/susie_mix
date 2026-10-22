# Open this file in RStudio and click Source (or run source() on its full path).
# Required for Figure 1: install.packages(c("ggplot2", "patchwork"))
# All other panels use base R. Existing saved results are read; no models are fit.

# Locate this script, so Source works even when RStudio's working directory differs.
.publication_script_file <- local({
  files <- Filter(function(x) !is.null(x) && length(x) == 1L,
                  lapply(sys.frames(), function(frame) frame$ofile))
  if (length(files)) tail(files, 1L)[[1L]] else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg)) sub("^--file=", "", arg[1]) else NA_character_
  }
})

generate_publication_figures <- function(
    project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR",
      if (!is.na(.publication_script_file))
        dirname(normalizePath(.publication_script_file, winslash = "/")) else getwd()),
    input_dir = file.path(project_dir, "input_data"),
    output_dir = file.path(project_dir, "figures")) {
  # Do not inherit stale manuscript autorun settings from another RStudio session.
  previous_options <- options(susie.manuscript.autorun = FALSE)
  on.exit(options(previous_options), add = TRUE)
  source(file.path(project_dir, "script/sim/plot_publication_simulations.R"), local = TRUE)
  source(file.path(project_dir, "script/sim/main_manuscript_fig1.R"), local = TRUE)
  missing <- c("ggplot2", "patchwork")[!vapply(c("ggplot2", "patchwork"),
    requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the missing packages: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))")
  panels <- make_publication_simulation_figures(project_dir, input_dir, output_dir)
  figure1 <- make_manuscript_simulation_figure(project_dir, input_dir,
    output_dir = file.path(output_dir, "figure1"))
  invisible(list(panels = panels, figure1 = figure1))
}

# Set this option to FALSE only when you want to call the function yourself.
if (isTRUE(getOption("susie.publication.autorun", TRUE))) {
  publication_figures <- generate_publication_figures()
}
