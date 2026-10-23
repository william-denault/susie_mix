run_loading_checks <- function() {
  root <- normalizePath(getwd(), winslash = "/")
  script <- file.path(root, "script/sim/sim_additive_slide_init.R")
  expressions <- parse(script)
  previous_root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
  previous_wd <- getwd()
  on.exit({
    setwd(previous_wd)
    if (is.na(previous_root)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else
      Sys.setenv(SUSIE_MIX_PROJECT_DIR = previous_root)
  })

  # source() discovers the project even when the working directory is elsewhere.
  Sys.unsetenv("SUSIE_MIX_PROJECT_DIR")
  setwd(dirname(root))
  sourced <- new.env(parent = globalenv())
  source(script, local = sourced)
  stopifnot(identical(sourced$init_project_dir, root),
    is.function(sourced$qc_filter_geno), is.function(sourced$cs_summary),
    is.function(sourced$run_additive_initialization))

  # Simulate direct Console evaluation: no ofile or --file, at top level.
  setwd(root)
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
  console <- new.env(parent = globalenv())
  console$interactive <- function() TRUE
  console$sys.nframe <- function() 0L
  console$sys.frames <- function() list()
  console$commandArgs <- function(trailingOnly = FALSE) character()
  eval(expressions[[1L]], console)
  # Nested source() calls have a real filename even when the outer code was pasted.
  rm("sys.frames", envir = console)
  for (expr in head(expressions, -1L)[-1L]) eval(expr, console)
  stopifnot(identical(console$init_project_dir, root),
    is.function(console$qc_filter_geno), is.function(console$summarize_metrics))
  console$run_additive_initialization <- function(...) stop("Unexpected automatic batch")
  console$summarize_additive_initialization <- function(...) stop("Unexpected automatic summary")
  eval(tail(expressions, 1L)[[1L]], console)

  # The current project directory is also a valid Console fallback.
  Sys.unsetenv("SUSIE_MIX_PROJECT_DIR")
  console$sys.frames <- function() list()
  eval(expressions[[1L]], console)
  stopifnot(identical(console$init_project_dir, root))

  # Preserve CLI dispatch without launching any real simulations.
  console$interactive <- function() FALSE
  console$commandArgs <- function(trailingOnly = FALSE) c("2", "3", "0.025")
  called <- NULL
  console$run_additive_initialization <- function(chunk, reps_per_chunk, pve)
    called <<- c(chunk, reps_per_chunk, pve)
  eval(tail(expressions, 1L)[[1L]], console)
  stopifnot(identical(called, c(2, 3, .025)))
  console$commandArgs <- function(trailingOnly = FALSE) "summarize"
  console$summarize_additive_experiment <- function() called <<- "summary"
  eval(tail(expressions, 1L)[[1L]], console)
  stopifnot(identical(called, "summary"))
  console$commandArgs <- function(trailingOnly = FALSE) "experiment"
  console$run_additive_initialization_experiment <- function() called <<- "experiment"
  eval(tail(expressions, 1L)[[1L]], console)
  stopifnot(identical(called, "experiment"))

  # Missing dependencies fail with actual paths, instead of NA-based paths.
  missing_project <- tempfile("additive_loading_validation_")
  dir.create(missing_project)
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = missing_project)
  failure <- tryCatch({source(script, local = new.env()); NULL}, error = conditionMessage)
  stopifnot(grepl("Cannot read simulation helper file", failure, fixed = TRUE),
    grepl("script/sim/simulation_metric_helpers.R", failure, fixed = TRUE),
    !grepl("NA/", failure, fixed = TRUE))
  cat("PASS: source from another directory, pasted Console setup, project fallback, no interactive auto-run, CLI batch/summary dispatch, and missing-helper errors.\n")
}
run_loading_checks()
