# source(.../regenerate_figures.R), or run with Rscript from any directory.
# Select the old two-method September panels with:
# options(susie.archive.dataset = "legacy")
local({
  files <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  entry <- if (any(nzchar(files))) tail(files[nzchar(files)], 1L) else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  root <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))
  dataset <- getOption("susie.archive.dataset", "slide_prior_v1")
  stopifnot(dataset %in% c("legacy", "slide_prior_v1"))
  previous_root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
  previous_options <- options(susie.sim.results_dir = if (dataset == "legacy")
    "simulation results" else "simulation results/slide_prior_v1",
    susie.sim.methods = if (dataset == "legacy") c("SuSiE", "SuSiE-mix") else
      c("SuSiE", "SuSiE-mix", "SuSiE-slide", "SuSiE-slide-prior", "SuSiE-init-slide"),
    susie.sim.reuse_saved_summaries = FALSE)
  on.exit({
    options(previous_options)
    if (is.na(previous_root)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else
      Sys.setenv(SUSIE_MIX_PROJECT_DIR = previous_root)
  })
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
  script <- if (dataset == "legacy") "legacy_source/script/sim/plot_simulations.R" else
    "script/sim/plot_simulations.R"
  if (dataset == "legacy") {
    # The historical script hard-coded the original Windows project path.
    # Override only that assignment, preserving the frozen source file.
    for (expression in parse(file.path(root, script))) {
      if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
          identical(expression[[2]], as.name("project_dir"))) expression[[3]] <- root
      eval(expression, envir = environment())
    }
  } else source(file.path(root, script), local = TRUE)
})
