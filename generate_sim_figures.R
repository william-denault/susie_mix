# Source this file or run Rscript --vanilla generate_sim_figures.R.
# Reads only the new n=250, 300-replicate, three-method archive.
local({
  files <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  entry <- if (any(nzchar(files))) tail(files[nzchar(files)], 1L) else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  root <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))
  previous_root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
  previous_options <- options(susie.sim.results_dir = "simulation results/n250_r300_v1",
    susie.sim.methods = c("SuSiE", "SuSiE-mix", "SuSiE-slide-prior"),
    susie.sim.reuse_saved_summaries = FALSE)
  on.exit({
    options(previous_options)
    if (is.na(previous_root)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else
      Sys.setenv(SUSIE_MIX_PROJECT_DIR = previous_root)
  })
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
  source(file.path(root, "script/sim/plot_simulations.R"), local = TRUE)
})
