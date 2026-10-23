# Redraw the current saved analysis, recomputing analytic CIs from compact counts.
# Old compact counts without CS-size sums of squares are rebuilt once.
# A first calibration run reads the original checkpoints for exact PIP means;
# subsequent redraws reuse its compact bin counts.
local({
  files <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  entry <- if (any(nzchar(files))) tail(files[nzchar(files)], 1L) else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR",
    dirname(dirname(dirname(normalizePath(entry, winslash = "/", mustWork = TRUE)))))
  old <- options(susie.sim.reuse_saved_summaries = TRUE)
  on.exit(options(old))
  source(file.path(root, "script/sim/plot_simulations.R"), local = TRUE)
})
