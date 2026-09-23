# Redraw the current saved analysis, recomputing analytic CIs from compact counts.
# Old compact counts without CS-size sums of squares are rebuilt once.
# A first calibration run reads the original checkpoints for exact PIP means;
# subsequent redraws reuse its compact bin counts.
local({
  root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR",
    "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix")
  old <- options(susie.sim.reuse_saved_summaries = TRUE)
  on.exit(options(old))
  source(file.path(root, "script/sim/plot_simulations.R"), local = TRUE)
})
