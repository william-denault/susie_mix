# Plot the saved slide-prior simulations in the same style as slide_v1.
# Run this whole file with source() in R/RStudio, or with Rscript.
# Only base R is needed; this does not run simulations or refit models.
#
# Input:  simulation results/slide_prior_v1/chunks
# Output: simulation results/slide_prior_v1/figures (PDF, PNG and summaries)
# Figures: coverage, purity, power, CS size, ROC, power-FDR and PIP calibration.

local({
  project_dir <- Sys.getenv("SUSIE_MIX_PROJECT_DIR",
    "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix")

  # FALSE reads all available chunks, including newly added results.
  # After a successful full run, TRUE redraws using the saved summaries.
  # Keep FALSE whenever the contents of chunks/ have changed.
  reuse_saved_summaries <- FALSE

  # Set the dataset and all five methods explicitly, even if the R session
  # still has options from plotting the older three-method slide_v1 results.
  old_options <- options(
    susie.sim.results_dir = "simulation results/slide_prior_v1",
    susie.sim.methods = c("SuSiE", "SuSiE-mix", "SuSiE-slide",
                          "SuSiE-slide-prior", "SuSiE-init-slide"),
    susie.sim.reuse_saved_summaries = reuse_saved_summaries
  )
  on.exit(options(old_options), add = TRUE)

  # Reuse the original layouts, colors, metrics and confidence intervals.
  source(file.path(project_dir, "script/sim/plot_simulations.R"), local = TRUE)
})
