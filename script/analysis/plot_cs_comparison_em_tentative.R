# Tentative four-panel comparison: additive SuSiE versus converged EM SuSiE-mix.
# RStudio: edit this settings block, then Source. Does not fit any models.
# CLI: Rscript --vanilla script/analysis/plot_cs_comparison_em_tentative.R [PROJECT_DIR] [ITERATION]
# See README_em_cs_tentative.md for panel definitions and input requirements.
em_tentative_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = 18L,                 # Pin the iteration identified as converged.
  association_threshold = 1e-8,
  minimum_mean_reads = 100,
  minimum_tissue_regions = 100L,
  output_dir = NULL,               # iteration_018/descriptive_results/cs_comparison_tentative/
  plot_data_file = NULL,           # Optional exported tentative_plot_data.rds, for local replotting.
  show_plot = interactive()
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: plot_cs_comparison_em_tentative.R [PROJECT_DIR] [ITERATION]")
  if (length(args) >= 1L) em_tentative_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_tentative_settings$iteration <- args[2]
}
source(file.path(em_tentative_settings$project_dir, "script/analysis/em_analysis_utils.R"))
source(file.path(em_tentative_settings$project_dir, "script/analysis/em_cs_agreement_utils.R"))
source(file.path(em_tentative_settings$project_dir, "script/analysis/em_cs_comparison_tentative_utils.R"))
if (is.null(em_tentative_settings$plot_data_file)) {
  em_tentative_context <- em_analysis_context(em_tentative_settings$project_dir,
                                             em_tentative_settings$iteration)
  em_tentative_summaries <- em_load_summary(em_tentative_context)
  em_tentative_result <- summarize_em_cs_tentative(
    em_tentative_summaries$res, em_tentative_summaries$cs, em_tentative_context$priors,
    association_threshold = em_tentative_settings$association_threshold,
    minimum_mean_reads = em_tentative_settings$minimum_mean_reads,
    minimum_tissue_regions = em_tentative_settings$minimum_tissue_regions)
  em_tentative_result$summary_generated_at <- em_tentative_summaries$metadata$generated_at
} else {
  em_tentative_result <- readRDS(em_tentative_settings$plot_data_file)
  if (!identical(em_tentative_result$schema, "em_cs_tentative_v1")) stop("Unsupported plot-data file.")
  if (as.character(em_tentative_settings$iteration) != as.character(em_tentative_result$iteration))
    stop("The requested iteration differs from the exported plot data.")
  message("Replotting exported data; cohort settings come from the export.")
}
em_tentative_output <- em_tentative_settings$output_dir
if (is.null(em_tentative_output)) em_tentative_output <- file.path(
  em_tentative_settings$project_dir, "results_em",
  sprintf("iteration_%03d", em_tentative_result$iteration),
  "descriptive_results", "cs_comparison_tentative")
write_em_cs_tentative(em_tentative_result, em_tentative_output,
                      show_plot = em_tentative_settings$show_plot)
