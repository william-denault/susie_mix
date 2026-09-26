# RStudio: edit these settings and Source. Uses the same iteration as the EM summary.
em_descriptive_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest",
  association_threshold = 1e-8,
  minimum_mean_reads = 100,
  tss_bandwidth_kb = 10,         # Same Gaussian bandwidth for both models.
  tss_plot_limit_kb = 200
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: descriptive_results_em.R [PROJECT_DIR] [latest|ITERATION]")
  if (length(args) >= 1L) em_descriptive_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_descriptive_settings$iteration <- args[2]
}
source(file.path(em_descriptive_settings$project_dir, "script/analysis/em_analysis_utils.R"))
source(file.path(em_descriptive_settings$project_dir, "script/analysis/descriptive_results_weighted.R"))
em_context <- em_analysis_context(em_descriptive_settings$project_dir, em_descriptive_settings$iteration)
em_summaries <- em_load_summary(em_context)
em_descriptive_tables <- run_weighted_descriptive_results(
  file.path(em_context$summary_dir, "res_summary.RData"),
  file.path(em_context$summary_dir, "res_cs_summary.RData"),
  file.path(em_context$iteration_dir, "descriptive_results"),
  association_threshold = em_descriptive_settings$association_threshold,
  minimum_mean_reads = em_descriptive_settings$minimum_mean_reads,
  project_dir = em_context$project_dir,
  tss_bandwidth_kb = em_descriptive_settings$tss_bandwidth_kb,
  tss_plot_limit_kb = em_descriptive_settings$tss_plot_limit_kb)
