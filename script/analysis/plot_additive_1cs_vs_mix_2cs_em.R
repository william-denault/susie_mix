# RStudio: edit these settings and Source. Uses the latest completed EM iteration.
# Run generate_summary_results_em.R first for that iteration.
# Each selected gene/tissue gets additive and EM SNP/PIP panels plus expression
# at each EM CS lead adjusted for the other EM CS lead's predictor coding.
em_add1_mix2_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest",         # Set to 12 to pin iteration 12.
  datadir = "/project2/mstephens/gtex",
  association_threshold = 1e-8,
  minimum_mean_reads = 100
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: plot_additive_1cs_vs_mix_2cs_em.R [PROJECT_DIR] [latest|ITERATION]")
  if (length(args) >= 1L) em_add1_mix2_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_add1_mix2_settings$iteration <- args[2]
}
source(file.path(em_add1_mix2_settings$project_dir, "script/analysis/em_analysis_utils.R"), local = TRUE)
source(file.path(em_add1_mix2_settings$project_dir, "script/analysis/em_add1_mix2_utils.R"), local = TRUE)
em_add1_mix2_results <- do.call(em_plot_add1_mix2, em_add1_mix2_settings)
