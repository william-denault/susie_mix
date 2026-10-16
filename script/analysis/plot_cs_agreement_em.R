# RStudio: edit settings and Source. Compares SuSiE with the selected final EM
# fit, using both EM summary tables. Does not refit models or load GTEx data.
# CLI: Rscript --vanilla script/analysis/plot_cs_agreement_em.R [PROJECT_DIR] [latest|ITERATION]
em_agreement_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest",
  association_threshold = 1e-8, # Same cohort as descriptive_results_em.R.
  minimum_mean_reads = 100,
  minimum_tissue_regions = 100L, # Rank by rates, with this minimum denominator.
  output_dir = NULL,             # Default: selected iteration/descriptive_results/cs_agreement/.
  show_plot = interactive()
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: plot_cs_agreement_em.R [PROJECT_DIR] [latest|ITERATION]")
  if (length(args) >= 1L) em_agreement_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_agreement_settings$iteration <- args[2]
}
source(file.path(em_agreement_settings$project_dir, "script/analysis/em_analysis_utils.R"))
source(file.path(em_agreement_settings$project_dir, "script/analysis/em_cs_agreement_utils.R"))
em_agreement_context <- em_analysis_context(em_agreement_settings$project_dir, em_agreement_settings$iteration)
em_agreement_summaries <- em_load_summary(em_agreement_context)
em_agreement_results <- summarize_em_cs_agreement(
  em_agreement_summaries$res, em_agreement_summaries$cs,
  association_threshold = em_agreement_settings$association_threshold,
  minimum_mean_reads = em_agreement_settings$minimum_mean_reads,
  minimum_tissue_regions = em_agreement_settings$minimum_tissue_regions)
em_agreement_output_dir <- em_agreement_settings$output_dir
if (is.null(em_agreement_output_dir)) em_agreement_output_dir <- file.path(
  em_agreement_context$iteration_dir, "descriptive_results", "cs_agreement")
write_em_cs_agreement(em_agreement_results, em_agreement_output_dir,
                      show_plot = em_agreement_settings$show_plot)
