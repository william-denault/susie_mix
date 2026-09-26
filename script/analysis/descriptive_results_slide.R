# RStudio: edit settings and Source. CLI: Rscript this_file.R [PROJECT_DIR]
slide_descriptive_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  association_threshold = 1e-8,
  minimum_mean_reads = 100,
  tss_bandwidth_kb = 10,         # Same Gaussian bandwidth for both models.
  tss_plot_limit_kb = 200
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: descriptive_results_slide.R [PROJECT_DIR]")
  if (length(args)) slide_descriptive_settings$project_dir <- args[1]
}
source(file.path(slide_descriptive_settings$project_dir, "script/analysis/slide_analysis_utils.R"))
source(file.path(slide_descriptive_settings$project_dir, "script/analysis/slide_descriptive_utils.R"))
slide_descriptive_tables <- do.call(run_slide_descriptive_results, slide_descriptive_settings)
