# RStudio: edit these settings and Source. CLI: Rscript this_file.R [PROJECT_DIR] [latest|12]
em_summary_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest"
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: generate_summary_results_em.R [PROJECT_DIR] [latest|ITERATION]")
  if (length(args) >= 1L) em_summary_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_summary_settings$iteration <- args[2]
}
source(file.path(em_summary_settings$project_dir, "script/analysis/em_analysis_utils.R"))
em_summary_results <- em_generate_summary(em_summary_settings$project_dir, em_summary_settings$iteration)
