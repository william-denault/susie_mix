# RStudio: edit settings and Source. CLI: Rscript this_file.R [PROJECT_DIR]
slide_summary_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  delta_tolerance = 1e-6,
  gtf_file = "/project2/mstephens/gtex/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz"
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: generate_summary_results_slide.R [PROJECT_DIR]")
  if (length(args)) slide_summary_settings$project_dir <- args[1]
}
source(file.path(slide_summary_settings$project_dir, "script/analysis/slide_analysis_utils.R"))
slide_summary_results <- do.call(generate_summary_results_slide, slide_summary_settings)
