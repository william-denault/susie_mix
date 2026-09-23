# RStudio: edit settings and Source. CLI: Rscript this_file.R [PROJECT_DIR] [endpoints|direction]
slide_one_cs_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  expected_coding = "dominant",
  coding_selection = "endpoints", # "direction" also includes partial dominant effects.
  association_threshold = 5e-8,
  minimum_mean_reads = 100,
  both_one_cs = FALSE,            # TRUE also requires one additive CS.
  top_n = 20L
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: finding_interesting_dominant_1cs_slide.R [PROJECT_DIR] [endpoints|direction]")
  if (length(args) >= 1L) slide_one_cs_settings$project_dir <- args[1]
  if (length(args) >= 2L) slide_one_cs_settings$coding_selection <- args[2]
}
source(file.path(slide_one_cs_settings$project_dir, "script/analysis/slide_analysis_utils.R"))
source(file.path(slide_one_cs_settings$project_dir, "script/analysis/slide_one_cs_utils.R"))
slide_one_cs_results <- do.call(plot_one_cs_slide, slide_one_cs_settings)
