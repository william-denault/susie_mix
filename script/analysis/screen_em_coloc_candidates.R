# RStudio: edit settings and Source. Uses the latest COMPLETED EM fits, not history max.
# CLI: Rscript screen_em_coloc_candidates.R [PROJECT_DIR] [latest|ITERATION] [GTEX_DIR]
em_coloc_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest",
  datadir = "/project2/mstephens/gtex",
  association_threshold = 5e-8,
  shift_bp = 10000,
  compute_ld = TRUE,
  unshifted_control_fraction = .01, # Reproducible 1% of eligible unshifted gene-tissue pairs.
  export_tsv = TRUE
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 3L) stop("Usage: screen_em_coloc_candidates.R [PROJECT_DIR] [latest|ITERATION] [GTEX_DIR]")
  if (length(args) >= 1L) em_coloc_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_coloc_settings$iteration <- args[2]
  if (length(args) >= 3L) em_coloc_settings$datadir <- args[3]
}
source(file.path(em_coloc_settings$project_dir, "script/analysis/em_analysis_utils.R"))
source(file.path(em_coloc_settings$project_dir, "script/analysis/em_coloc_screen_utils.R"))
em_coloc_screen <- do.call(em_screen_coloc, em_coloc_settings)
