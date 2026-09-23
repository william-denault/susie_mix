# RStudio: edit and Source. Produces the top-20 distance bars and PIP scatter.
# Set expression_plots = TRUE to also produce four-panel expression examples.
em_one_cs_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest",
  expected_coding = "recessive",
  association_threshold = 5e-8,
  minimum_mean_reads = 100,
  both_one_cs = FALSE,
  expression_plots = FALSE,
  datadir = "/project2/mstephens/gtex",
  top_n = 20L
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  expression_flag <- "--expression-plots" %in% args
  args <- setdiff(args, "--expression-plots")
  if (length(args) > 2L) stop("Usage: finding_interesting_recessive_1cs_em.R [PROJECT_DIR] [latest|ITERATION] [--expression-plots]")
  if (length(args) >= 1L) em_one_cs_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_one_cs_settings$iteration <- args[2]
  if (expression_flag) em_one_cs_settings$expression_plots <- TRUE
}
source(file.path(em_one_cs_settings$project_dir, "script/analysis/em_analysis_utils.R"))
source(file.path(em_one_cs_settings$project_dir, "script/analysis/em_one_cs_utils.R"))
em_one_cs_results <- do.call(em_plot_one_cs, em_one_cs_settings)
