# RStudio: edit and Source. Generates the original four-panel SNP/PIP and
# expression plots, plus the lead-distance overview. Requires the GTEx inputs.
em_one_cs_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest",
  expected_coding = "dominant",  # Also accepts "recessive".
  association_threshold = 5e-8,
  minimum_mean_reads = 100,
  both_one_cs = FALSE,           # TRUE restricts additive SuSiE to one CS too.
  expression_plots = TRUE,       # FALSE requests only saved-fit summaries/overview.
  datadir = "/project2/mstephens/gtex",
  top_n = 20L
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  expression_flag <- "--expression-plots" %in% args
  summary_only <- "--summary-only" %in% args
  if (expression_flag && summary_only) stop("Choose only one of --expression-plots and --summary-only.")
  args <- setdiff(args, c("--expression-plots", "--summary-only"))
  if (length(args) > 2L) stop("Usage: finding_interesting_dominant_1cs_em.R [PROJECT_DIR] [latest|ITERATION] [--summary-only]")
  if (length(args) >= 1L) em_one_cs_settings$project_dir <- args[1]
  if (length(args) >= 2L) em_one_cs_settings$iteration <- args[2]
  if (expression_flag) em_one_cs_settings$expression_plots <- TRUE
  if (summary_only) em_one_cs_settings$expression_plots <- FALSE
}
source(file.path(em_one_cs_settings$project_dir, "script/analysis/em_analysis_utils.R"))
source(file.path(em_one_cs_settings$project_dir, "script/analysis/em_one_cs_utils.R"))
em_one_cs_results <- do.call(em_plot_one_cs, em_one_cs_settings)
