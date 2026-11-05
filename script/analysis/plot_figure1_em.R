# Figure 1: additive SuSiE versus the latest COMPLETED EM SuSiE-mix iteration.
# RStudio on RCC: edit this block if needed, then click Source.
# No models are refitted. See README_figure1_em.md for inputs and interpretation.
figure1_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  datadir = Sys.getenv("GTEX_DATA_DIR", "/project2/mstephens/gtex"),
  iteration = "latest",             # Resolved once; use 18L to pin iteration 18.
  association_threshold = 1e-8,       # Minimum marginal additive P, strictly below.
  minimum_mean_reads = 100,          # Inclusive; at least one CS in either model.
  example_gene = "CYP4F11",
  example_tissue = "Adrenal Gland",
  example_coding = "recessive",
  output_dir = NULL,                # New run folder inside iteration_NNN/descriptive_results/figure1.
  plot_data_file = NULL,            # Set to an exported figure1_plot_data.rds to redraw without GTEx.
  show_plot = interactive()
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: plot_figure1_em.R [PROJECT_DIR] [latest|ITERATION]")
  if (length(args) >= 1L) figure1_settings$project_dir <- args[1L]
  if (length(args) >= 2L) figure1_settings$iteration <- args[2L]
}
source(file.path(figure1_settings$project_dir, "script/analysis/em_analysis_utils.R"))
source(file.path(figure1_settings$project_dir, "script/analysis/em_cs_agreement_utils.R"))
source(file.path(figure1_settings$project_dir, "script/analysis/figure1_em_utils.R"))
figure1_result <- run_figure1_em(figure1_settings)
