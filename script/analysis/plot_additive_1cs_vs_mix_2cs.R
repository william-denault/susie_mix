# plot_additive_1cs_vs_mix_2cs.R
#
# For every significant, sufficiently expressed gene-tissue pair for which
# additive SuSiE has one credible set and SuSiE-mix has two credible sets,
# create a four-panel figure:
#
#   1. additive SuSiE PIPs;
#   2. SuSiE-mix PIPs, with coding-block boundaries;
#   3. expression versus the CS1 lead SNP, adjusted for the CS2 lead predictor;
#   4. expression versus the CS2 lead SNP, adjusted for the CS1 lead predictor.
#
# The phenotype, genotype QC, mixed codings, and tissue-specific predictor
# filtering are reconstructed with the same rules as the current workhorse.
# Each title also reports the one-CS plots' ELBO + KL likelihood comparison
# and min(|additive lead - mixed CS1 lead|, |additive lead - mixed CS2 lead|)
# in GRCh38 base pairs. The audit CSV retains all three pairwise distances.

library(data.table)
library(matrixStats)
library(susieR)


# ============================================================
# User settings
# ============================================================

project_dir <- "/project2/mstephens/wdenault/susie_mix"
datadir <- "/project2/mstephens/gtex"

summary_file <- file.path(project_dir, "res_summary.RData")
results_dir <- file.path(project_dir, "results")
temp_dir <- file.path(project_dir, "temp_plink")

plot_dir <- file.path(
  project_dir,
  "plot",
  "additive_1cs_vs_mix_2cs"
)

plot_summary_file <- file.path(
  plot_dir,
  "additive_1cs_vs_mix_2cs_plot_summary.csv"
)

association_threshold <- 1e-8
minimum_mean_reads <- 100

# These values match scan_tissue_attempt/workhorse.R.
min_maf_plink <- 0.00
min_maf <- 0.05
hwe_thresh <- 1e-8
min_n_rec <- 5
cis_window <- 5e5
min_samples <- 50


# ============================================================
# Shared workhorse-compatible reconstruction functions
# ============================================================

source(
  file.path(
    project_dir,
    "script",
    "analysis",
    "fit_mix_expression_plot_utils.R"
  )
)


# Shared renderer also serves the EM entry point.
source(file.path(project_dir, "script/analysis/add1_mix2_plot_utils.R"))
summary_environment <- new.env(parent = emptyenv())
loaded_objects <- load(summary_file, envir = summary_environment)
if (!"res_summary" %in% loaded_objects) {
  stop("res_summary.RData does not contain an object named res_summary.")
}
cases <- select_add1_mix2_cases(
  summary_environment$res_summary, association_threshold, minimum_mean_reads
)
plot_summary <- run_add1_mix2_plots(
  cases = cases, plot_dir = plot_dir, project_dir = project_dir, datadir = datadir,
  results_dir = results_dir, temp_dir = temp_dir,
  summary_filename = basename(plot_summary_file),
  min_maf_plink = min_maf_plink, min_maf = min_maf, hwe_thresh = hwe_thresh,
  min_n_rec = min_n_rec, cis_window = cis_window, min_samples = min_samples
)
