# RStudio: paste the frozen job path printed by download_em_coloc_gwas.R and Source.
# CLI: Rscript run_em_gwas_coloc.R PROJECT_DIR JOB_DIR
gwas_coloc_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
    "/project2/mstephens/wdenault/susie_mix"),
  job_dir = if (exists("em_gwas_job_dir")) em_gwas_job_dir else "SET_FROZEN_GWAS_JOB_PATH",
  modes = c("abf", "susie"), # ABF is a separately labelled one-signal sensitivity analysis.
  traits = NULL,
  ld_bfile = NULL,       # GRCh38 PLINK prefix; use ancestry-matched unrelated donors.
  ld_keep_file = NULL,   # FID/IID file selecting the reference donors, if needed.
  ld_ancestry = NULL,    # e.g. "EUR", only after checking donor ancestry.
  ld_build = NULL,       # "GRCh38"; no automatic assumptions about GTEx.
  plink = "/project2/mstephens/gtex/plink2",
  min_mass = .9,
  priors = data.frame(p1 = 1e-4, p2 = 1e-4, p12 = c(1e-6, 5e-6, 1e-5)),
  h4_threshold = .8,
  min_delta = .1,
  effect_prior_sd = list() # Optional named raw-effect-scale SDs, e.g. glucose_2h = ...
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: run_em_gwas_coloc.R PROJECT_DIR JOB_DIR")
  gwas_coloc_settings$project_dir <- args[1]
  gwas_coloc_settings$job_dir <- args[2]
}
source(file.path(gwas_coloc_settings$project_dir, "script/analysis/em_gwas_coloc_utils.R"))
em_gwas_coloc_output <- do.call(run_em_gwas_coloc, gwas_coloc_settings[names(gwas_coloc_settings) != "project_dir"])
