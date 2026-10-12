# RStudio: edit settings and Source. Run screening first. Full GWAS files can be GBs;
# run this on RCC with adequate disk space. Existing cached downloads are verified.
gwas_download_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
    "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest",
  screen_dir = NULL, # NULL selects latest completed screen of the completed iteration.
  job_dir = NULL,    # To resume downloads, paste the printed frozen GWAS job path here.
  cache_dir = NULL,  # NULL -> PROJECT_DIR/gwas_cache. Retains original full summary files.
  chain_grch37 = NULL, # Optional local hg19ToHg38.over.chain.gz; otherwise downloads from UCSC.
  chain_grch36 = NULL, # Optional local hg18ToHg38.over.chain.gz, needed for birth length.
  python = "python3",
  traits = NULL     # NULL -> all 32; or c("ad", "ldl") for a first run.
)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: download_em_coloc_gwas.R [PROJECT_DIR] [latest|ITERATION]")
  if (length(args) >= 1L) gwas_download_settings$project_dir <- args[1]
  if (length(args) >= 2L) gwas_download_settings$iteration <- args[2]
}
local({
  s <- gwas_download_settings
  source(file.path(s$project_dir, "script/analysis/em_analysis_utils.R"))
  source(file.path(s$project_dir, "script/analysis/em_gwas_coloc_utils.R"))
  job <- s$job_dir
  if (is.null(job)) job <- prepare_em_gwas_coloc(s$project_dir, s$iteration, s$screen_dir)
  if (is.null(s$cache_dir)) s$cache_dir <- file.path(s$project_dir, "gwas_cache")
  options <- c(shQuote(file.path(s$project_dir, "script/analysis/download_coloc_gwas.py")),
    "--job-dir", shQuote(job), "--sources", shQuote(file.path(job, "sources.json")),
    "--cache", shQuote(s$cache_dir),
    if (!is.null(s$chain_grch37)) c("--chain-grch37", shQuote(s$chain_grch37)),
    if (!is.null(s$chain_grch36)) c("--chain-grch36", shQuote(s$chain_grch36)),
    if (!is.null(s$traits)) c("--traits", shQuote(paste(s$traits, collapse = ","))))
  status <- system2(s$python, options)
  assign("em_gwas_job_dir", job, envir = .GlobalEnv)
  message("Use this job_dir in run_em_gwas_coloc.R: ", job)
  if (status != 0L) warning("Some downloads failed; inspect gwas/download_status.json. Set job_dir above to resume.")
})
