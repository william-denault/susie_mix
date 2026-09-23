# Rscript --vanilla script/analysis/tests/test_slide_source.R
# RStudio Jobs wraps source/eval in additional calls whose frames lack ofile.
helper <- normalizePath("script/analysis/slide_analysis_utils.R", winslash = "/", mustWork = TRUE)
check_loaded <- function(env) {
  for (name in c("generate_summary_results_slide", "summarize_credible_sets",
                 "select_tss_disagreement", "plot_one_cs_lead_distance_overview")) {
    stopifnot(is.function(get(name, envir = env, inherits = FALSE)))
  }
}
run_job <- function(file, env) {
  stopifnot(is.null(get0("ofile", envir = environment(), inherits = FALSE)))
  source(file, local = env)
}
job_env <- new.env(parent = globalenv())
run_job(helper, job_env)
check_loaded(job_env)

# The outer source file lives elsewhere; its directory is not the helper's.
outer_script <- tempfile("slide-source-", fileext = ".R")
writeLines("source(helper, local = TRUE)", outer_script)
outer_env <- new.env(parent = globalenv())
outer_env$helper <- helper
old_wd <- getwd()
tryCatch({
  setwd(tempdir())
  run_job(outer_script, outer_env)
  check_loaded(outer_env)
}, finally = {
  setwd(old_wd)
  unlink(outer_script)
})

direct_env <- new.env(parent = globalenv())
source(helper, local = direct_env)
check_loaded(direct_env)
cat("PASS: slide helpers load directly, inside a Jobs-style wrapper, and through an outer source in another directory.\n")
