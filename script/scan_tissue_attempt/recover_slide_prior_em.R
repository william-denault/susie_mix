# CLI used by the bounded Slurm coordinator and isolated workers.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: recover_slide_prior_em.R MODE PROJECT ITERATION [ARGS]")
mode <- args[1]; project <- args[2]
for (file in c("em_utils.R", "slide_prior_em_utils.R", "slide_prior_recovery.R"))
  source(file.path(project, "script/scan_tissue_attempt", file))
iteration <- sre_int(args[3], "iteration", 0L)
iteration_dir <- file.path(project, "results_slide_prior_em", sprintf("iteration_%03d", iteration))
if (mode == "new") {
  target <- Sys.getenv("SUSIE_SLIDE_TARGET_ITERATION", "")
  if (nzchar(target)) {
    target <- sre_int(target, "target iteration", 0L)
    history <- read.csv(file.path(project, "results_slide_prior_em/prior_history.csv"))
    if (max(history$iteration) >= target) stop("Target iteration already prepared; resume it instead of creating another.")
  }
  source(file.path(project, "script/scan_tissue_attempt/prepare_slide_prior_em_iteration.R"))
  launch <- spe_prepare_iteration(project, "new")
  ctx <- sre_context(project, launch$iteration_dir, check_package = TRUE)
  sre_audit(ctx)
  cat(ctx$iteration, "\n", sep = "")
} else {
  ctx <- sre_context(project, iteration_dir, check_package = mode %in% c("plan", "gene"),
                      allow_scope_adoption = mode == "exclude-sex-mt")
  if (mode == "exclude-sex-mt") {
    sre_exclude_sex_mt(ctx)
  } else if (mode == "audit") {
    sre_audit(ctx, if (length(args) >= 4L) args[4] else 10L,
              if (length(args) >= 5L) args[5] else 3L)
  } else if (mode == "plan") {
    if (length(args) != 4L) stop("plan requires capacity")
    pilot <- NULL
    if (Sys.getenv("SUSIE_SLIDE_PILOT") == "1" && nzchar(Sys.getenv("SUSIE_SLIDE_PILOT_GENES")))
      pilot <- trimws(readLines(Sys.getenv("SUSIE_SLIDE_PILOT_GENES"), warn = FALSE))
    batch <- sre_plan_batch(ctx, args[4], pilot,
                            automatic = Sys.getenv("SUSIE_SLIDE_AUTOMATIC") == "1",
                            genes_per_task = Sys.getenv("SUSIE_SLIDE_GENES_PER_TASK", "4"))
    cat(batch$status, "\n", sep = "")
    if (batch$status == "batch") cat(batch$path, "\n", batch$tasks, "\n", sep = "")
  } else if (mode == "genes") {
    if (length(args) != 5L) stop("genes requires batch directory and array task")
    cat(paste(sre_task_genes(ctx, args[4], args[5]), collapse = "\n"), "\n", sep = "")
  } else if (mode == "gene") {
    if (length(args) != 6L) stop("gene requires batch directory, array task and gene")
    if (!args[6] %in% sre_task_genes(ctx, args[4], args[5])) stop("Gene is outside this task.")
    source(file.path(project, "script/scan_tissue_attempt/workhorse_slide_prior_em.R"))
    sre_run_gene(ctx, args[6])
  } else if (mode == "verify") {
    sre_verify(ctx)
  } else stop("Unknown recovery mode: ", mode)
}
