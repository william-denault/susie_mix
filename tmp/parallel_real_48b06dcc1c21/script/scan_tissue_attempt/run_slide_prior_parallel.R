args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: run_slide_prior_parallel.R MODE PROJECT ITERATION [ARGS]")
mode <- args[1]; project <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
for (f in c("em_utils.R", "slide_prior_em_utils.R", "slide_prior_recovery.R",
            "prepare_slide_prior_em_iteration.R", "slide_prior_parallel.R", "slide_prior_run_control.R"))
  source(file.path(project, "script/scan_tissue_attempt", f))
history <- read.csv(file.path(project, "results_slide_prior_em/prior_history.csv"))
latest <- max(history$iteration)
iteration <- if (args[3] == "latest") latest else sre_int(args[3], "iteration", 0L)
if (iteration != latest) stop("Resume the latest prepared iteration (", latest, "); earlier priors are frozen.")
ctx <- sre_context(project, file.path(project, "results_slide_prior_em", sprintf("iteration_%03d", iteration)),
                   check_package = mode %in% c("plan", "preflight", "gene"))
if (mode == "set-target") {
  if (length(args) != 4L) stop("set-target requires the new absolute target iteration")
  spt_set_target(project, iteration, args[4])
} else if (mode == "audit") {
  if (length(args) != 3L) stop("audit takes PROJECT latest")
  spt_audit(project)
} else if (mode == "preflight") {
  for (f in c("workhorse_slide_prior_em.R", "workhorse_utils.R", "get_gene_annotations.R"))
    if (!file.exists(file.path(project, "script/scan_tissue_attempt", f))) stop("Missing worker source: ", f)
  message("Preflight passed: iteration ", iteration, ", ", nrow(ctx$manifest), " genes; frozen package/settings match.")
} else if (mode == "plan") {
  if (length(args) != 5L) stop("plan requires capacity and queue snapshot")
  target <- spt_target(project, iteration, Sys.getenv("SUSIE_SLIDE_TARGET_ITERATION", as.character(iteration)))
  message("Run target iteration: ", target)
  spp_init(ctx)
  plan <- spp_plan(ctx, args[4], args[5])
  if (plan$status == "complete" && iteration < target) {
    launch <- spe_prepare_iteration(project, "new", completion_reader = spp_require_complete)
    ctx <- sre_context(project, launch$iteration_dir, check_package = TRUE)
    spp_init(ctx)
    # No old worker remains when completion was certified. Its old registrations
    # belong to the old iteration and are not adopted by the new iteration.
    plan <- spp_plan(ctx, args[4], args[5])
  }
  cat(plan$status, "\n", plan$iteration, "\n", sep = "")
  if (plan$status == "batch") cat(plan$batch, "\n", plan$tasks, "\n", plan$action, "\n", sep = "")
} else if (mode == "task") {
  if (length(args) != 5L) stop("task requires batch and array index")
  w <- spp_task(ctx, args[4], args[5])
  cat(paste(w$gene, w$action, sep = "\t"), sep = "\n")
  cat("\n")
} else if (mode == "gene") {
  if (length(args) != 6L) stop("gene requires batch, task and gene")
  w <- spp_batch(ctx, args[4]); w <- w[w$task == sre_int(args[5], "task") & w$gene == args[6], , drop = FALSE]
  if (nrow(w) != 1L) stop("Gene is not assigned to this task.")
  if (w$action == "fit") source(file.path(project, "script/scan_tissue_attempt/workhorse_slide_prior_em.R"))
  spp_worker_gene(ctx, w$gene, w$action)
} else stop("Unknown parallel mode.")
