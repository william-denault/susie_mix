# Runtime controls and a read-only metadata audit. No regional fitting changes.
spt_file <- function(project) file.path(project, "results_slide_prior_em/run_control.rds")
spt_target <- function(project, iteration, fallback = Sys.getenv("SUSIE_SLIDE_TARGET_ITERATION", "20")) {
  path <- spt_file(project)
  if (!file.exists(path)) return(sre_int(fallback, "target iteration", iteration))
  saved <- readRDS(path)
  if (!is.list(saved) || !identical(saved$schema, 1L)) stop("Invalid run_control.rds.")
  sre_int(saved$target_iteration, "persisted target iteration", iteration)
}

spt_set_target <- function(project, iteration, target) {
  target <- sre_int(target, "new target iteration", iteration)
  # Only extend: this command cannot truncate a running experiment.
  old <- spt_target(project, iteration, Sys.getenv("SUSIE_SLIDE_TARGET_ITERATION", as.character(max(20L, iteration))))
  if (target < old) stop("set-target only extends the run; current target is ", old, ".")
  if (target != old || !file.exists(spt_file(project)))
    em_atomic_write(list(schema = 1L, target_iteration = target, previous_target = old,
                        updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), spt_file(project))
  message("Saved target iteration ", spt_target(project, iteration), ". Already prepared: ", iteration,
          ". The next planning check reads this override; do not submit another controller.")
  invisible(target)
}

spt_audit <- function(project, sample_size = 16L) {
  root <- file.path(project, "results_slide_prior_em")
  history <- read.csv(file.path(root, "prior_history.csv"), stringsAsFactors = FALSE)
  if (!nrow(history) || anyNA(history$iteration) || any(history$iteration != floor(history$iteration)) ||
      !identical(sort(unique(as.integer(history$iteration))), 0L:max(history$iteration)))
    stop("Invalid iteration history.")
  latest <- as.integer(max(history$iteration))
  reference <- reference_cohort <- NULL
  for (iteration in 0:latest) {
    d <- file.path(root, sprintf("iteration_%03d", iteration))
    priors <- read.csv(file.path(d, "priors.csv"), stringsAsFactors = FALSE)
    spe_validate_priors(priors)
    if (!sre_same_priors(history[history$iteration == iteration, , drop = FALSE], priors))
      stop("History/snapshot mismatch at ", iteration)
    settings <- readRDS(file.path(d, "settings.rds"))
    if (is.null(reference)) reference <- settings
    if (!identical(settings, reference)) stop("Frozen fitting settings changed at ", iteration)
    if (!identical(spe_settings(), settings[names(spe_settings())])) stop("Current code differs from frozen settings.")
    if (iteration > 0L) {
      prev <- history[history$iteration == iteration - 1L, , drop = FALSE]
      counts <- read.csv(file.path(d, "pooled_counts.csv"), stringsAsFactors = FALSE)
      if (!setequal(counts$tissue, priors$tissue) || anyDuplicated(counts$tissue) ||
          !setequal(prev$tissue, priors$tissue)) stop("Tissue universe changed at ", iteration)
      for (t in priors$tissue) {
        expected <- spe_update_weights(as.numeric(counts[counts$tissue == t, spe_columns()]), spe_prior(prev, t))
        if (max(abs(expected - spe_prior(priors, t))) > 1e-10) stop("M-step count/weight mismatch at ", iteration, "/", t)
      }
      audit <- read.csv(file.path(d, "source_audit.csv"), stringsAsFactors = FALSE)
      cohort <- audit[audit$issue == "fit", c("gene", "tissue"), drop = FALSE]
      cohort <- cohort[order(cohort$gene, cohort$tissue), ]; rownames(cohort) <- NULL
      if (anyDuplicated(cohort) || !nrow(cohort)) stop("Invalid fit cohort at ", iteration)
      if (is.null(reference_cohort)) reference_cohort <- cohort
      if (!identical(cohort, reference_cohort)) stop("Successful fit cohort changed at ", iteration)
      if (sum(priors$n_fits) != nrow(cohort)) stop("Fit count/history mismatch at ", iteration)
      message("Update ", iteration, ": ", nrow(cohort), " gene-tissue fits; normalized pooled counts match all ",
              nrow(priors), " tissue priors. max |weight change| = ", signif(max(priors$prior_max_change), 5))
    }
  }
  ctx <- sre_context(project, file.path(root, sprintf("iteration_%03d", latest)))
  files <- sort(list.files(file.path(spp_root(ctx), "summaries"), "\\.rds$", full.names = TRUE))
  # Spread a bounded, deterministic sample across filenames. Read ONLY compact
  # completed summaries; a running worker can make a sample temporarily stale.
  checked <- stale <- 0L
  if (length(files)) {
    selected <- unique(round(seq(1, length(files), length.out = min(sample_size, length(files)))))
    for (file in files[selected]) {
      gene <- sub("\\.rds$", "", basename(file))
      if (!gene %in% ctx$manifest$gene) stop("Unknown gene summary: ", gene)
      saved <- readRDS(file)
      if (!spp_valid_summary(saved, ctx, gene, spp_stamp(sre_result(ctx, gene)), spp_stamp(spp_previous_file(ctx, gene))) ||
          !sre_ready(saved$row$status)) {
        stale <- stale + 1L
      } else checked <- checked + 1L
    }
  }
  message("Metadata audit passed through prepared iteration ", latest, "; current summary sample: ", checked,
          " valid completed, ", stale, " stale/unfinished. This is not a full posterior audit.")
  message("Run target: ", spt_target(project, latest), ". Frozen tol: ", reference$tol,
          "; grid points: ", length(reference$delta_grid), "; L: ", reference$L, ".")
  invisible(list(iteration = latest, checked = checked, stale = stale))
}
