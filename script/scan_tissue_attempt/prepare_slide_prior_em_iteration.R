# Iteration 000 fits uniform priors; iterations 001+ update and then refit.
spe_prepare_iteration <- function(project_dir, mode = "new", n_chunks = 298L,
                                  package_version = spe_check_packages()) {
  if (!mode %in% c("new", "resume")) stop("Mode must be new or resume.")
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  root <- file.path(project_dir, "results_slide_prior_em")
  history_file <- file.path(root, "prior_history.csv")
  history <- previous <- NULL
  latest <- -1L
  settings <- spe_settings()
  settings$package_version <- package_version
  if (file.exists(history_file)) {
    history <- read.csv(history_file, stringsAsFactors = FALSE)
    if (!nrow(history) || !is.numeric(history$iteration) || anyNA(history$iteration) ||
        any(history$iteration < 0 | history$iteration != floor(history$iteration)) ||
        !identical(sort(unique(as.integer(history$iteration))), 0L:max(history$iteration)))
      stop("Invalid slider iteration history.")
    for (i in unique(history$iteration)) spe_validate_priors(history[history$iteration == i, ])
    latest <- max(history$iteration)
    previous <- history[history$iteration == latest, , drop = FALSE]
  }
  if (mode == "resume" && latest < 0) stop("No slider iteration exists to resume.")
  if (latest >= 0) {
    previous_dir <- file.path(root, sprintf("iteration_%03d", latest))
    snapshot <- read.csv(file.path(previous_dir, "priors.csv"), stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(previous, snapshot, check.attributes = FALSE)))
      stop("Slider history differs from the frozen priors.csv.")
    saved_settings <- readRDS(file.path(previous_dir, "settings.rds"))
    if (!identical(settings, saved_settings[names(settings)]))
      stop("Slider fitting settings/package version changed; inspect before continuing this run.")
    manifest <- read.csv(file.path(previous_dir, "manifest.csv"), stringsAsFactors = FALSE)
    pending <- em_pending_chunks(previous_dir, manifest)
    if (mode == "resume") {
      if (!length(pending)) stop("Latest iteration is complete; use new to continue.")
      return(list(iteration_dir = previous_dir, chunks = pending))
    }
    if (length(pending)) stop("Slider iteration ", latest, " is unfinished. Wait or use resume.")
    files <- em_check_result_files(file.path(previous_dir, "results"), manifest)
    message("Pooling slider counts from ", length(files), " gene files.")
    pooled <- spe_pool(files, previous)
    priors <- pooled$priors
  } else {
    # Reuse only the existing scan's gene universe and tissue labels. Mix
    # weights and posteriors are not valid initializations for this model.
    mix_history <- read.csv(file.path(project_dir, "results_em", "prior_history.csv"),
                            stringsAsFactors = FALSE)
    if (!nrow(mix_history) || !is.numeric(mix_history$iteration) ||
        any(!is.finite(mix_history$iteration))) stop("Cannot identify the mix scan cohort.")
    mix_latest <- max(mix_history$iteration)
    tissues <- mix_history$tissue[mix_history$iteration == mix_latest]
    manifest <- read.csv(file.path(project_dir, "results_em",
                                  sprintf("iteration_%03d", mix_latest), "manifest.csv"),
                          stringsAsFactors = FALSE)
    w <- matrix(1 / 17, length(tissues), 17, dimnames = list(NULL, spe_columns()))
    priors <- data.frame(tissue = tissues, w, n_fits = 0L, n_active = 0L,
                         free_mass = 0, prior_max_change = 0, m_step_q_gain = 0,
                         source_elbo = NA_real_)
    pooled <- list(audit = data.frame(gene = character(), tissue = character(), issue = character()),
                   counts = data.frame(tissue = tissues, w * 0))
    settings$cohort_source <- file.path("results_em", sprintf("iteration_%03d", mix_latest))
    message("Initializing the discrete slider with uniform weights; cohort from ", settings$cohort_source)
  }
  # Keep the initial provenance while comparing only actual fitting settings.
  if (latest >= 0) settings <- readRDS(file.path(previous_dir, "settings.rds"))
  manifest <- em_rechunk_manifest(manifest, n_chunks)
  if (any(grepl("[/\\\\]", manifest$gene)) || any(manifest$gene %in% c(".", "..")))
    stop("Invalid gene filenames in manifest.")
  iteration <- latest + 1L
  priors <- data.frame(iteration = iteration, source_iteration = latest, priors,
                       update_method = "slider_prior_active_counts_v1",
                       created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  spe_validate_priors(priors)
  iteration_dir <- file.path(root, sprintf("iteration_%03d", iteration))
  if (dir.exists(iteration_dir)) stop("Unrecorded iteration directory exists: ", iteration_dir)
  for (subdir in c("results", "completed", "logs"))
    dir.create(file.path(iteration_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  em_atomic_write(manifest, file.path(iteration_dir, "manifest.csv"), csv = TRUE)
  em_atomic_write(priors, file.path(iteration_dir, "priors.csv"), csv = TRUE)
  em_atomic_write(settings, file.path(iteration_dir, "settings.rds"))
  em_atomic_write(data.frame(weight_column = spe_columns(), delta = spe_grid()),
                  file.path(iteration_dir, "grid.csv"), csv = TRUE)
  em_atomic_write(pooled$audit, file.path(iteration_dir, "source_audit.csv"), csv = TRUE)
  em_atomic_write(pooled$counts, file.path(iteration_dir, "pooled_counts.csv"), csv = TRUE)
  em_atomic_write(rbind(history, priors), history_file, csv = TRUE)
  message("Prepared slider iteration ", iteration, " for ", nrow(priors), " tissues and ",
          length(unique(manifest$chunk)), " chunks. Priors are frozen before submission.")
  list(iteration_dir = iteration_dir, chunks = sort(unique(manifest$chunk)))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || length(args) > 2L)
    stop("Usage: prepare_slide_prior_em_iteration.R PROJECT_DIR [new|resume]")
  for (file in c("em_utils.R", "slide_prior_em_utils.R"))
    source(file.path(args[1], "script/scan_tissue_attempt", file))
  launch <- spe_prepare_iteration(args[1], if (length(args) == 2L) args[2] else "new")
  cat(launch$iteration_dir, "\n", sep = "")
  cat(paste(launch$chunks, collapse = ","), "\n", sep = "")
}
