# Estimate priors before any fine mapping, or resume an unfinished iteration.
# CLI: Rscript --vanilla prepare_em_iteration.R PROJECT_DIR [new|resume]
em_prepare_iteration <- function(project_dir, mode = "new", gene_annotations = NULL) {
  if (!mode %in% c("new", "resume")) stop("Mode must be new or resume.")
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  em_dir <- file.path(project_dir, "results_em")
  history_file <- file.path(em_dir, "prior_history.csv")
  history <- previous <- NULL
  latest <- 0L
  if (file.exists(history_file)) {
    history <- read.csv(history_file, stringsAsFactors = FALSE)
    if (!"iteration" %in% names(history) || !nrow(history) ||
        !is.numeric(history$iteration) || anyNA(history$iteration) ||
        any(history$iteration < 1 | history$iteration != floor(history$iteration)) ||
        !identical(sort(unique(as.integer(history$iteration))), seq_len(max(history$iteration)))) {
      stop("Invalid iteration numbers in ", history_file)
    }
    for (i in unique(history$iteration)) em_validate_priors(history[history$iteration == i, ])
    latest <- max(history$iteration)
    previous <- history[history$iteration == latest, , drop = FALSE]
  }
  if (mode == "resume" && latest == 0L) stop("No EM iteration exists to resume.")

  if (latest > 0L) {
    previous_dir <- file.path(em_dir, sprintf("iteration_%03d", latest))
    manifest <- read.csv(file.path(previous_dir, "manifest.csv"), stringsAsFactors = FALSE)
    snapshot <- read.csv(file.path(previous_dir, "priors.csv"), stringsAsFactors = FALSE)
    # Older snapshots predate some diagnostics. Compare every original
    # snapshot column; absent fields in history must stay empty, apart from
    # the legacy PIP-method label added for snapshots without a method tag.
    extra <- setdiff(names(previous), names(snapshot))
    legacy <- !"update_method" %in% names(snapshot)
    extras_valid <- !length(extra) || all(vapply(extra, function(column) {
      values <- previous[[column]]
      all(is.na(values) | (legacy & column == "update_method" & values == "legacy_pip_share"))
    }, logical(1)))
    if (!all(names(snapshot) %in% names(previous)) || !extras_valid ||
        !isTRUE(all.equal(previous[names(snapshot)], snapshot, check.attributes = FALSE))) {
      stop("Latest prior history rows differ from the iteration's priors.csv.")
    }
    pending <- em_pending_chunks(previous_dir, manifest)
    if (mode == "resume") {
      if (!length(pending)) stop("Latest iteration is complete; submit without 'resume' for the next one.")
      return(list(iteration_dir = previous_dir, chunks = pending))
    }
    if (length(pending)) {
      stop("Iteration ", latest, " is unfinished (", length(pending),
           " chunks). Wait for the jobs, or use: sbatch em_susie_mix resume")
    }
    source_dir <- file.path(previous_dir, "results")
    fit_name <- "weighted_fit_mix"
  } else {
    manifest <- em_read_manifest(file.path(project_dir, "data/temp_index"))
    source_dir <- file.path(project_dir, "results")
    fit_name <- "susie_mix"
  }

  # Exclude non-autosomal genes before opening results, so even successful old
  # X/Y/MT fits cannot enter the M-step. Keep chunks stable and audit exclusions.
  manifest <- em_annotate_manifest(manifest, project_dir, gene_annotations)
  included <- em_is_autosome(manifest$chromosome)
  if (!any(included)) stop("No autosomal genes in the EM manifest.")
  files <- em_check_result_files(source_dir, manifest[included, , drop = FALSE],
                                allowed_extra_genes = manifest$gene[!included])
  message("Reading all ", length(files), " autosomal gene results from ", source_dir,
          "; excluding ", sum(!included), " non-autosomal genes.")
  pooled <- em_estimate_coding_priors(files, fit_name, previous)
  excluded <- data.frame(gene = manifest$gene[!included], tissue = rep("", sum(!included)),
                         issue = rep("excluded_chromosome", sum(!included)),
                         message = sprintf("Chromosome %s; analysis restricted to autosomes 1-22",
                                           manifest$chromosome[!included]))
  if (nrow(excluded)) pooled$audit <- rbind(pooled$audit, excluded)
  iteration <- latest + 1L
  priors <- data.frame(iteration = iteration, source_iteration = latest,
                       source_fit = fit_name, pooled$priors,
                       analysis_chromosomes = "autosomes_1_22",
                       n_source_files = length(files),
                       n_source_chromosome_exclusions = sum(!included),
                       n_source_gene_errors = sum(pooled$audit$issue == "gene_error"),
                       n_source_tissue_errors = sum(pooled$audit$issue == "tissue_error"),
                       created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       stringsAsFactors = FALSE)
  iteration_dir <- file.path(em_dir, sprintf("iteration_%03d", iteration))
  if (dir.exists(iteration_dir)) {
    stop("Unrecorded iteration directory already exists: ", iteration_dir,
         ". Inspect it before preparing this iteration again.")
  }
  for (subdir in c("results", "completed", "logs")) {
    dir.create(file.path(iteration_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  em_atomic_write(manifest, file.path(iteration_dir, "manifest.csv"), csv = TRUE)
  em_atomic_write(priors, file.path(iteration_dir, "priors.csv"), csv = TRUE)
  em_atomic_write(pooled$audit, file.path(iteration_dir, "source_audit.csv"), csv = TRUE)
  em_atomic_write(pooled$component_counts, file.path(iteration_dir, "component_counts.rds"))
  em_atomic_write(em_bind_history(history, priors), history_file, csv = TRUE)
  message("Prepared iteration ", iteration, " for ", nrow(priors), " tissues; ",
          nrow(pooled$audit), " source audit entries. Priors saved before job submission.")
  list(iteration_dir = iteration_dir, chunks = sort(unique(manifest$chunk)))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || length(args) > 2L) stop("Usage: prepare_em_iteration.R PROJECT_DIR [new|resume]")
  source(file.path(args[1], "script/scan_tissue_attempt/em_utils.R"))
  launch <- em_prepare_iteration(args[1], if (length(args) == 2L) args[2] else "new")
  # Exactly two stdout lines, consumed by job/em_susie_mix.
  cat(launch$iteration_dir, "\n", paste(launch$chunks, collapse = ","), "\n", sep = "")
}
