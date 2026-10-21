# Helpers for post-EM analysis. Original fits are used only as comparators;
# weighted fits always come from the selected, fully completed EM iteration.
em_analysis_context <- function(project_dir, iteration = "latest") {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  source(file.path(project_dir, "script/scan_tissue_attempt/em_utils.R"), local = TRUE)
  em_dir <- file.path(project_dir, "results_em")
  history <- read.csv(file.path(em_dir, "prior_history.csv"), stringsAsFactors = FALSE)
  ids <- history$iteration
  if (!is.numeric(ids) || !length(ids) || anyNA(ids) ||
      any(ids < 1 | ids != floor(ids))) stop("Invalid EM prior history iterations.")
  candidates <- sort(unique(ids), decreasing = TRUE)
  latest <- identical(as.character(iteration), "latest")
  if (!latest) {
    requested <- suppressWarnings(as.numeric(iteration))
    if (length(requested) != 1L || !is.finite(requested) ||
        !requested %in% candidates) stop("Requested EM iteration is absent from prior history.")
    candidates <- requested
  }
  rejected <- character()
  for (i in candidates) {
    context <- tryCatch({
      iteration_dir <- file.path(em_dir, sprintf("iteration_%03d", as.integer(i)))
      manifest <- read.csv(file.path(iteration_dir, "manifest.csv"), stringsAsFactors = FALSE)
      if (!all(c("gene", "chunk") %in% names(manifest)) || !nrow(manifest) ||
          anyNA(manifest[c("gene", "chunk")]) || anyDuplicated(manifest$gene) ||
          any(!nzchar(manifest$gene) | grepl("[/\\\\]", manifest$gene)) ||
          !is.numeric(manifest$chunk) || any(manifest$chunk < 1 | manifest$chunk != floor(manifest$chunk))) {
        stop("Invalid iteration manifest.")
      }
      pending <- em_pending_chunks(iteration_dir, manifest)
      if (length(pending)) stop(length(pending), " unfinished chunks")
      for (chunk in unique(manifest$chunk)) {
        marker <- readRDS(file.path(iteration_dir, "completed", sprintf("chunk_%03d.done", chunk)))
        if (!isTRUE(marker$iteration == i) || !isTRUE(marker$chunk == chunk) ||
            !isTRUE(marker$n_genes == sum(manifest$chunk == chunk))) {
          stop("Invalid completion marker for chunk ", chunk)
        }
      }
      result_files <- em_check_result_files(file.path(iteration_dir, "results"), manifest)
      priors <- read.csv(file.path(iteration_dir, "priors.csv"), stringsAsFactors = FALSE)
      em_validate_priors(priors)
      expected <- history[history$iteration == i, , drop = FALSE]
      if (!all(names(priors) %in% names(expected)) ||
          !isTRUE(all.equal(priors, expected[names(priors)], check.attributes = FALSE))) {
        stop("Iteration priors do not match prior_history.csv.")
      }
      list(project_dir = project_dir, iteration = as.integer(i),
           iteration_dir = iteration_dir, result_files = result_files,
           baseline_dir = file.path(project_dir, "results"), manifest = manifest,
           priors = priors, summary_dir = file.path(iteration_dir, "summary"))
    }, error = identity)
    if (!inherits(context, "error")) {
      if (length(rejected)) warning(paste(rejected, collapse = "; "), call. = FALSE)
      message("Using completed EM iteration ", i, ": ", context$iteration_dir)
      return(context)
    }
    rejected <- c(rejected, paste0("Iteration ", i, ": ", conditionMessage(context)))
  }
  stop("No requested completed EM iteration is available. ", paste(rejected, collapse = "; "))
}

em_join_gene <- function(file, context) {
  gene <- sub("\\.rds$", "", basename(file))
  current <- readRDS(file)
  if (!is.list(current) || !isTRUE(attr(current, "em_iteration", exact = TRUE) == context$iteration)) {
    stop("Missing or wrong em_iteration attribute for ", gene)
  }
  if (!is.null(current[["error"]])) return(current)
  if (length(current) && (is.null(names(current)) || anyDuplicated(names(current)) ||
                          anyNA(names(current)) || any(!nzchar(names(current))))) {
    stop("EM tissue names must be unique and nonempty for ", gene)
  }
  baseline <- readRDS(file.path(context$baseline_dir, basename(file)))
  if (!is.list(baseline) || !is.null(baseline[["error"]])) {
    stop("Original gene result is unavailable for comparison: ", gene)
  }
  joined <- list()
  errors <- attr(current, "tissue_errors", exact = TRUE)
  for (tissue in names(current)) {
    x <- current[[tissue]]
    original <- baseline[[tissue]]
    if (!is.list(x) || !is.list(original) || !is_successful_tissue_result(original)) {
      stop("Missing EM or original tissue result: ", gene, " / ", tissue)
    }
    fit <- x[["weighted_fit_mix"]]
    if (!is.matrix(fit$alpha) || !isTRUE(fit$converged)) {
      stop("EM fit is missing or unconverged: ", gene, " / ", tissue)
    }
    validate_predictor_map(fit, x$mix_predictor_map, "EM mix_predictor_map")
    columns <- c("predictor_index", "predictor_name", "snp", "coding")
    # Index-based overlap is meaningful only with the same biological predictors.
    if (!isTRUE(all.equal(original$mix_predictor_map[columns], x$mix_predictor_map[columns],
                          check.attributes = FALSE))) {
      stop("Original and EM predictor maps differ: ", gene, " / ", tissue)
    }
    for (field in c("n_ind", "mean_read", "median_read", "mean_phe", "median_phe")) {
      if (!is.null(original[[field]]) && !is.null(x[[field]]) &&
          !isTRUE(all.equal(original[[field]], x[[field]], tolerance = 1e-8))) {
        stop("Original and EM sample/phenotype metadata differ (", field, "): ", gene, " / ", tissue)
      }
    }
    expected_prior <- unlist(context$priors[match(tissue, context$priors$tissue),
                                           c("pi_add", "pi_rec", "pi_dom")], use.names = FALSE)
    observed_prior <- x$mix_coding_prior[c("additive", "recessive", "dominant")]
    if (length(observed_prior) != 3L || anyNA(expected_prior) ||
        !isTRUE(all.equal(unname(observed_prior), expected_prior, tolerance = 1e-10))) {
      stop("Saved EM coding prior differs from iteration priors: ", gene, " / ", tissue)
    }
    for (field in c("weighted_fit_mix", "weighted_fit_mix_lead_snp_tss_distance",
                    "weighted_mix_prior_weights", "mix_coding_prior", "em_warm_started")) {
      # Remove a baseline field when absent in EM; never retain an old weighted fit.
      original[field] <- list(x[[field]])
    }
    original$em_iteration <- context$iteration
    joined[[tissue]] <- original
  }
  absent <- setdiff(names(baseline), c(names(current), names(errors)))
  for (tissue in absent) errors[[tissue]] <- "Absent from selected EM iteration; original weighted fit was not substituted."
  attr(joined, "tissue_errors") <- errors
  joined
}

em_generate_summary <- function(project_dir, iteration = "latest") {
  context <- em_analysis_context(project_dir, iteration)
  missing_baseline <- !file.exists(file.path(context$baseline_dir, basename(context$result_files)))
  if (any(missing_baseline)) stop(sum(missing_baseline), " original gene files are missing from ", context$baseline_dir)
  source(file.path(context$project_dir, "script/analysis/generate_summary_results.R"), local = TRUE)
  # The join helper resolves summary validation functions in this environment.
  join <- em_join_gene
  environment(join) <- environment()
  output <- generate_summary_results(
    file.path(context$iteration_dir, "results"), context$summary_dir,
    result_files = context$result_files, result_reader = function(file) join(file, context))
  metadata <- list(iteration = context$iteration, iteration_dir = context$iteration_dir,
                   baseline_dir = context$baseline_dir, priors = context$priors,
                   n_files = length(context$result_files),
                   generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  saveRDS(metadata, file.path(context$summary_dir, "em_summary_metadata.rds"))
  write.csv(context$priors, file.path(context$summary_dir, "em_priors_used.csv"), row.names = FALSE)
  if (nrow(output$res_errors) || nrow(output$cs_errors)) {
    warning("Summary includes omitted/error records; inspect res_errors.csv and res_cs_errors.csv.", call. = FALSE)
  }
  invisible(output)
}

em_load_summary <- function(context) {
  metadata_file <- file.path(context$summary_dir, "em_summary_metadata.rds")
  if (!file.exists(metadata_file)) stop("Run generate_summary_results_em.R for this iteration first.")
  metadata <- readRDS(metadata_file)
  if (!identical(metadata$iteration, context$iteration) ||
      !isTRUE(all.equal(metadata$priors, context$priors))) stop("EM summary provenance mismatch; regenerate summaries.")
  env <- new.env(parent = emptyenv())
  load(file.path(context$summary_dir, "res_summary.RData"), envir = env)
  load(file.path(context$summary_dir, "res_cs_summary.RData"), envir = env)
  res <- env$res_summary
  if (!nrow(res)) stop("No successful EM comparisons in summary; inspect res_errors.csv.")
  if (!"em_iteration" %in% names(res) || anyNA(res$em_iteration) ||
      any(res$em_iteration != context$iteration) || anyDuplicated(res[c("gene", "tissue")]) ||
      !all(res$has_weighted_fit_mix %in% TRUE)) stop("Invalid or mixed-iteration EM summary.")
  cs <- env$res_cs_summary
  if (nrow(cs) && (!"em_iteration" %in% names(cs) || anyNA(cs$em_iteration) ||
                  any(cs$em_iteration != context$iteration))) stop("Invalid or mixed-iteration EM CS summary.")
  if (nrow(env$res_errors) || nrow(env$cs_errors)) {
    warning("Selected EM summary contains omissions; see its res_errors.csv and res_cs_errors.csv.", call. = FALSE)
  }
  list(res = res, cs = cs, metadata = metadata)
}
