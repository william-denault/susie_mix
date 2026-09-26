# Two additive causal SNPs, 2.5% and 5% PVE: SuSiE, slide, and slide-initialized SuSiE.
# Source/paste to load functions, or run: Rscript --vanilla sim_additive_slide_init.R CHUNK [REPS] [PVE]
init_project_dir <- local({
  root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", "")
  if (!nzchar(root)) {
    source_files <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
    source_files <- source_files[!is.na(source_files) & nzchar(source_files)]
    cli_files <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    script_file <- if (length(source_files)) tail(source_files, 1L) else
      if (length(cli_files)) cli_files[1L] else ""
    root <- if (nzchar(script_file))
      dirname(dirname(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)))) else getwd()
  }
  if (!dir.exists(root)) stop("Project directory does not exist: ", root,
    ". Set SUSIE_MIX_PROJECT_DIR to your susie_mix project.", call. = FALSE)
  normalizePath(root, winslash = "/", mustWork = TRUE)
})
.init_helpers <- file.path(init_project_dir, c(
  "script/scan_tissue_attempt/workhorse_utils.R", "script/sim/simulation_metric_helpers.R",
  "script/sim/additive_init_genotypes.R", "script/sim/sim_workhorse.R"))
.init_missing <- .init_helpers[!file.exists(.init_helpers) | file.access(.init_helpers, 4L) != 0L]
if (length(.init_missing)) stop("Cannot read simulation helper file(s):\n",
  paste(.init_missing, collapse = "\n"),
  "\nSet SUSIE_MIX_PROJECT_DIR to the project root and copy the required R files there.", call. = FALSE)
for (.init_helper in .init_helpers) source(.init_helper, local = TRUE)
rm(.init_helpers, .init_missing, .init_helper)
init_methods <- c("SuSiE", "SuSiE-slide", "SuSiE-init-slide")

# PLINK workspace size affects execution, not the genotype/phenotype design.
# Keep recording its actual value, while allowing memory-only resume changes.
init_comparable_settings <- function(settings) {
  if (identical(settings$genotype_inputs$mode, "gtex_fallback"))
    settings$genotype_inputs$plink_memory <- NULL
  settings
}

check_init_packages <- function() {
  packages <- c("susieR", "susieSlide", "data.table", "matrixStats")
  for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install required package: ", pkg)
  formal_names <- names(formals(susieR::susie))
  init_arg <- if ("model_init" %in% formal_names) "model_init" else
    if ("s_init" %in% formal_names) "s_init" else stop("This susieR version has no fitted-model initialization argument.")
  if (!"min_obs" %in% names(formals(susieSlide::susie))) stop("Install the genotype-slider version of susieSlide.")
  list(init_arg = init_arg,
       versions = setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), ""), packages))
}

simulate_additive_init_data <- function(seed, genotype_dir = NULL, n = 500L, pve = .05,
                                       genotype_source = NULL) {
  stopifnot(length(seed) == 1L, is.finite(seed), seed >= 0, seed == floor(seed),
            n >= 3, n == as.integer(n), pve > 0, pve < 1)
  if (is.null(genotype_source)) {
    if (is.null(genotype_dir)) genotype_source <- prepare_additive_genotypes() else
      genotype_source <- prepare_additive_genotypes(genotype_dir = genotype_dir)
  }
  if (identical(genotype_source$mode, "raw")) {
    # A changed region pool changes the experiment, even with the same seed.
    current_files <- list.files(genotype_source$directory, "\\.raw$", full.names = TRUE)
    if (!identical(normalizePath(current_files, winslash = "/", mustWork = TRUE), genotype_source$files))
      stop("The original simulation genotype file list changed during this run.")
    directory <- genotype_source$directory
    region <- NULL
  } else if (identical(genotype_source$mode, "gtex")) {
    # The same seed selects the same real gene at both PVEs and on reruns.
    set.seed(seed)
    index <- sample.int(nrow(genotype_source$loci), 1L)
    region <- extract_additive_region(genotype_source, index)
    on.exit(region$cleanup(), add = TRUE)
    directory <- region$directory
  } else stop("Unknown genotype source mode.")
  data <- sim_mix(pve = pve, n = n, L_add = 2L, L_rec = 0L, L_dom = 0L,
    L_prec = 0L, L_pdom = 0L, seed = seed, all_additive = FALSE,
    temp_dir = directory, return_data = TRUE)
  data$genotype_mode <- genotype_source$settings$mode
  data$gene <- region$gene
  data$region <- region$region
  data$causal_r <- cor(data$X[, data$true_pos])[1, 2]
  data
}

fit_additive_init <- function(data, fit_L = 10L, max_iter = 1000L, tol = NULL) {
  api <- check_init_packages()
  stopifnot(fit_L >= 1, fit_L == as.integer(fit_L), max_iter >= 1, max_iter == as.integer(max_iter))
  args <- list(X = data$X, y = data$y, L = fit_L, standardize = TRUE,
    estimate_prior_method = "optim", coverage = .95, min_abs_corr = .5,
    max_iter = max_iter)
  # Leave the package defaults untouched, exactly as in sim_mix().
  if (!is.null(tol)) args$tol <- tol
  fits <- warnings <- setNames(vector("list", 3L), init_methods)
  seconds <- setNames(numeric(3L), init_methods)
  for (method in init_methods) {
    fit_args <- args
    fun <- susieR::susie
    if (method == "SuSiE-slide") {
      fun <- susieSlide::susie
      fit_args$min_obs <- 5L
    } else if (method == "SuSiE-init-slide") {
      # Refit the ADDITIVE model on the original X/y, initialized by the full
      # slide fit. The new additive fit does not optimize or retain its deltas.
      fit_args[[api$init_arg]] <- fits[["SuSiE-slide"]]
    }
    messages <- character()
    elapsed <- system.time(fit <- withCallingHandlers(do.call(fun, fit_args),
      warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }))
    if (length(fit$pip) != ncol(data$X) || any(!is.finite(fit$pip))) stop("Invalid PIPs from ", method)
    fits[[method]] <- fit
    warnings[[method]] <- unique(messages)
    seconds[method] <- unname(elapsed["elapsed"])
  }
  list(fits = fits, warnings = warnings, seconds = seconds, init_argument = api$init_arg)
}

compact_additive_init <- function(data, fitted, pve, save_fits = FALSE) {
  metrics <- do.call(rbind, lapply(init_methods, function(method) {
    fit <- fitted$fits[[method]]
    counts <- cs_summary(fit$sets, fit$sets$cs, data$true_pos)
    data.frame(scenario = "Additive only", pve = pve, K = 2L,
      configuration = "add2", seed = data$seed, method = method,
      converged = isTRUE(fit$converged), as.list(counts),
      elbo = tail(fit$elbo, 1L), niter = fit$niter, sigma2 = fit$sigma2,
      seconds = unname(fitted$seconds[method]), row.names = NULL)
  }))
  fit_summaries <- lapply(fitted$fits, function(fit)
    fit[intersect(c("pip", "sets", "elbo", "niter", "converged", "V", "sigma2", "delta_cs"), names(fit))])
  out <- list(seed = data$seed, raw_file = data$raw_file, n = nrow(data$X), p = ncol(data$X),
    genotype_mode = data$genotype_mode, gene = data$gene, region = data$region,
    true_pos = data$true_pos, causal_snps = data$causal_snps,
    beta_standardized = data$beta_standardized, genetic_variance = data$genetic_variance,
    phenotype_variance = data$phenotype_variance, causal_r = data$causal_r,
    metrics = metrics, fits = fit_summaries, warnings = fitted$warnings,
    init_argument = fitted$init_argument,
    slide_delta_causal = fitted$fits[["SuSiE-slide"]]$delta[, data$true_pos, drop = FALSE])
  if (save_fits) out$full_fits <- fitted$fits
  out
}

run_additive_initialization <- function(chunk = 1L, reps_per_chunk = 100L,
    genotype_dir = Sys.getenv("SUSIE_MIX_GENOTYPE_DIR", file.path(init_project_dir, "temp_plink")),
    output_dir = file.path(init_project_dir, "simulation results/additive_slide_init_v3", sprintf("pve_%g", pve)),
    seed_base = 1000000L, n = 500L, pve = .05, fit_L = 10L,
    max_iter = 1000L, tol = NULL, save_fits = FALSE, genotype_source = NULL) {
  stopifnot(chunk >= 1, chunk == as.integer(chunk), reps_per_chunk >= 1,
    reps_per_chunk == as.integer(reps_per_chunk), seed_base >= 0,
    seed_base + chunk * reps_per_chunk <= .Machine$integer.max,
    length(pve) == 1L, is.finite(pve), pve > 0, pve < 1)
  api <- check_init_packages() # Fail once on missing packages, before the loop.
  old_threads <- data.table::setDTthreads(1L)
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  if (is.null(genotype_source)) genotype_source <- prepare_additive_genotypes(genotype_dir = genotype_dir)
  settings <- list(schema = "additive_slide_init_v3", chunk = chunk, reps_per_chunk = reps_per_chunk,
    seed_base = seed_base, K = 2L, n = n, pve = pve, fit_L = fit_L,
    max_iter = max_iter, tol = tol, save_fits = save_fits, init_argument = api$init_arg,
    package_versions = api$versions, genotype_inputs = genotype_source$settings)
  chunk_dir <- file.path(output_dir, "chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(chunk_dir, sprintf("additive_init_chunk%03d.rds", chunk))
  results <- vector("list", reps_per_chunk)
  if (file.exists(output)) {
    old <- readRDS(output)
    if (!identical(init_comparable_settings(old$settings), init_comparable_settings(settings)))
      stop("Checkpoint settings, input files or packages changed; use a new output_dir.")
    results <- old$results
    if (!is.list(results) || length(results) != reps_per_chunk) stop("Invalid checkpoint length.")
  }
  for (i in seq_len(reps_per_chunk)) {
    seed <- seed_base + (chunk - 1L) * reps_per_chunk + i
    if (!is.null(results[[i]]) && is.null(results[[i]]$error)) {
      if (results[[i]]$seed != seed) stop("Checkpoint seed mismatch.")
      next
    }
    fatal_error <- NULL
    results[[i]] <- tryCatch({
      data <- simulate_additive_init_data(seed, n = n, pve = pve, genotype_source = genotype_source)
      fitted <- fit_additive_init(data, fit_L, max_iter, tol)
      compact_additive_init(data, fitted, pve, save_fits)
    }, error = function(e) {
      if (inherits(e, "init_plink_memory_error")) fatal_error <<- e
      list(seed = seed, error = conditionMessage(e))
    })
    temporary <- paste0(output, ".tmp")
    saveRDS(list(settings = settings, results = results), temporary)
    if (!file.rename(temporary, output)) stop("Could not replace checkpoint: ", output)
    message("PVE ", 100 * pve, "%, chunk ", chunk, ": ", i, "/", reps_per_chunk, " (seed ", seed, ")",
      if (!is.null(results[[i]]$error)) paste0(" ERROR: ", results[[i]]$error) else "")
    # Save the failed seed first, then stop instead of repeating a resource error.
    if (!is.null(fatal_error)) stop(fatal_error)
  }
  n_errors <- sum(vapply(results, function(x) !is.null(x$error), logical(1)))
  message("Saved ", output, "; ", n_errors, " failed replicates. Repeating the job retries failures and skips successes.")
  if (n_errors) stop("Some replicates failed; see the saved error messages before rerunning.")
  invisible(output)
}

summarize_additive_initialization <- function(
    output_dir = file.path(init_project_dir, "simulation results/additive_slide_init_v3/pve_0.05")) {
  files <- list.files(file.path(output_dir, "chunks"), "^additive_init_chunk[0-9]+\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No initialization simulation checkpoints in ", output_dir)
  results <- list(); reference <- NULL
  for (file in files) {
    saved <- readRDS(file)
    design <- init_comparable_settings(saved$settings); design$chunk <- NULL
    if (is.null(reference)) reference <- design
    if (!identical(design, reference)) stop("Incompatible checkpoint settings: ", file)
    results <- c(results, Filter(Negate(is.null), saved$results))
  }
  seeds <- vapply(results, function(x) as.numeric(x$seed), 0)
  if (anyDuplicated(seeds)) stop("Duplicate simulation seeds across checkpoints.")
  good <- vapply(results, function(x) is.null(x$error), logical(1))
  errors <- data.frame(seed = seeds[!good],
    error = vapply(results[!good], function(x) x$error, ""))
  write.csv(errors, file.path(output_dir, "errors.csv"), row.names = FALSE)
  if (!any(good)) {
    empty <- data.frame()
    for (name in c("replicate_metrics", "metric_summary", "method_comparison", "optimization_comparison"))
      write.csv(empty, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
    message("No successful replicates in ", output_dir, "; see errors.csv.")
    return(invisible(list(metrics = empty, summary = empty, differences = empty,
                          optimization = empty, errors = errors)))
  }
  metrics <- do.call(rbind, lapply(results[good], `[[`, "metrics"))
  analysis <- summarize_metrics(metrics, methods = init_methods)
  convergence <- aggregate(list(n_converged = as.integer(metrics$converged)), metrics["method"], sum)
  analysis$summary$n_converged <- convergence$n_converged[match(analysis$summary$method, convergence$method)]
  optimization <- do.call(rbind, lapply(results[good], function(x) {
    a <- x$metrics[x$metrics$method == "SuSiE", ]
    b <- x$metrics[x$metrics$method == "SuSiE-init-slide", ]
    data.frame(seed = x$seed, pve = a$pve, causal_r = x$causal_r,
      additive_elbo = a$elbo, initialized_additive_elbo = b$elbo,
      additive_elbo_gain = b$elbo - a$elbo,
      extra_recovered = b$recovered - a$recovered,
      additive_converged = a$converged, initialized_additive_converged = b$converged)
  }))
  write.csv(metrics, file.path(output_dir, "replicate_metrics.csv"), row.names = FALSE)
  write.csv(analysis$summary, file.path(output_dir, "metric_summary.csv"), row.names = FALSE)
  write.csv(analysis$differences, file.path(output_dir, "method_comparison.csv"), row.names = FALSE)
  write.csv(optimization, file.path(output_dir, "optimization_comparison.csv"), row.names = FALSE)
  print(analysis$summary[c("method", "n_replicates", "n_converged", "coverage", "power", "purity", "cs_size")], row.names = FALSE)
  message("Summarized ", sum(good), " successful replicates; ", sum(!good), " errors. Files in ", output_dir)
  invisible(list(metrics = metrics, summary = analysis$summary, differences = analysis$differences,
                 optimization = optimization, errors = errors))
}

summarize_additive_experiment <- function(
    output_dir = file.path(init_project_dir, "simulation results/additive_slide_init_v3"),
    pves = c(.025, .05)) {
  analyses <- lapply(pves, function(pve) {
    directory <- file.path(output_dir, sprintf("pve_%g", pve))
    if (!length(list.files(file.path(directory, "chunks"), "\\.rds$"))) return(NULL)
    out <- summarize_additive_initialization(directory)
    out$errors$pve <- rep(pve, nrow(out$errors))
    out
  })
  analyses <- Filter(Negate(is.null), analyses)
  if (!length(analyses)) stop("No completed simulation checkpoints in ", output_dir)
  names_by_file <- c(metrics = "replicate_metrics.csv", summary = "metric_summary.csv",
    differences = "method_comparison.csv", optimization = "optimization_comparison.csv", errors = "errors.csv")
  combined <- setNames(lapply(names(names_by_file), function(nm)
    do.call(rbind, lapply(analyses, `[[`, nm))), names(names_by_file))
  for (nm in names(names_by_file)) write.csv(combined[[nm]], file.path(output_dir, names_by_file[[nm]]), row.names = FALSE)
  message("Combined PVE summaries saved in ", output_dir)
  invisible(combined)
}

run_additive_initialization_experiment <- function(
    pves = c(.025, .05), chunks = 1:4, reps_per_chunk = 100L,
    output_dir = file.path(init_project_dir, "simulation results/additive_slide_init_v3"),
    genotype_dir = Sys.getenv("SUSIE_MIX_GENOTYPE_DIR", file.path(init_project_dir, "temp_plink")),
    n = 500L, fit_L = 10L, max_iter = 1000L, tol = NULL, save_fits = FALSE) {
  stopifnot(length(pves) > 0L, all(is.finite(pves)), all(pves > 0 & pves < 1), !anyDuplicated(pves),
    length(chunks) > 0L, all(is.finite(chunks)), all(chunks >= 1 & chunks == floor(chunks)), !anyDuplicated(chunks),
    length(reps_per_chunk) == 1L, is.finite(reps_per_chunk),
    reps_per_chunk >= 1L, reps_per_chunk == floor(reps_per_chunk))
  check_init_packages()
  source <- prepare_additive_genotypes(genotype_dir = genotype_dir)
  failures <- character()
  for (pve in pves) for (chunk in chunks) {
    tryCatch(run_additive_initialization(chunk = chunk, reps_per_chunk = reps_per_chunk,
      output_dir = file.path(output_dir, sprintf("pve_%g", pve)),
      n = n, pve = pve, fit_L = fit_L, max_iter = max_iter, tol = tol,
      save_fits = save_fits, genotype_source = source), error = function(e) {
        if (inherits(e, "init_plink_memory_error")) stop(e)
        detail <- paste0("PVE ", pve, ", chunk ", chunk, ": ", conditionMessage(e))
        failures <<- c(failures, detail)
        message(detail)
      })
  }
  results <- summarize_additive_experiment(output_dir, pves)
  if (length(failures)) stop("Some batches reported errors; completed replicates are saved.\n",
    paste(failures, collapse = "\n"), call. = FALSE)
  invisible(results)
}

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) && args[1L] == "summarize") {
    summarize_additive_experiment()
  } else if (length(args) && args[1L] == "experiment") {
    run_additive_initialization_experiment()
  } else {
    run_additive_initialization(chunk = if (length(args)) as.integer(args[1L]) else 1L,
      reps_per_chunk = if (length(args) > 1L) as.integer(args[2L]) else 100L,
      pve = if (length(args) > 2L) as.numeric(args[3L]) else .05)
  }
}
