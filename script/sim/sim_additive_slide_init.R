# Two additive causal SNPs, 5% PVE: SuSiE, SuSiE-slide, and slide-initialized SuSiE.
# Source to load functions, or run: Rscript --vanilla sim_additive_slide_init.R CHUNK [REPS]
.init_sources <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
.init_file <- if (any(nzchar(.init_sources))) tail(.init_sources[nzchar(.init_sources)], 1L) else
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.init_file <- normalizePath(.init_file, winslash = "/", mustWork = TRUE)
init_project_dir <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", dirname(dirname(dirname(.init_file))))
source(file.path(dirname(.init_file), "../scan_tissue_attempt/workhorse_utils.R"), local = TRUE)
source(file.path(dirname(.init_file), "simulation_metric_helpers.R"), local = TRUE)
rm(.init_sources, .init_file)
init_methods <- c("SuSiE", "SuSiE-slide", "SuSiE-init-slide")

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

simulate_additive_init_data <- function(seed, genotype_dir, n = 500L, pve = .05) {
  stopifnot(length(seed) == 1L, is.finite(seed), seed >= 0, seed == floor(seed),
            n >= 3, n == as.integer(n), pve > 0, pve < 1)
  set.seed(seed)
  files <- list.files(genotype_dir, pattern = "\\.raw$", full.names = TRUE)
  if (!length(files)) stop("No PLINK .raw files found in: ", genotype_dir)
  raw_file <- files[sample.int(length(files), 1)]
  raw <- data.table::fread(raw_file, data.table = FALSE)
  if (ncol(raw) < 8L || !"IID" %in% names(raw)) stop("Expected PLINK .raw input with six metadata columns.")
  X <- as.matrix(raw[, -(1:6), drop = FALSE])
  storage.mode(X) <- "double"
  rownames(X) <- raw$IID
  X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
  if (!all(X %in% 0:2)) stop("Expected hard-call additive genotypes 0/1/2.")
  X <- X[, matrixStats::colSds(X) > 0, drop = FALSE]
  # Same minor-allele orientation, MAF/HWE and donor filtering as sim_mix().
  X <- qc_filter_geno(X, hwe_thresh = 1e-8, maf_min = .05)$X
  if (nrow(X) < n) stop("Not enough donors in ", raw_file)
  X <- X[sample.int(nrow(X), n), , drop = FALSE]
  X <- X[, colSums(X) >= 5 & matrixStats::colSds(X) > 0, drop = FALSE]
  if (ncol(X) < 2L) stop("Fewer than two eligible SNPs in ", raw_file)
  storage.mode(X) <- "double"
  true_pos <- sample.int(ncol(X), 2L)
  Z <- scale(X[, true_pos, drop = FALSE])
  beta <- sample(c(-1, 1), 2L, replace = TRUE)
  g <- drop(Z %*% beta)
  if (!is.finite(var(g)) || var(g) <= 0) stop("Zero genetic variance for this causal pair.")
  beta <- beta * sqrt(pve / var(g))
  g <- drop(Z %*% beta)
  y <- g + rnorm(n, sd = sqrt(1 - pve))
  list(X = X, y = y, true_pos = true_pos, causal_snps = colnames(X)[true_pos],
       raw_file = raw_file, beta_standardized = beta, genetic_variance = var(g),
       phenotype_variance = var(y), causal_r = cor(X[, true_pos])[1, 2], seed = seed)
}

fit_additive_init <- function(data, fit_L = 10L, max_iter = 1000L, tol = 1e-3) {
  api <- check_init_packages()
  stopifnot(fit_L >= 1, fit_L == as.integer(fit_L), max_iter >= 1, max_iter == as.integer(max_iter))
  args <- list(X = data$X, y = data$y, L = fit_L, standardize = TRUE,
    estimate_prior_method = "optim", coverage = .95, min_abs_corr = .5,
    max_iter = max_iter, tol = tol, verbose = FALSE)
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
    output_dir = file.path(init_project_dir, "simulation results/additive_slide_init_v1"),
    seed_base = 1000000L, n = 500L, pve = .05, fit_L = 10L,
    max_iter = 1000L, tol = 1e-3, save_fits = FALSE) {
  stopifnot(chunk >= 1, chunk == as.integer(chunk), reps_per_chunk >= 1,
    reps_per_chunk == as.integer(reps_per_chunk), seed_base >= 0,
    seed_base + chunk * reps_per_chunk <= .Machine$integer.max)
  api <- check_init_packages() # Fail once on missing packages, before the loop.
  genotype_dir <- normalizePath(genotype_dir, winslash = "/", mustWork = TRUE)
  files <- list.files(genotype_dir, pattern = "\\.raw$", full.names = TRUE)
  if (!length(files)) stop("No PLINK .raw files found in: ", genotype_dir)
  info <- file.info(files)
  settings <- list(schema = "additive_slide_init_v1", chunk = chunk, reps_per_chunk = reps_per_chunk,
    seed_base = seed_base, K = 2L, n = n, pve = pve, fit_L = fit_L,
    max_iter = max_iter, tol = tol, save_fits = save_fits, init_argument = api$init_arg,
    package_versions = api$versions, genotype_dir = genotype_dir,
    genotype_files = data.frame(file = basename(files), size = info$size, mtime = as.numeric(info$mtime)))
  chunk_dir <- file.path(output_dir, "chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(chunk_dir, sprintf("additive_init_chunk%03d.rds", chunk))
  results <- vector("list", reps_per_chunk)
  if (file.exists(output)) {
    old <- readRDS(output)
    if (!identical(old$settings, settings)) stop("Checkpoint settings, input files or packages changed; use a new output_dir.")
    results <- old$results
    if (!is.list(results) || length(results) != reps_per_chunk) stop("Invalid checkpoint length.")
  }
  for (i in seq_len(reps_per_chunk)) {
    seed <- seed_base + (chunk - 1L) * reps_per_chunk + i
    if (!is.null(results[[i]]) && is.null(results[[i]]$error)) {
      if (results[[i]]$seed != seed) stop("Checkpoint seed mismatch.")
      next
    }
    results[[i]] <- tryCatch({
      data <- simulate_additive_init_data(seed, genotype_dir, n, pve)
      fitted <- fit_additive_init(data, fit_L, max_iter, tol)
      compact_additive_init(data, fitted, pve, save_fits)
    }, error = function(e) list(seed = seed, error = conditionMessage(e)))
    temporary <- paste0(output, ".tmp")
    saveRDS(list(settings = settings, results = results), temporary)
    if (!file.rename(temporary, output)) stop("Could not replace checkpoint: ", output)
    message("Chunk ", chunk, ": ", i, "/", reps_per_chunk, " (seed ", seed, ")",
      if (!is.null(results[[i]]$error)) paste0(" ERROR: ", results[[i]]$error) else "")
  }
  n_errors <- sum(vapply(results, function(x) !is.null(x$error), logical(1)))
  message("Saved ", output, "; ", n_errors, " failed replicates. Repeating the job retries failures and skips successes.")
  if (n_errors) stop("Some replicates failed; see the saved error messages before rerunning.")
  invisible(output)
}

summarize_additive_initialization <- function(
    output_dir = file.path(init_project_dir, "simulation results/additive_slide_init_v1")) {
  files <- list.files(file.path(output_dir, "chunks"), "^additive_init_chunk[0-9]+\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No initialization simulation checkpoints in ", output_dir)
  results <- list(); reference <- NULL
  for (file in files) {
    saved <- readRDS(file)
    design <- saved$settings; design$chunk <- NULL
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
  if (!any(good)) stop("No successful replicates; see errors.csv.")
  metrics <- do.call(rbind, lapply(results[good], `[[`, "metrics"))
  analysis <- summarize_metrics(metrics, methods = init_methods)
  convergence <- aggregate(list(n_converged = as.integer(metrics$converged)), metrics["method"], sum)
  analysis$summary$n_converged <- convergence$n_converged[match(analysis$summary$method, convergence$method)]
  optimization <- do.call(rbind, lapply(results[good], function(x) {
    a <- x$metrics[x$metrics$method == "SuSiE", ]
    b <- x$metrics[x$metrics$method == "SuSiE-init-slide", ]
    data.frame(seed = x$seed, causal_r = x$causal_r,
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

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) && args[1L] == "summarize") {
    summarize_additive_initialization()
  } else {
    run_additive_initialization(chunk = if (length(args)) as.integer(args[1L]) else 1L,
      reps_per_chunk = if (length(args) > 1L) as.integer(args[2L]) else 100L)
  }
}
