#!/usr/bin/env Rscript

# Simulate GTEx-like eQTL summary statistics from real-LD genotype matrices,
# then compare in-sample, ADSP-like, and UKB-like LD references with and without
# the population-bias correction.

default_config <- list(
  input_dir = "/home/gw/Documents/susie_ash_test/chat_test",
  output_dir = "/home/gw/Documents/susie_ash_test/chat_test/lambda_pop_sim_results",
  susie_repo = "/home/gw/GIT/susieR",
  simxqtl_repo = "/home/gw/GIT/simxQTL",
  n_reps = 20L,
  seed = 20260501L,
  p_max = 2500L,
  n_ref = 500L,
  L = 10L,
  max_iter = 80L,
  coverage = 0.95,
  min_abs_corr = 0.5,
  ld_proxy_threshold = 0.8,
  h2g = 0.15,
  n_sparse = 3L,
  n_oligogenic = 0L,
  n_inf = 0L,
  prop_h2_sparse = 1.00,
  prop_h2_oligogenic = 0.00,
  prop_h2_infinitesimal = 0.00,
  adsp_delta = 0.02,
  ukb_delta = 0.35,
  lambda_zero_tol = 0.01,
  include_map_qc = TRUE,
  verbose_fit = FALSE,
  verbose_reps = 1L,
  smoke = FALSE
)

vlog <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  message(msg)
  flush.console()
}

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  cfg <- defaults
  if (!length(args)) {
    return(cfg)
  }
  for (arg in args) {
    if (!grepl("^--", arg)) {
      stop("Arguments must be --name=value; got: ", arg)
    }
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[1])
    value <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    if (!key %in% names(cfg)) {
      stop("Unknown argument --", kv[1])
    }
    old <- cfg[[key]]
    if (is.logical(old)) {
      cfg[[key]] <- tolower(value) %in% c("true", "t", "1", "yes", "y")
    } else if (is.integer(old)) {
      cfg[[key]] <- as.integer(value)
    } else if (is.numeric(old)) {
      cfg[[key]] <- as.numeric(value)
    } else {
      cfg[[key]] <- value
    }
  }
  if (isTRUE(cfg$smoke)) {
    cfg$n_reps <- min(cfg$n_reps, 1L)
    cfg$p_max <- min(cfg$p_max, 300L)
    cfg$n_ref <- min(cfg$n_ref, 200L)
    cfg$L <- min(cfg$L, 6L)
    cfg$max_iter <- min(cfg$max_iter, 20L)
  }
  cfg
}

load_local_packages <- function(cfg) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required to load local susieR and simxQTL repos.")
  }
  pkgload::load_all(cfg$susie_repo, quiet = TRUE)
  pkgload::load_all(cfg$simxqtl_repo, quiet = TRUE)
}

write_json <- function(x, file) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, file, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    capture.output(str(x), file = file)
  }
}

standardize_matrix <- function(X) {
  X <- as.matrix(X)
  X <- scale(X, center = TRUE, scale = TRUE)
  X[is.na(X)] <- 0
  storage.mode(X) <- "double"
  X
}

select_variant_window <- function(G, p_max, seed) {
  G <- as.matrix(G)
  p <- ncol(G)
  if (p <= p_max) {
    return(standardize_matrix(G))
  }
  set.seed(seed)
  start <- sample.int(p - p_max + 1L, 1L)
  standardize_matrix(G[, start:(start + p_max - 1L), drop = FALSE])
}

read_genotype_files <- function(input_dir) {
  files <- list.files(input_dir, pattern = "[.]rds$", recursive = TRUE,
                      full.names = TRUE)
  ok <- logical(length(files))
  dims <- vector("list", length(files))
  for (i in seq_along(files)) {
    obj <- tryCatch(readRDS(files[i]), error = function(e) NULL)
    X <- if (!is.null(obj$G)) obj$G else obj$X
    ok[i] <- is.matrix(X) || is.data.frame(X)
    if (ok[i]) {
      dims[[i]] <- dim(X)
    }
  }
  data.frame(
    file = files[ok],
    n = vapply(dims[ok], `[`, numeric(1), 1L),
    p = vapply(dims[ok], `[`, numeric(1), 2L),
    stringsAsFactors = FALSE
  )
}

make_reference_panel <- function(G, n_ref, delta, seed) {
  set.seed(seed)
  n <- nrow(G)
  idx <- sample.int(n, n_ref, replace = n_ref > n)
  X <- G[idx, , drop = FALSE]
  if (delta > 0) {
    E <- matrix(rnorm(n_ref * ncol(G)), n_ref, ncol(G))
    X <- sqrt(1 - delta) * X + sqrt(delta) * E
  }
  standardize_matrix(X)
}

make_z_scores <- function(X, y) {
  z <- calc_z(X, y, center = TRUE, scale = FALSE)
  z[!is.finite(z)] <- 0
  as.numeric(z)
}

fit_susie_rss <- function(z, X_ref, n_target, cfg, method, verbose_fit = FALSE) {
  R_mismatch <- switch(method,
                   no_finite_no_bias = "none",
                   finite_only = "none",
                   bias_map = "map",
                   bias_map_qc = "map_qc",
                   stop("Unknown method: ", method))
  R_finite <- switch(method,
                     no_finite_no_bias = NULL,
                     finite_only = TRUE,
                     bias_map = TRUE,
                     bias_map_qc = TRUE,
                     stop("Unknown method: ", method))
  args <- list(
    z = z,
    X = X_ref,
    n = n_target,
    L = cfg$L,
    coverage = cfg$coverage,
	    min_abs_corr = cfg$min_abs_corr,
	    max_iter = cfg$max_iter,
	    R_finite = R_finite,
	    R_mismatch = R_mismatch,
	    estimate_residual_variance = FALSE,
	    verbose = isTRUE(verbose_fit)
	  )
  if (isTRUE(verbose_fit)) {
    return(do.call(susie_rss, args))
  }
  withCallingHandlers(do.call(susie_rss, args),
                      message = function(m) invokeRestart("muffleMessage"),
                      warning = function(w) invokeRestart("muffleWarning"))
}

extract_lambda_table <- function(fit, rep_id, panel, method) {
  diag <- fit$R_finite_diagnostics
  lb <- diag$lambda_bias
  if (is.null(lb)) {
    lb <- rep(NA_real_, nrow(fit$alpha))
  }
  bc <- diag$B_corrected
  if (is.null(bc)) {
    bc <- rep(NA_real_, length(lb))
  }
  data.frame(
    rep = rep_id,
    panel = panel,
    method = method,
    effect = seq_along(lb),
    R_finite_B = if (is.null(diag$B)) NA_real_ else as.numeric(diag$B),
    lambda_pop = as.numeric(lb),
    B_corrected = as.numeric(bc),
    stringsAsFactors = FALSE
  )
}

max_abs_ld_to_causal <- function(X_target, idx, causal) {
  if (!length(idx) || !length(causal)) {
    return(0)
  }
  idx <- intersect(idx, seq_len(ncol(X_target)))
  causal <- intersect(causal, seq_len(ncol(X_target)))
  if (!length(idx) || !length(causal)) {
    return(0)
  }
  C <- crossprod(X_target[, idx, drop = FALSE],
                 X_target[, causal, drop = FALSE]) / (nrow(X_target) - 1)
  max(abs(C))
}

causal_detected_by_cs <- function(X_target, cs, causal, ld_threshold) {
  if (!length(cs) || !length(causal)) {
    return(integer(0))
  }
  detected <- integer(0)
  for (j in causal) {
    if (j %in% cs) {
      detected <- c(detected, j)
    } else {
      C <- crossprod(X_target[, cs, drop = FALSE],
                     X_target[, j, drop = FALSE]) / (nrow(X_target) - 1)
      if (max(abs(C)) >= ld_threshold) {
        detected <- c(detected, j)
      }
    }
  }
  unique(detected)
}

extract_cs_metrics <- function(fit, X_target, causal, rep_id, panel, method,
                               ld_threshold) {
  cs_list <- fit$sets$cs
  if (is.null(cs_list) || !length(cs_list)) {
    return(data.frame(
      rep = rep_id, panel = panel, method = method, cs_name = NA_character_,
      cs_size = 0L, exact_hit = FALSE, proxy_hit = FALSE,
      max_abs_ld_to_causal = NA_real_, detected_causal = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  rows <- vector("list", length(cs_list))
  for (i in seq_along(cs_list)) {
    cs <- as.integer(cs_list[[i]])
    detected <- causal_detected_by_cs(X_target, cs, causal, ld_threshold)
    rows[[i]] <- data.frame(
      rep = rep_id,
      panel = panel,
      method = method,
      cs_name = names(cs_list)[i],
      cs_size = length(cs),
      exact_hit = any(cs %in% causal),
      proxy_hit = length(detected) > 0,
      max_abs_ld_to_causal = max_abs_ld_to_causal(X_target, cs, causal),
      detected_causal = paste(detected, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_fit <- function(fit, X_target, causal, rep_id, panel, method,
                          elapsed, err = NULL, ld_threshold = 0.8) {
  if (!is.null(err)) {
    return(data.frame(
      rep = rep_id, panel = panel, method = method, status = "error",
      error = conditionMessage(err), elapsed_sec = elapsed,
      n_cs = NA_integer_, cs_tp_exact = NA_integer_, cs_tp_proxy = NA_integer_,
      cs_fp_proxy = NA_integer_, cs_fdr_proxy = NA_real_,
      causal_recall_proxy = NA_real_, top1_is_causal = NA,
      top1_ld_proxy = NA, max_pip_causal = NA_real_,
      mean_lambda_pop = NA_real_, max_lambda_pop = NA_real_,
      nonzero_lambda_pop = NA_integer_, R_finite_B = NA_real_,
      mean_B_corrected = NA_real_,
      max_per_variable_penalty = NA_real_, Q_art = NA_real_,
      artifact_flag = NA, mode_label = NA_character_, converged = NA,
      stringsAsFactors = FALSE
    ))
  }

  cs_rows <- extract_cs_metrics(fit, X_target, causal, rep_id, panel, method,
                                ld_threshold)
  has_cs <- !all(is.na(cs_rows$cs_name))
  detected <- integer(0)
  if (has_cs) {
    det_str <- cs_rows$detected_causal[nzchar(cs_rows$detected_causal)]
    detected <- unique(as.integer(unlist(strsplit(paste(det_str, collapse = ";"),
                                             ";", fixed = TRUE))))
    detected <- detected[!is.na(detected)]
  }
  top1 <- which.max(fit$pip)
  diag <- fit$R_finite_diagnostics
  lb <- diag$lambda_bias
  penalty <- diag$per_variable_penalty
  data.frame(
    rep = rep_id,
    panel = panel,
    method = method,
    status = "ok",
    error = NA_character_,
    elapsed_sec = elapsed,
    n_cs = if (has_cs) nrow(cs_rows) else 0L,
    cs_tp_exact = if (has_cs) sum(cs_rows$exact_hit) else 0L,
    cs_tp_proxy = if (has_cs) sum(cs_rows$proxy_hit) else 0L,
    cs_fp_proxy = if (has_cs) sum(!cs_rows$proxy_hit) else 0L,
    cs_fdr_proxy = if (has_cs) mean(!cs_rows$proxy_hit) else NA_real_,
    causal_recall_proxy = length(intersect(detected, causal)) / length(causal),
    top1_is_causal = top1 %in% causal,
    top1_ld_proxy = max_abs_ld_to_causal(X_target, top1, causal) >= ld_threshold,
    max_pip_causal = max(fit$pip[causal]),
    R_finite_B = if (is.null(diag$B)) NA_real_ else as.numeric(diag$B),
    mean_lambda_pop = if (is.null(lb)) NA_real_ else mean(lb),
    max_lambda_pop = if (is.null(lb)) NA_real_ else max(lb),
    nonzero_lambda_pop = if (is.null(lb)) NA_integer_ else sum(lb > 0),
    mean_B_corrected = if (is.null(diag$B_corrected)) NA_real_ else mean(diag$B_corrected),
    max_per_variable_penalty = if (is.null(penalty)) NA_real_ else max(penalty),
    Q_art = if (is.null(diag$Q_art)) NA_real_ else diag$Q_art,
    artifact_flag = if (is.null(diag$artifact_flag)) NA else isTRUE(diag$artifact_flag),
    mode_label = if (is.null(diag$mode_label)) NA_character_ else as.character(diag$mode_label),
    converged = isTRUE(fit$converged),
    stringsAsFactors = FALSE
  )
}

write_ai_readme <- function(cfg, out_dir) {
  lines <- c(
    "# Lambda-pop LD-bias simulation outputs",
    "",
    "Primary files for AI parsing:",
    "",
    "- `per_fit_metrics.csv`: one row per replicate, panel, and method. Key methods are `no_finite_no_bias`, `finite_only`, `bias_map`, and `bias_map_qc`; key columns are `R_finite_B`, `max_lambda_pop`, `mean_lambda_pop`, `mean_B_corrected`, `causal_recall_proxy`, `cs_fdr_proxy`, `top1_is_causal`, and `max_pip_causal`.",
    "- `per_effect_lambda.csv`: per-effect `R_finite_B`, `lambda_pop`, and `B_corrected` estimates.",
    "- `cs_metrics.csv`: one row per credible set, with exact and LD-proxy truth labels.",
    "- `replicate_metadata.csv`: source file, dimensions, causal indices, and realized h2.",
    "- `aggregate_summary.csv`: mean/median summaries grouped by panel and method.",
    "- `run_config.json`: exact simulation settings.",
    "",
    "Expected checks:",
    "",
    "1. In-sample LD with `method == bias_map` should have `lambda_pop` equal or very close to zero.",
    "2. `ADSP_like` should have smaller `lambda_pop` and larger `B_corrected` than `UKB_like`.",
    "3. `bias_map` should reduce false credible sets under biased LD without losing too much causal recall.",
    "",
    paste0("Configured ADSP delta = ", cfg$adsp_delta,
           "; UKB delta = ", cfg$ukb_delta,
           "; LD proxy threshold = ", cfg$ld_proxy_threshold, ".")
  )
  writeLines(lines, file.path(out_dir, "AI_PARSE_ME.md"))
}

aggregate_metrics <- function(metrics) {
  ok <- metrics[metrics$status == "ok", , drop = FALSE]
  if (!nrow(ok)) {
    return(data.frame())
  }
  groups <- unique(ok[, c("panel", "method")])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- ok$panel == groups$panel[i] & ok$method == groups$method[i]
    x <- ok[idx, , drop = FALSE]
    rows[[i]] <- data.frame(
      panel = groups$panel[i],
      method = groups$method[i],
      n_ok = nrow(x),
      mean_max_lambda_pop = mean(x$max_lambda_pop, na.rm = TRUE),
      median_max_lambda_pop = median(x$max_lambda_pop, na.rm = TRUE),
      mean_lambda_pop = mean(x$mean_lambda_pop, na.rm = TRUE),
      R_finite_B = mean(x$R_finite_B, na.rm = TRUE),
      mean_B_corrected = mean(x$mean_B_corrected, na.rm = TRUE),
      mean_causal_recall_proxy = mean(x$causal_recall_proxy, na.rm = TRUE),
      mean_cs_fdr_proxy = mean(x$cs_fdr_proxy, na.rm = TRUE),
      mean_n_cs = mean(x$n_cs, na.rm = TRUE),
      mean_max_pip_causal = mean(x$max_pip_causal, na.rm = TRUE),
      artifact_flag_rate = mean(x$artifact_flag, na.rm = TRUE),
      mean_Q_art = mean(x$Q_art, na.rm = TRUE),
      top1_causal_rate = mean(x$top1_is_causal, na.rm = TRUE),
      top1_proxy_rate = mean(x$top1_ld_proxy, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

metric_value <- function(metrics, rep_i, panel, method, column) {
  idx <- metrics$rep == rep_i & metrics$panel == panel & metrics$method == method
  if (!any(idx)) {
    return(NA_real_)
  }
  value <- metrics[idx, column][1]
  if (is.logical(value)) {
    return(as.numeric(value))
  }
  as.numeric(value)
}

mean_metric <- function(metrics, panel, method, column) {
  idx <- metrics$status == "ok" & metrics$panel == panel & metrics$method == method
  if (!any(idx)) {
    return(NA_real_)
  }
  mean(as.numeric(metrics[idx, column]), na.rm = TRUE)
}

fmt_metric <- function(x, digits = 3) {
  if (!is.finite(x)) {
    return("NA")
  }
  formatC(x, format = "fg", digits = digits)
}

print_progress_line <- function(metrics, rep_i, cfg) {
  in_lambda <- metric_value(metrics, rep_i, "in_sample", "bias_map",
                            "max_lambda_pop")
  adsp_lambda <- metric_value(metrics, rep_i, "ADSP_like", "bias_map",
                              "max_lambda_pop")
  ukb_lambda <- metric_value(metrics, rep_i, "UKB_like", "bias_map",
                             "max_lambda_pop")
  adsp_recall <- metric_value(metrics, rep_i, "ADSP_like", "bias_map",
                              "causal_recall_proxy")
  ukb_recall <- metric_value(metrics, rep_i, "UKB_like", "bias_map",
                             "causal_recall_proxy")
  adsp_fdr <- metric_value(metrics, rep_i, "ADSP_like", "bias_map",
                           "cs_fdr_proxy")
  ukb_fdr <- metric_value(metrics, rep_i, "UKB_like", "bias_map",
                          "cs_fdr_proxy")

  mean_in_lambda <- mean_metric(metrics, "in_sample", "bias_map",
                                "max_lambda_pop")
  mean_adsp_lambda <- mean_metric(metrics, "ADSP_like", "bias_map",
                                  "max_lambda_pop")
  mean_ukb_lambda <- mean_metric(metrics, "UKB_like", "bias_map",
                                 "max_lambda_pop")
  mean_adsp_recall <- mean_metric(metrics, "ADSP_like", "bias_map",
                                  "causal_recall_proxy")
  mean_ukb_recall <- mean_metric(metrics, "UKB_like", "bias_map",
                                 "causal_recall_proxy")

  in_ok <- is.finite(in_lambda) && in_lambda <= cfg$lambda_zero_tol
  adsp_lt_ukb <- is.finite(adsp_lambda) && is.finite(ukb_lambda) &&
    adsp_lambda < ukb_lambda
  mean_adsp_lt_ukb <- is.finite(mean_adsp_lambda) &&
    is.finite(mean_ukb_lambda) && mean_adsp_lambda < mean_ukb_lambda

  message(
    "PROGRESS rep=", rep_i,
    " current: in_lambda=", fmt_metric(in_lambda),
    " in_zero=", in_ok,
    " ADSP_lambda=", fmt_metric(adsp_lambda),
    " UKB_lambda=", fmt_metric(ukb_lambda),
    " ADSP<UKB=", adsp_lt_ukb,
    " ADSP_recall=", fmt_metric(adsp_recall),
    " UKB_recall=", fmt_metric(ukb_recall),
    " ADSP_FDR=", fmt_metric(adsp_fdr),
    " UKB_FDR=", fmt_metric(ukb_fdr),
    " | running_mean: in_lambda=", fmt_metric(mean_in_lambda),
    " ADSP_lambda=", fmt_metric(mean_adsp_lambda),
    " UKB_lambda=", fmt_metric(mean_ukb_lambda),
    " ADSP<UKB=", mean_adsp_lt_ukb,
    " ADSP_recall=", fmt_metric(mean_adsp_recall),
    " UKB_recall=", fmt_metric(mean_ukb_recall)
  )
}

run_one_fit <- function(z, X_ref, n_target, cfg, method, rep_id, panel,
                        X_target, causal, verbose_fit = FALSE) {
  vlog("  fit start  rep=", rep_id, " panel=", panel, " method=", method,
       " (n_ref=", nrow(X_ref), ", p=", ncol(X_ref), ")")
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(fit_susie_rss(z, X_ref, n_target, cfg, method, verbose_fit),
                  error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error")) {
    vlog("  fit ERROR  rep=", rep_id, " panel=", panel, " method=", method,
         " in ", fmt_metric(elapsed), "s: ", conditionMessage(fit))
    return(list(
      metric = summarize_fit(NULL, X_target, causal, rep_id, panel, method,
                             elapsed, fit, cfg$ld_proxy_threshold),
      lambda = data.frame(),
      cs = data.frame(),
      fit = NULL
    ))
  }
  metric <- summarize_fit(fit, X_target, causal, rep_id, panel, method,
                          elapsed, NULL, cfg$ld_proxy_threshold)
  vlog("  fit done   rep=", rep_id, " panel=", panel, " method=", method,
       " in ", fmt_metric(elapsed), "s | lambda=", fmt_metric(metric$max_lambda_pop),
       " Bcorr=", fmt_metric(metric$mean_B_corrected),
       " n_cs=", metric$n_cs, " recall=", fmt_metric(metric$causal_recall_proxy),
       " Q_art=", fmt_metric(metric$Q_art))
  list(
    metric = metric,
    lambda = extract_lambda_table(fit, rep_id, panel, method),
    cs = extract_cs_metrics(fit, X_target, causal, rep_id, panel, method,
                            cfg$ld_proxy_threshold),
    fit = fit
  )
}

main <- function() {
  cfg <- parse_args(default_config)
  set.seed(cfg$seed)
  vlog("=== lambda_pop LD-bias simulation ===")
  vlog("smoke=", cfg$smoke, "  n_reps=", cfg$n_reps, "  L=", cfg$L,
       "  max_iter=", cfg$max_iter, "  p_max=", cfg$p_max,
       "  n_ref=", cfg$n_ref)
  vlog("h2g=", cfg$h2g, "  n_sparse=", cfg$n_sparse,
       "  ADSP_delta=", cfg$adsp_delta, "  UKB_delta=", cfg$ukb_delta,
       "  include_map_qc=", cfg$include_map_qc,
       "  verbose_fit=", cfg$verbose_fit, "  verbose_reps=", cfg$verbose_reps)
  vlog("output_dir=", cfg$output_dir)
  vlog("Loading susieR and simxQTL via pkgload::load_all...")
  load_local_packages(cfg)
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_json(cfg, file.path(cfg$output_dir, "run_config.json"))
  write_ai_readme(cfg, cfg$output_dir)

  vlog("Indexing genotype RDS files in ", cfg$input_dir, " ...")
  geno_index <- read_genotype_files(cfg$input_dir)
  if (!nrow(geno_index)) {
    stop("No RDS files with G or X matrices found in ", cfg$input_dir)
  }
  geno_index <- geno_index[order(-geno_index$n, -geno_index$p, geno_index$file), ]
  write.csv(geno_index, file.path(cfg$output_dir, "genotype_file_index.csv"),
            row.names = FALSE)
  vlog("Found ", nrow(geno_index), " usable genotype files; using top ",
       min(cfg$n_reps, nrow(geno_index)))

  chosen <- geno_index$file[seq_len(min(cfg$n_reps, nrow(geno_index)))]
  all_metrics <- list()
  all_lambda <- list()
  all_cs <- list()
  all_meta <- list()
  compact_fits <- list()

  for (rep_i in seq_along(chosen)) {
    rep_seed <- as.integer((cfg$seed + rep_i * 1000L) %% 1000000L)
    if (rep_seed <= 0L) {
      rep_seed <- rep_i
    }
    vlog("Replicate ", rep_i, "/", length(chosen), ": ", chosen[rep_i])
    obj <- readRDS(chosen[rep_i])
    G0 <- if (!is.null(obj$G)) obj$G else obj$X
    G <- select_variant_window(G0, cfg$p_max, rep_seed)
    vlog("  rep_seed=", rep_seed, "  genotype dim after window=",
         nrow(G), "x", ncol(G))

    sim <- generate_cis_qtl_data(
      G = G,
      h2g = cfg$h2g,
      prop_h2_sparse = cfg$prop_h2_sparse,
      prop_h2_oligogenic = cfg$prop_h2_oligogenic,
      prop_h2_infinitesimal = cfg$prop_h2_infinitesimal,
      n_sparse = cfg$n_sparse,
      n_oligogenic = cfg$n_oligogenic,
      n_inf = cfg$n_inf,
      standardize = TRUE,
      independent = TRUE,
      seed = rep_seed
    )
    X_target <- standardize_matrix(sim$G)
    y <- as.numeric(scale(sim$y, center = TRUE, scale = TRUE))
    z <- make_z_scores(X_target, y)
    causal <- sort(unique(sim$sparse_indices))
    causal <- causal[causal >= 1L & causal <= ncol(X_target)]
    if (!length(causal)) {
      causal <- sort(unique(which(sim$beta != 0)))
    }
    vlog("  simulated y: h2g_real=", fmt_metric(sim$h2g),
         "  causal=[", paste(causal, collapse = ","), "]")

    panels <- list(
      in_sample = X_target,
      ADSP_like = make_reference_panel(X_target, cfg$n_ref, cfg$adsp_delta,
                                       rep_seed + 11L),
      UKB_like = make_reference_panel(X_target, cfg$n_ref, cfg$ukb_delta,
                                      rep_seed + 23L)
    )
    methods <- c("no_finite_no_bias", "finite_only", "bias_map")
    if (isTRUE(cfg$include_map_qc)) {
      methods <- c(methods, "bias_map_qc")
    }

    all_meta[[rep_i]] <- data.frame(
      rep = rep_i,
      source_file = chosen[rep_i],
      n_target = nrow(X_target),
      p = ncol(X_target),
      n_causal_sparse = length(causal),
      causal_sparse = paste(causal, collapse = ";"),
      h2g_realized = sim$h2g,
      h2_sparse_realized = sim$h2_sparse,
      h2_oligogenic_realized = sim$h2_oligogenic,
      h2_infinitesimal_realized = sim$h2_infinitesimal,
      stringsAsFactors = FALSE
    )

    for (panel_name in names(panels)) {
      for (method in methods) {
        verbose_fit <- isTRUE(cfg$verbose_fit) &&
          rep_i <= cfg$verbose_reps &&
          method %in% c("bias_map", "bias_map_qc")
        if (verbose_fit) {
          message("VERBOSE_FIT_BEGIN rep=", rep_i, " panel=", panel_name,
                  " method=", method)
        }
        res <- run_one_fit(z, panels[[panel_name]], nrow(X_target), cfg,
                           method, rep_i, panel_name, X_target, causal,
                           verbose_fit = verbose_fit)
        if (verbose_fit) {
          message("VERBOSE_FIT_END rep=", rep_i, " panel=", panel_name,
                  " method=", method)
        }
        key <- paste(rep_i, panel_name, method, sep = "__")
        all_metrics[[key]] <- res$metric
        all_lambda[[key]] <- res$lambda
        all_cs[[key]] <- res$cs
        compact_fits[[key]] <- list(
          pip = if (!is.null(res$fit)) res$fit$pip else NULL,
          sets = if (!is.null(res$fit)) res$fit$sets else NULL,
          R_finite_diagnostics =
            if (!is.null(res$fit)) res$fit$R_finite_diagnostics else NULL
        )
      }
    }

    metrics_so_far <- do.call(rbind, all_metrics)
    write.csv(metrics_so_far, file.path(cfg$output_dir, "per_fit_metrics.csv"),
              row.names = FALSE)
    write.csv(do.call(rbind, all_meta),
              file.path(cfg$output_dir, "replicate_metadata.csv"),
              row.names = FALSE)
    if (length(all_lambda)) {
      write.csv(do.call(rbind, all_lambda),
                file.path(cfg$output_dir, "per_effect_lambda.csv"),
                row.names = FALSE)
    }
    if (length(all_cs)) {
      write.csv(do.call(rbind, all_cs), file.path(cfg$output_dir, "cs_metrics.csv"),
                row.names = FALSE)
    }
    write.csv(aggregate_metrics(metrics_so_far),
              file.path(cfg$output_dir, "aggregate_summary.csv"),
              row.names = FALSE)
    print_progress_line(metrics_so_far, rep_i, cfg)
  }

  final_metrics <- do.call(rbind, all_metrics)
  final <- list(
    config = cfg,
    metrics = final_metrics,
    lambda = if (length(all_lambda)) do.call(rbind, all_lambda) else data.frame(),
    cs = if (length(all_cs)) do.call(rbind, all_cs) else data.frame(),
    metadata = do.call(rbind, all_meta),
    aggregate = aggregate_metrics(final_metrics),
    compact_fits = compact_fits
  )
  saveRDS(final, file.path(cfg$output_dir, "all_results.rds"))
  write_json(list(
    completed_at = as.character(Sys.time()),
    output_dir = normalizePath(cfg$output_dir, mustWork = FALSE),
    n_reps_requested = cfg$n_reps,
    n_reps_completed = length(chosen),
    files = c("AI_PARSE_ME.md", "run_config.json", "genotype_file_index.csv",
              "replicate_metadata.csv", "per_fit_metrics.csv",
              "per_effect_lambda.csv", "cs_metrics.csv",
              "aggregate_summary.csv", "all_results.rds")
  ), file.path(cfg$output_dir, "manifest.json"))
  vlog("Done. Results written to: ", cfg$output_dir)
}

main()
