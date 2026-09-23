# Plot the saved simulations. No simulations or SuSiE fits are run here.
# Run this whole file with source() in R or RStudio. Only base R is needed.

# ------------------------------------------------------------
# Settings to change if needed
# ------------------------------------------------------------

project_dir <- Sys.getenv("SUSIE_MIX_PROJECT_DIR",
                         "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix")
source(file.path(project_dir, "script/sim/simulation_design.R"), local = TRUE)
source(file.path(project_dir, "script/sim/simulation_plot_helpers.R"), local = TRUE)
source(file.path(project_dir, "script/sim/simulation_metric_helpers.R"), local = TRUE)
chunk_dir <- file.path(project_dir, "simulation results/slide_v1/chunks")
output_dir <- file.path(project_dir, "simulation results/slide_v1/figures")
reuse_saved_summaries <- isTRUE(getOption("susie.sim.reuse_saved_summaries", FALSE))

pve_values <- c(0.05, 0.10, 0.20, 0.30, 0.40)
n_value <- 500
fit_L <- 10                   # SuSiE's fitted upper bound, NOT the true count.
target_coverage <- 0.95
exclude_nonconverged <- FALSE # TRUE excludes a replicate if ANY fit did not converge.
interval_level <- 0.95
proportion_ci_n <- "denominator" # n_cs for coverage/purity; n_causal for power.
                                # Set "replicates" to use the simulation count instead.

file_pattern <- "\\.RData$"
max_reps_per_file <- Inf      # For a quick preview, change this to e.g. 5.
roc_max_fpr <- 0.25           # Display range, as in Supplementary Figure 6.
fdr_max <- 0.25               # Display range for power versus empirical FDR.
write_roc_all_L <- TRUE       # Also pool all causal counts into one curve per method.
write_roc_pages <- FALSE      # Optional extra PDFs collecting the per-L figures.
write_png <- TRUE            # PDF is always written; PNG is useful for previews.

method_names <- c("SuSiE", "SuSiE-mix", "SuSiE-slide")
method_colors <- c("#2489FF", "#D81B60", "#009E73")
scenario_rows <- lapply(1:3, function(k) {
  vapply(combn(5L, k, simplify = FALSE), function(ids) {
    counts <- integer(5L)
    counts[ids] <- 1L
    sim_scenario_name(counts, display = TRUE)
  }, "")
})
pure_rows <- scenario_rows[[1L]]
mixed_rows <- c(scenario_rows[[2L]], scenario_rows[[3L]])
# Limit each figure to five scenario rows, including all ten pairs and triples.
scenario_groups <- list(pure = pure_rows,
  pairs_1 = scenario_rows[[2L]][1:5], pairs_2 = scenario_rows[[2L]][6:10],
  triples_1 = scenario_rows[[3L]][1:5], triples_2 = scenario_rows[[3L]][6:10])

# ROC counts are exact at these common PIP thresholds. Dense tails preserve
# resolution near zero/one without keeping every SNP from every run in memory.
pip_thresholds <- sort(unique(c(
  0, 10^seq(-12, -2, length.out = 301), seq(0.01, 0.99, by = 0.001),
  1 - 10^seq(-2, -12, length.out = 301), 1
)))

# ------------------------------------------------------------
# Small helpers: scenarios, PIP counts, and saved CS summaries
# ------------------------------------------------------------

scenario_name <- function(a, r, d, pr = 0, pd = 0) {
  sim_scenario_name(c(a, r, d, pr, pd), display = TRUE)
}

pip_counts <- function(pip, truth, thresholds = pip_thresholds) {
  if (!length(pip) || anyNA(pip) || any(!is.finite(pip)) ||
      any(pip < -1e-10 | pip > 1 + 1e-10)) stop("Invalid saved SNP PIPs.")
  if (!length(truth) || anyNA(truth) ||
      any(truth < 1 | truth > length(pip) | truth != as.integer(truth))) {
    stop("Causal SNP indices do not match the saved PIPs.")
  }
  truth <- unique(truth)
  bin <- findInterval(pmin(1, pmax(0, pip)), thresholds)
  all_bins <- tabulate(bin, nbins = length(thresholds))
  causal_bins <- tabulate(bin[truth], nbins = length(thresholds))
  # All equal PIPs enter together: a discovery means PIP >= threshold.
  tp <- rev(cumsum(rev(causal_bins)))
  fp <- rev(cumsum(rev(all_bins - causal_bins)))
  cbind(tp = c(tp, 0), fp = c(fp, 0)) # Final row is threshold = Inf.
}

# This is pooled empirical FDP (often labelled empirical FDR in power-FDR plots),
# not FPR and not the unweighted mean of each replicate's FDP.
add_fdr <- function(counts) {
  discoveries <- counts$tp + counts$fp
  counts$fdr <- ifelse(discoveries > 0, counts$fp / discoveries, 0)
  counts
}

# Stack all L values by summing TP/FP at each common threshold, then recalculate
# rates. Averaging the per-L rates would give different weights to the SNPs.
pool_curve_counts <- function(counts) {
  by <- c("scenario", "pve", "method", "threshold")
  pooled <- aggregate(counts[c("tp", "fp")], counts[by], sum)
  pooled$tpr <- pooled$fpr <- NA_real_
  groups <- split(seq_len(nrow(pooled)),
                  interaction(pooled[c("scenario", "pve", "method")], drop = TRUE))
  for (ids in groups) {
    baseline <- ids[pooled$threshold[ids] == 0]
    if (length(baseline) != 1 || pooled$tp[baseline] <= 0 || pooled$fp[baseline] <= 0) {
      stop("Pooled curves require threshold zero with all causal and noncausal SNPs.")
    }
    pooled$tpr[ids] <- pooled$tp[ids] / pooled$tp[baseline]
    pooled$fpr[ids] <- pooled$fp[ids] / pooled$fp[baseline]
  }
  add_fdr(pooled)
}

# ------------------------------------------------------------
# Read one checkpoint at a time
# ------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (reuse_saved_summaries) {
  replicates <- readRDS(file.path(output_dir, "replicate_metrics.rds"))
  if (!"cs_size_sum_sq" %in% names(replicates)) {
    cat("Older summaries lack CS-size sums of squares; rebuilding them once from checkpoints for Gaussian intervals.\n")
    reuse_saved_summaries <- FALSE
  }
}
if (!reuse_saved_summaries) {
files <- list.files(chunk_dir, pattern = file_pattern, full.names = TRUE)
if (!length(files)) stop("No matching .RData files found in: ", chunk_dir)

# These fields also identify the condition for error-only saved replicates.
file_info <- sim_parse_checkpoints(files)
file_info <- file_info[file_info$n == n_value & file_info$fit_L == fit_L &
                       file_info$pve %in% pve_values, , drop = FALSE]
if (!nrow(file_info)) stop("No chunks match n_value, fit_L and pve_values.")

# Prefer the larger checkpoint when old/new files repeat the same condition
# and seed. The seed key below removes overlaps, including partial overlaps.
file_info <- file_info[order(-file_info$requested_reps,
                            -as.numeric(file.info(file_info$file)$mtime)), ]
seen <- new.env(hash = TRUE, parent = emptyenv())
curve_counts <- list()
replicate_rows <- list()
audit_rows <- list()
calibration_state <- new.env(hash = TRUE, parent = emptyenv())

for (f in seq_len(nrow(file_info))) {
  info <- file_info[f, ]
  cat("Reading", f, "of", nrow(file_info), ":", basename(info$file), "\n")
  saved <- new.env(parent = emptyenv())
  load(info$file, envir = saved)
  if (!exists("results", envir = saved, inherits = FALSE) || !is.list(saved$results)) {
    stop("No results list in ", info$file)
  }
  scenario <- scenario_name(info$L_add, info$L_rec, info$L_dom, info$L_prec, info$L_pdom)
  K <- sum(as.numeric(info[sim_count_columns]))
  if (!K %in% 1:5) stop("Expected 1 to 5 true causal SNPs in ", info$file)
  configuration <- sim_configuration(as.numeric(info[sim_count_columns]))
  audit <- data.frame(file = basename(info$file), scenario = scenario,
                      configuration = configuration, pve = info$pve, K = K,
                      saved = length(saved$results), examined = 0L, included = 0L,
                      errors = 0L, duplicates = 0L, nonconverged = 0L,
                      excluded_nonconverged = 0L, error_messages = "")
  file_rows <- list()
  row_number <- 0L
  file_seeds <- rep(NA_real_, length(saved$results))
  file_calibration <- lapply(method_names, function(m) matrix(0, length(saved$results), 31))

  for (o in seq_len(min(length(saved$results), max_reps_per_file))) {
    x <- saved$results[[o]]
    audit$examined <- audit$examined + 1L
    if (is.null(x$seed) || length(x$seed) != 1 || !is.finite(x$seed)) {
      stop("Missing seed in ", basename(info$file), ", replicate ", o)
    }
    # Including the causal allocation prevents combining distinct simulations
    # that intentionally share a seed across PVE and architecture conditions.
    key <- paste(configuration, info$pve, info$n, info$fit_L, x$seed, sep = "|")
    if (exists(key, envir = seen, inherits = FALSE)) {
      audit$duplicates <- audit$duplicates + 1L
      next
    }
    if (!is.null(x$error)) {
      audit$errors <- audit$errors + 1L
      audit$error_messages <- paste(unique(c(audit$error_messages, x$error)), collapse = " | ")
      # A successful copy in another checkpoint may still be included.
      next
    }
    s <- x$settings
    if (is.null(s$L_prec)) s$L_prec <- 0L
    if (is.null(s$L_pdom)) s$L_pdom <- 0L
    if (is.null(s) || !isTRUE(all.equal(as.numeric(c(s$L_add, s$L_rec, s$L_dom,
                                                   s$L_prec, s$L_pdom, s$n, s$L, s$pve)),
                                      as.numeric(info[c(sim_count_columns,
                                                        "n", "fit_L", "pve")])))) {
      stop("Saved settings disagree with filename in ", basename(info$file))
    }
    if (isTRUE(s$all_additive)) stop("Paired all-additive controls need a separate plot.")
    if (!isTRUE(all.equal(c(s$min_maf, s$hwe_thresh, s$min_n_rec), c(.05, 1e-8, 5)))) {
      stop("Unexpected genotype-QC settings; do not silently pool these runs.")
    }
    if (length(x$susie_pip) != length(x$susie_mix_pip_snp) ||
        !length(x$susie_mix_pip_snp)) {
      stop("Missing mixed SNP-level PIPs. Do not substitute coding-level PIPs.")
    }
    if (length(unique(x$true_pos)) != K) stop("Causal count disagrees with settings.")
    if (!identical(as.integer(table(factor(x$causal_coding, levels = names(sim_effect_delta)))),
                   as.integer(unlist(s[sim_count_columns]))))
      stop("Causal coding labels disagree with settings.")
    if ("SuSiE-slide" %in% method_names) {
      if (!identical(s$schema_version, sim_schema_version) ||
          !isTRUE(all.equal(c(s$delta_prec, s$delta_pdom, s$slide_min_obs), c(-.5, .5, 5))))
        stop("Missing or unexpected slider design. Use new three-method checkpoints.")
      if (length(x$susie_slide_pip) != length(x$susie_pip))
        stop("Missing or incompatible slider SNP PIPs.")
      if (!identical(x$causal_delta, unname(sim_effect_delta[x$causal_coding])))
        stop("Saved causal delta values disagree with the design.")
    }
    mt <- match(method_names, x$metrics$method)
    if (anyNA(mt)) stop("Missing method in the saved metrics.")
    converged <- x$metrics$converged[mt]
    if (length(converged) != length(method_names) || anyNA(converged)) stop("Missing convergence status.")
    assign(key, TRUE, envir = seen)
    if (!all(converged)) {
      audit$nonconverged <- audit$nonconverged + 1L
      if (exclude_nonconverged) {
        audit$excluded_nonconverged <- audit$excluded_nonconverged + 1L
        next
      }
    }

    cs <- list(SuSiE = x$susie_cs$cs, `SuSiE-mix` = x$cs_mix_as_additive_indices,
               `SuSiE-slide` = x$susie_slide_cs$cs)[method_names]
    sets <- list(SuSiE = x$susie_cs, `SuSiE-mix` = x$susie_mix_cs,
                 `SuSiE-slide` = x$susie_slide_cs)[method_names]
    pips <- list(SuSiE = x$susie_pip, `SuSiE-mix` = x$susie_mix_pip_snp,
                 `SuSiE-slide` = x$susie_slide_pip)[method_names]
    file_seeds[o] <- x$seed
    for (m in seq_along(method_names)) {
      file_calibration[[m]][o, ] <- calibration_bin_totals(pips[[m]], x$true_pos)
      summary <- cs_summary(sets[[m]], cs[[m]], x$true_pos)
      if (summary["n_cs"] > 0 &&
          !isTRUE(all.equal(sets[[m]]$requested_coverage, target_coverage))) {
        stop("The requested CS coverage differs from target_coverage.")
      }
      row_number <- row_number + 1L
      file_rows[[row_number]] <- data.frame(
        scenario = scenario, pve = info$pve, K = K, configuration = configuration,
        seed = x$seed, method = method_names[m], converged = converged[m],
        as.list(summary), check.names = FALSE
      )
      group_key <- paste(scenario, info$pve, K, method_names[m], sep = "|")
      counts <- pip_counts(pips[[m]], x$true_pos)
      if (is.null(curve_counts[[group_key]])) {
        curve_counts[[group_key]] <- counts
      } else {
        curve_counts[[group_key]] <- curve_counts[[group_key]] + counts
      }
    }
    audit$included <- audit$included + 1L
  }
  audit_rows[[f]] <- audit
  included_ids <- which(!is.na(file_seeds))
  for (m in seq_along(method_names)) {
    calibration_merge(calibration_state,
      data.frame(scenario = scenario, pve = info$pve, K = K, method = method_names[m]),
      file_seeds[included_ids], file_calibration[[m]][included_ids, , drop = FALSE])
  }
  # Combine each file's small summaries now; do not retain thousands of
  # individual data.frame objects or any of its large saved PIP vectors.
  replicate_rows[[f]] <- if (length(file_rows)) do.call(rbind, file_rows) else NULL
  rm(saved, file_rows)
}

audit <- do.call(rbind, audit_rows)
write.csv(audit, file.path(output_dir, "file_audit.csv"), row.names = FALSE)
if (!sum(audit$included)) stop("No usable simulations; see file_audit.csv.")
replicates <- do.call(rbind, replicate_rows)
rm(replicate_rows, seen)

# ------------------------------------------------------------
# Save compact counts for analytic intervals and later redraws
# ------------------------------------------------------------

group_vars <- c("scenario", "pve", "K", "method")

# Every causal allocation at a given total K is pooled. If jobs are incomplete,
# allocations with more completed replicates contribute more observations.
# Save the detailed counts so that imbalance is visible.
configuration_counts <- aggregate(list(n_replicates = replicates$seed),
                                  replicates[c(group_vars, "configuration")], length)
write.csv(configuration_counts, file.path(output_dir, "configuration_counts.csv"), row.names = FALSE)
saveRDS(replicates, file.path(output_dir, "replicate_metrics.rds"))

curve_rows <- lapply(names(curve_counts), function(key) {
  fields <- strsplit(key, "|", fixed = TRUE)[[1]]
  counts <- curve_counts[[key]]
  total_causal <- counts[1, "tp"]
  total_null <- counts[1, "fp"]
  if (total_null == 0) stop("ROC is undefined without noncausal SNPs.")
  data.frame(scenario = fields[1], pve = as.numeric(fields[2]),
             K = as.integer(fields[3]), method = fields[4],
             threshold = c(pip_thresholds, Inf),
             tp = counts[, "tp"], fp = counts[, "fp"],
             tpr = counts[, "tp"] / total_causal,
             fpr = counts[, "fp"] / total_null)
})
roc_table <- add_fdr(do.call(rbind, curve_rows))
saveRDS(roc_table, file.path(output_dir, "roc_counts.rds"))
pooled_roc_table <- pool_curve_counts(roc_table)
saveRDS(pooled_roc_table, file.path(output_dir, "roc_counts_all_L.rds"))
rm(curve_counts, curve_rows)
calibration_groups <- as.list(calibration_state)
saveRDS(list(selection = calibration_selection(replicates, method_names),
             signature = calibration_source_signature(file_info$file), groups = calibration_groups),
        file.path(output_dir, "pip_calibration_seed_counts.rds"))
rm(calibration_state)
} else {
  cat("Refreshing figures from saved summaries; no models are fitted.\n")
  audit <- read.csv(file.path(output_dir, "file_audit.csv"), stringsAsFactors = FALSE)
  roc_table <- readRDS(file.path(output_dir, "roc_counts.rds"))
  pooled_roc_table <- readRDS(file.path(output_dir, "roc_counts_all_L.rds"))
  calibration_groups <- collect_calibration(file.path(chunk_dir, audit$file), replicates,
    method_names, file.path(output_dir, "pip_calibration_seed_counts.rds"))
}
cat("Calculating", round(100 * interval_level), "% analytic normal intervals (no bootstrap)...\n")
metric_analysis <- summarize_metrics(replicates, level = interval_level,
  methods = method_names, proportion_n = proportion_ci_n)
summary_table <- metric_analysis$summary
write.csv(summary_table, file.path(output_dir, "metric_summary.csv"), row.names = FALSE)
write.csv(metric_analysis$differences, file.path(output_dir, "method_comparison.csv"), row.names = FALSE)
calibration_table <- calibration_summary(calibration_groups)
pooled_calibration_table <- calibration_summary(calibration_groups, pool_K = TRUE)
write.csv(calibration_table, file.path(output_dir, "pip_calibration.csv"), row.names = FALSE)
write.csv(pooled_calibration_table, file.path(output_dir, "pip_calibration_all_L.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Draw the panels in the style of the supplementary figures
# ------------------------------------------------------------

figure_data <- function(metric, scenarios, only_K = NULL, methods = method_names) {
  d <- if (metric == "pip_calibration" && is.null(only_K)) pooled_calibration_table else
    if (metric == "pip_calibration") calibration_table else
    if (metric %in% c("roc", "power_fdr") && is.null(only_K)) pooled_roc_table else
    if (metric %in% c("roc", "power_fdr")) roc_table else summary_table
  d <- d[d$scenario %in% scenarios & d$pve %in% pve_values &
           d$method %in% methods, , drop = FALSE]
  if (!is.null(only_K)) d <- d[d$K %in% only_K, , drop = FALSE]
  d
}

draw_figure <- function(metric, scenarios, only_K = NULL, y_limits = NULL,
                        methods = method_names) {
  method_ids <- match(methods, method_names)
  stopifnot(length(method_ids) > 0, !anyNA(method_ids))
  nr <- length(scenarios)
  nc <- length(pve_values)
  panels <- matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE)
  layout(cbind(panels, nr * nc + seq_len(nr)), widths = c(rep(1, nc), 1.15))
  par(oma = c(6, 3.5, 3, 0.3), mar = c(2.3, 3.15, 1.7, 0.65),
      mgp = c(1.3, 0.4, 0), tcl = -0.2, family = "sans", cex = 0.9, xpd = FALSE)
  is_roc <- metric == "roc"
  is_fdr <- metric == "power_fdr"
  is_calibration <- metric == "pip_calibration"
  is_curve <- is_roc || is_fdr
  is_pip <- is_curve || is_calibration
  d_figure <- figure_data(metric, scenarios, only_K, methods)
  # Every panel in this output shares the range of all displayed estimates.
  ylim <- if (is.null(y_limits)) figure_y_limits(metric, d_figure, c(0, fdr_max)) else y_limits
  yticks <- if (metric %in% c("power", "roc", "pip_calibration")) seq(0, 1, .25) else panel_y_ticks(ylim)

  for (i in seq_along(scenarios)) {
    min_K <- length(strsplit(scenarios[i], " + ", fixed = TRUE)[[1L]])
    for (j in seq_along(pve_values)) {
      d <- d_figure[d_figure$scenario == scenarios[i] & d_figure$pve == pve_values[j], , drop = FALSE]
      if (is_calibration) {
        xlim <- c(0, 1)
        xticks <- seq(0, 1, .25)
      } else if (is_curve) {
        xlim <- c(0, if (is_fdr) fdr_max else roc_max_fpr)
        xticks <- pretty(xlim, n = 5)
        xticks <- xticks[xticks >= 0 & xticks <= xlim[2]]
      } else {
        xlim <- c(min_K - 0.35, 5.35)
        xticks <- min_K:5
      }
      plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
           axes = FALSE, xlab = "", ylab = "")
      abline(v = xticks, h = yticks, col = "#DEDEDE", lwd = 0.8)
      if (metric == "coverage") abline(h = target_coverage, lty = 2, lwd = 1.2)
      if (is_roc || is_calibration) abline(0, 1, col = "#888888", lty = 2)
      axis(1, at = xticks, cex.axis = 0.9, col = "#777777")
      axis(2, at = yticks, labels = format(signif(yticks, 3), trim = TRUE),
           las = 1, cex.axis = .78, col = "#777777")
      par(xpd = NA)
      if (i == 1) mtext(paste0("PVE = ", round(100 * pve_values[j]), "%"),
                        side = 3, line = 0.55, font = 2)
      par(xpd = FALSE)
      if (!nrow(d)) {
        label <- if (!is.null(only_K) && only_K < min_K) "Not applicable" else "No saved results"
        text(mean(xlim), mean(ylim), label, col = "#777777", cex = .8)
        next
      }
      for (m in method_ids) {
        dm <- d[d$method == method_names[m], , drop = FALSE]
        if (is_calibration) {
          dm <- dm[order(dm$bin), ]
          ok <- is.finite(dm$lower) & is.finite(dm$upper) & is.finite(dm$mean_pip)
          if (any(ok)) segments(dm$mean_pip[ok], dm$lower[ok], dm$mean_pip[ok], dm$upper[ok],
                               col = adjustcolor(method_colors[m], alpha.f = .6))
          shown <- is.finite(dm$mean_pip) & is.finite(dm$frequency)
          if (any(shown)) points(dm$mean_pip[shown], dm$frequency[shown],
                                 pch = c(16, 17, 15)[m], cex = .75,
                                 col = method_colors[m], xpd = NA)
        } else if (is_curve) {
          # The all-L table already pools the underlying TP/FP counts.
          # FDR need not increase monotonically: preserve threshold order,
          # rather than sorting FDR or reporting an optimized envelope.
          dm <- dm[order(dm$threshold, decreasing = TRUE), ]
          seg <- visible_curve_segments(if (is_fdr) dm$fdr else dm$fpr, dm$tpr, xlim)
          if (nrow(seg)) segments(seg[, "x0"], pmax(ylim[1], pmin(ylim[2], seg[, "y0"])),
                                 seg[, "x1"], pmax(ylim[1], pmin(ylim[2], seg[, "y1"])),
                                 col = method_colors[m], lty = 1, lwd = 1.5)
        } else {
          xpos <- dm$K + seq(-.12, .12, length.out = length(method_names))[m]
          lo <- dm[[paste0(metric, "_lo")]]
          hi <- dm[[paste0(metric, "_hi")]]
          if (!is.null(lo) && !is.null(hi)) {
            ok <- is.finite(lo) & is.finite(hi) & hi >= ylim[1] & lo <= ylim[2]
            bar_color <- adjustcolor(method_colors[m], alpha.f = .65)
            # Clip explicitly: some Windows raster devices mishandle a segment
            # outside the plot followed by symbols drawn on the boundary.
            if (any(ok)) segments(xpos[ok], pmax(lo[ok], ylim[1]), xpos[ok], pmin(hi[ok], ylim[2]), col = bar_color)
            lo_inside <- ok & lo >= ylim[1] & lo <= ylim[2]
            hi_inside <- ok & hi >= ylim[1] & hi <= ylim[2]
            if (any(lo_inside)) segments(xpos[lo_inside] - .035, lo[lo_inside], xpos[lo_inside] + .035, lo[lo_inside], col = bar_color)
            if (any(hi_inside)) segments(xpos[hi_inside] - .035, hi[hi_inside], xpos[hi_inside] + .035, hi[hi_inside], col = bar_color)
          }
          shown <- is.finite(dm[[metric]]) & dm[[metric]] >= ylim[1] & dm[[metric]] <= ylim[2]
          if (any(shown)) points(xpos[shown], dm[[metric]][shown],
                                 pch = 16, cex = 0.95, col = method_colors[m], xpd = NA)
        }
      }
    }
  }
  # A narrow column at the right labels each scenario, without crowding panels.
  par(mar = c(0, 0, 0, 0), xpd = NA)
  for (s in scenarios) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    label <- gsub(" + ", "\n+ ", s, fixed = TRUE)
    text(.05, .5, label, adj = c(0, .5), font = 2, cex = .85, xpd = NA)
  }
  titles <- c(coverage = "Credible-set coverage", purity = "Credible-set purity",
              power = "Causal SNP recovery by credible sets", cs_size = "Credible-set size",
              roc = "Detection of causal SNPs using PIPs", power_fdr = "Power versus empirical FDR using PIPs",
              pip_calibration = "PIP calibration")
  title_text <- paste0(titles[metric], "  |  n = ", n_value)
  if (!is.null(only_K)) title_text <- paste0(title_text, "  |  L = ", only_K, " causal SNP",
                                          if (only_K == 1) "" else "s")
  if (is_pip && is.null(only_K)) title_text <- paste0(title_text, "  |  All L pooled")
  if (length(methods) == 1L) title_text <- paste0(title_text, "  |  ", methods)
  mtext(title_text, side = 3, outer = TRUE, line = 1.2, font = 2, cex = 1.1)
  mtext(if (is_calibration) "Mean PIP within bin" else if (is_fdr) "Empirical FDR" else if (is_roc) "False positive rate" else "Number of causal SNPs",
        side = 1, outer = TRUE, line = 0.5)
  ylab <- c(coverage = "Coverage", purity = "Mean minimum absolute correlation",
            power = "Power", cs_size = "Mean number of unique SNPs per CS",
            roc = "True positive rate (power)", power_fdr = "Power (true positive rate)",
            pip_calibration = "Fraction of SNPs in bin that are causal")
  mtext(ylab[metric], side = 2, outer = TRUE, line = 1.8)

  # Draw a common legend in the outer bottom margin.
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend(.5, .045, legend = methods, col = method_colors[method_ids],
         pch = if (is_curve) NA else if (is_calibration) c(16, 17, 15)[method_ids] else 16, lty = if (is_curve) 1 else NA,
         lwd = 1.5, horiz = TRUE, xjust = .5, yjust = .5, bty = "n", cex = .95)
  if (is_calibration) {
    text(.5, .015, "Ten equal-width PIP bins; error bars: +/-2 empirical SE across seed blocks", cex = .85)
  } else if (!is_curve && paste0(metric, "_lo") %in% names(summary_table)) {
    text(.5, .015, paste0(round(100 * interval_level),
      if (metric == "cs_size") "% Gaussian CI: mean +/- z * SD/sqrt(n CS)" else
        paste0("% normal CI: p +/- z * sqrt(p(1-p)/n); n = ",
          if (proportion_ci_n == "replicates") "simulation replicates" else
            if (metric == "power") "true causal SNPs" else "reported CS"),
      if (metric != "power") "; intervals may be clipped by shared y-axis limits" else ""), cex = .85)
  }
}

save_figure <- function(metric, scenarios, name, only_K = NULL, methods = method_names) {
  height <- 2.0 * length(scenarios) + 1.4
  pdf(file.path(output_dir, paste0(name, ".pdf")), width = 16, height = height,
      useDingbats = FALSE)
  tryCatch(draw_figure(metric, scenarios, only_K, methods = methods), finally = dev.off())
  if (write_png) {
    png(file.path(output_dir, paste0(name, ".png")), width = 16, height = height,
        units = "in", res = 180,
        type = if (capabilities("cairo")) "cairo-png" else getOption("bitmapType"))
    tryCatch(draw_figure(metric, scenarios, only_K, methods = methods), finally = dev.off())
  }
}

for (metric in c("coverage", "purity", "power", "cs_size")) {
  for (group in names(scenario_groups))
    save_figure(metric, scenario_groups[[group]], paste0(metric, "_", group))
}

save_roc_figures <- function(metric = "roc") {
  for (group in names(scenario_groups)) {
    scenarios <- scenario_groups[[group]]
    min_K <- length(strsplit(scenarios[1L], " + ", fixed = TRUE)[[1L]])
    ks <- min_K:5
    for (k in ks) {
      save_figure(metric, scenarios, paste0(metric, "_", group, "_L", k), only_K = k)
    }
    if (write_roc_all_L) {
      save_figure(metric, scenarios, paste0(metric, "_", group, "_all_L"))
    }
    if (write_roc_pages) {
      # Also share the range across pages of the optional combined PDF.
      y_limits <- figure_y_limits(metric, figure_data(metric, scenarios, only_K = ks), c(0, fdr_max))
      pdf(file.path(output_dir, paste0(metric, "_", group, "_by_L.pdf")),
          width = 16, height = 2.0 * length(scenarios) + 1.4, useDingbats = FALSE)
      tryCatch(for (k in ks) draw_figure(metric, scenarios, only_K = k, y_limits = y_limits), finally = dev.off())
    }
  }
}
save_roc_figures()
save_roc_figures("power_fdr")

save_calibration_figures <- function() {
  # Calibration exports always pool all L and show one method per figure.
  for (group in names(scenario_groups)) {
    for (method in method_names) {
      method_suffix <- tolower(gsub("[^[:alnum:]]+", "_", method))
      save_figure("pip_calibration", scenario_groups[[group]],
        paste0("pip_calibration_", group, "_all_L_", method_suffix), methods = method)
    }
  }
  # Remove only obsolete generated calibration figures, after all replacements
  # have been written successfully. Preserve CSV summaries and the RDS cache.
  obsolete <- list.files(output_dir,
    pattern = "^pip_calibration_(pure|pairs_[12]|triples_[12])_(L[1-5]|by_L|all_L)\\.(pdf|png)$",
    full.names = TRUE)
  if (length(obsolete)) {
    removed <- file.remove(obsolete)
    cat("Removed", sum(removed), "obsolete calibration plot files.\n")
    if (!all(removed)) warning("Some obsolete calibration plots could not be removed; check whether they are open.")
  }
}
save_calibration_figures()

cat("\nFigures and summary tables saved in:", output_dir, "\n")
cat("Included", sum(audit$included), "unique simulations; skipped",
    sum(audit$duplicates), "duplicate copies; recorded", sum(audit$errors), "error entries.\n")
cat("Replicates with a nonconverged fit:", sum(audit$nonconverged),
    "; excluded:", sum(audit$excluded_nonconverged), "\n")
cat("L in each ROC title is the TRUE causal count; fitted SuSiE L =", fit_L, ".\n")
cat("Check configuration_counts.csv for incomplete or unbalanced conditions.\n")
