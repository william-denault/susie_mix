# Plot the saved simulations. No simulations or SuSiE fits are run here.
# Run this whole file with source() in R or RStudio. Only base R is needed.

# ------------------------------------------------------------
# Settings to change if needed
# ------------------------------------------------------------

project_dir <- "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix"
chunk_dir <- file.path(project_dir, "simulation results/chunks")
output_dir <- file.path(project_dir, "simulation results/figures")

pve_values <- c(0.10, 0.20, 0.30, 0.40)
n_value <- 500
fit_L <- 10                   # SuSiE's fitted upper bound, NOT the true count.
target_coverage <- 0.95
exclude_nonconverged <- FALSE # TRUE excludes a replicate if EITHER fit failed.

file_pattern <- "\\.RData$"
max_reps_per_file <- Inf      # For a quick preview, change this to e.g. 5.
roc_max_fpr <- 0.25           # Display range, as in Supplementary Figure 6.
write_roc_pages <- TRUE       # Also write one ROC page per causal SNP count.
write_png <- TRUE            # PDF is always written; PNG is useful for previews.

method_names <- c("SuSiE", "SuSiE-mix")
method_colors <- c("#2489FF", "#D81B60")
pure_rows <- c("Additive only", "Dominant only", "Recessive only")
mixed_rows <- c("Additive + dominant", "Additive + recessive",
                "Recessive + dominant", "Additive + recessive + dominant")

# ROC counts are exact at these common PIP thresholds. Dense tails preserve
# resolution near zero/one without keeping every SNP from every run in memory.
pip_thresholds <- sort(unique(c(
  0, 10^seq(-12, -2, length.out = 301), seq(0.01, 0.99, by = 0.001),
  1 - 10^seq(-2, -12, length.out = 301), 1
)))

# ------------------------------------------------------------
# Small helpers: scenarios, PIP counts, and saved CS summaries
# ------------------------------------------------------------

scenario_name <- function(a, r, d) {
  if (a > 0 && r == 0 && d == 0) return(pure_rows[1])
  if (a == 0 && r == 0 && d > 0) return(pure_rows[2])
  if (a == 0 && r > 0 && d == 0) return(pure_rows[3])
  if (a > 0 && r == 0 && d > 0) return(mixed_rows[1])
  if (a > 0 && r > 0 && d == 0) return(mixed_rows[2])
  if (a == 0 && r > 0 && d > 0) return(mixed_rows[3])
  if (a > 0 && r > 0 && d > 0) return(mixed_rows[4])
  stop("A simulation must have at least one causal SNP.")
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

cs_summary <- function(sets, cs_snps, truth) {
  n_cs <- length(cs_snps)
  if (length(sets$cs) != n_cs) stop("Saved CS mapping has the wrong length.")
  hit <- vapply(cs_snps, function(cs) any(cs %in% truth), logical(1))
  purity <- sets$purity$min.abs.corr
  if (n_cs && (length(purity) != n_cs || any(!is.finite(purity)))) {
    stop("Missing or invalid minimum-correlation purity for a reported CS.")
  }
  c(n_cs = n_cs, covered_cs = sum(hit), purity_sum = sum(purity),
    recovered = sum(unique(truth) %in% unlist(cs_snps, use.names = FALSE)),
    n_causal = length(unique(truth)))
}

# ------------------------------------------------------------
# Read one checkpoint at a time
# ------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(chunk_dir, pattern = file_pattern, full.names = TRUE)
if (!length(files)) stop("No matching .RData files found in: ", chunk_dir)

# These fields also identify the condition for error-only saved replicates.
filename_pattern <- paste0(
  "^.*_add([0-9]+)_rec([0-9]+)_dom([0-9]+)_n([0-9]+)_L([0-9]+)",
  "_pve([^_]+)_seed([^_]+)_reps([0-9]+)_chunk([0-9]+)\\.RData$"
)
parts <- regmatches(basename(files), regexec(filename_pattern, basename(files)))
if (any(lengths(parts) != 10)) stop("Unrecognized chunk filename: ",
                                      basename(files[which(lengths(parts) != 10)[1]]))
file_info <- as.data.frame(do.call(rbind, lapply(parts, function(z) as.numeric(z[-1]))))
names(file_info) <- c("L_add", "L_rec", "L_dom", "n", "fit_L", "pve",
                      "seed_base", "requested_reps", "chunk")
file_info$file <- files
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

for (f in seq_len(nrow(file_info))) {
  info <- file_info[f, ]
  cat("Reading", f, "of", nrow(file_info), ":", basename(info$file), "\n")
  saved <- new.env(parent = emptyenv())
  load(info$file, envir = saved)
  if (!exists("results", envir = saved, inherits = FALSE) || !is.list(saved$results)) {
    stop("No results list in ", info$file)
  }
  scenario <- scenario_name(info$L_add, info$L_rec, info$L_dom)
  K <- info$L_add + info$L_rec + info$L_dom
  if (!K %in% 1:5) stop("Expected 1 to 5 true causal SNPs in ", info$file)
  configuration <- paste0("add", info$L_add, "_rec", info$L_rec, "_dom", info$L_dom)
  audit <- data.frame(file = basename(info$file), scenario = scenario,
                      configuration = configuration, pve = info$pve, K = K,
                      saved = length(saved$results), examined = 0L, included = 0L,
                      errors = 0L, duplicates = 0L, nonconverged = 0L,
                      excluded_nonconverged = 0L, error_messages = "")
  file_rows <- list()
  row_number <- 0L

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
    if (is.null(s) || !isTRUE(all.equal(as.numeric(c(s$L_add, s$L_rec, s$L_dom,
                                                   s$n, s$L, s$pve)),
                                      as.numeric(info[c("L_add", "L_rec", "L_dom",
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
    mt <- match(method_names, x$metrics$method)
    if (anyNA(mt)) stop("Missing method in the saved metrics.")
    converged <- x$metrics$converged[mt]
    if (length(converged) != 2 || anyNA(converged)) stop("Missing convergence status.")
    assign(key, TRUE, envir = seen)
    if (!all(converged)) {
      audit$nonconverged <- audit$nonconverged + 1L
      if (exclude_nonconverged) {
        audit$excluded_nonconverged <- audit$excluded_nonconverged + 1L
        next
      }
    }

    cs <- list(x$susie_cs$cs, x$cs_mix_as_additive_indices)
    sets <- list(x$susie_cs, x$susie_mix_cs)
    pips <- list(x$susie_pip, x$susie_mix_pip_snp)
    for (m in seq_along(method_names)) {
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
  # Combine each file's small summaries now; do not retain thousands of
  # individual data.frame objects or any of its large saved PIP vectors.
  replicate_rows[[f]] <- if (length(file_rows)) do.call(rbind, file_rows) else NULL
  rm(saved, file_rows)
}

audit <- do.call(rbind, audit_rows)
write.csv(audit, file.path(output_dir, "file_audit.csv"), row.names = FALSE)
if (!length(replicate_rows)) stop("No usable simulations; see file_audit.csv.")
replicates <- do.call(rbind, replicate_rows)
rm(replicate_rows, seen)

# ------------------------------------------------------------
# Pool counts, rather than averaging per-run CS coverage
# ------------------------------------------------------------

group_vars <- c("scenario", "pve", "K", "method")
count_vars <- c("n_cs", "covered_cs", "purity_sum", "recovered", "n_causal")
summary_table <- aggregate(replicates[count_vars], replicates[group_vars], sum)
rep_counts <- aggregate(list(n_replicates = replicates$seed), replicates[group_vars], length)
config_counts <- aggregate(list(n_configurations = replicates$configuration),
                           replicates[group_vars], function(z) length(unique(z)))
summary_table <- merge(merge(summary_table, rep_counts, by = group_vars),
                       config_counts, by = group_vars)
summary_table$coverage <- with(summary_table, ifelse(n_cs > 0, covered_cs / n_cs, NA_real_))
summary_table$purity <- with(summary_table, ifelse(n_cs > 0, purity_sum / n_cs, NA_real_))
summary_table$power <- with(summary_table, recovered / n_causal)

# Every causal allocation at a given total K is pooled. If jobs are incomplete,
# allocations with more completed replicates contribute more observations.
# Save the detailed counts so that imbalance is visible.
configuration_counts <- aggregate(list(n_replicates = replicates$seed),
                                  replicates[c(group_vars, "configuration")], length)
write.csv(configuration_counts, file.path(output_dir, "configuration_counts.csv"), row.names = FALSE)
write.csv(summary_table, file.path(output_dir, "metric_summary.csv"), row.names = FALSE)
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
roc_table <- do.call(rbind, curve_rows)
saveRDS(roc_table, file.path(output_dir, "roc_counts.rds"))
rm(curve_counts, curve_rows)

# ------------------------------------------------------------
# Draw the panels in the style of the supplementary figures
# ------------------------------------------------------------

draw_figure <- function(metric, scenarios, only_K = NULL) {
  nr <- length(scenarios)
  nc <- length(pve_values)
  panels <- matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE)
  layout(cbind(panels, nr * nc + seq_len(nr)), widths = c(rep(1, nc), 0.65))
  par(oma = c(4.5, 3.5, 3, 0.3), mar = c(2.0, 2.0, 1.7, 0.4),
      mgp = c(1.3, 0.4, 0), tcl = -0.2, family = "sans", cex = 0.9)
  is_roc <- metric == "roc"
  # Use the same coverage scale for pure and mixed figures. Start just below
  # the lowest plotted value, rounded down to 0.05; retain the 0.95 reference.
  if (metric == "coverage") {
    values <- summary_table$coverage[
      summary_table$scenario %in% c(pure_rows, mixed_rows) &
        summary_table$pve %in% pve_values
    ]
    values <- values[is.finite(values)]
    coverage_lower <- max(0, floor((min(c(values, target_coverage)) - .01) / .05) * .05)
    coverage_step <- if (1 - coverage_lower <= .30) .05 else .10
    coverage_ticks <- sort(unique(round(c(coverage_lower,
                                         seq(coverage_lower, 1, by = coverage_step), 1), 2)))
  }

  for (i in seq_along(scenarios)) {
    min_K <- if (scenarios[i] %in% pure_rows) 1 else if (scenarios[i] == mixed_rows[4]) 3 else 2
    for (j in seq_along(pve_values)) {
      if (is_roc) {
        xlim <- c(0, roc_max_fpr)
        xticks <- pretty(xlim, n = 5)
        xticks <- xticks[xticks >= 0 & xticks <= roc_max_fpr]
      } else {
        xlim <- c(min_K - 0.35, 5.35)
        xticks <- min_K:5
      }
      ylim <- if (metric == "coverage") c(coverage_lower, 1) else
        if (metric == "purity") c(0.5, 1.015) else c(0, 1.015)
      yticks <- if (metric == "coverage") coverage_ticks else
        if (metric == "purity") seq(.5, 1, .1) else seq(0, 1, .2)
      plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
           axes = FALSE, xlab = "", ylab = "")
      abline(v = xticks, h = yticks, col = "#DEDEDE", lwd = 0.8)
      if (metric == "coverage") abline(h = target_coverage, lty = 2, lwd = 1.2)
      if (is_roc) abline(0, 1, col = "#888888", lty = 2)
      axis(1, at = xticks, cex.axis = 0.9, col = "#777777")
      if (j == 1) axis(2, at = yticks, las = 1, cex.axis = 0.9, col = "#777777")
      if (i == 1) mtext(paste0("PVE = ", round(100 * pve_values[j]), "%"),
                        side = 3, line = 0.55, font = 2)

      d <- if (is_roc) roc_table else summary_table
      d <- d[d$scenario == scenarios[i] & d$pve == pve_values[j], , drop = FALSE]
      if (!is.null(only_K)) d <- d[d$K == only_K, , drop = FALSE]
      if (!nrow(d)) {
        label <- if (!is.null(only_K) && only_K < min_K) "Not applicable" else "No saved results"
        text(mean(xlim), mean(ylim), label, col = "#777777", cex = .8)
        next
      }
      for (m in seq_along(method_names)) {
        dm <- d[d$method == method_names[m], , drop = FALSE]
        if (is_roc) {
          for (k in sort(unique(dm$K))) {
            dk <- dm[dm$K == k, ]
            # Threshold order retains vertical segments and tied-score jumps.
            dk <- dk[order(dk$threshold, decreasing = TRUE), ]
            lines(dk$fpr, dk$tpr, col = method_colors[m],
                  lty = if (is.null(only_K)) k else 1, lwd = 1.5)
          }
        } else {
          points(dm$K + c(-.07, .07)[m], dm[[metric]],
                 pch = 16, cex = 0.95, col = method_colors[m],
                 xpd = metric == "coverage")
        }
      }
    }
  }
  # A narrow column at the right labels each scenario, without crowding panels.
  par(mar = c(0, 0, 0, 0))
  for (s in scenarios) {
    plot.new()
    label <- gsub(" + ", "\n+ ", s, fixed = TRUE)
    text(.05, .5, label, adj = c(0, .5), font = 2, cex = .85)
  }
  titles <- c(coverage = "Credible-set coverage", purity = "Credible-set purity",
              power = "Causal SNP recovery by credible sets", roc = "Detection of causal SNPs using PIPs")
  title_text <- paste0(titles[metric], "  |  n = ", n_value)
  if (!is.null(only_K)) title_text <- paste0(title_text, "  |  ", only_K, " causal SNPs")
  mtext(title_text, side = 3, outer = TRUE, line = 1.2, font = 2, cex = 1.1)
  mtext(if (is_roc) "False positive rate" else "Number of causal SNPs",
        side = 1, outer = TRUE, line = 0.5)
  ylab <- c(coverage = "Coverage", purity = "Mean minimum absolute correlation",
            power = "Power", roc = "True positive rate (power)")
  mtext(ylab[metric], side = 2, outer = TRUE, line = 1.8)

  # Draw a common legend in the outer bottom margin.
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend(.5, .045, legend = method_names, col = method_colors,
         pch = if (is_roc) NA else 16, lty = if (is_roc) 1 else NA,
         lwd = 1.5, horiz = TRUE, xjust = .5, yjust = .5, bty = "n", cex = .95)
  if (is_roc && is.null(only_K)) {
    ks <- if (all(scenarios %in% pure_rows)) 1:5 else 2:5
    legend(.5, .015, legend = paste("L =", ks), lty = ks, horiz = TRUE,
           xjust = .5, yjust = .5, bty = "n", cex = .85, seg.len = 2.8)
  }
}

save_figure <- function(metric, scenarios, name) {
  height <- 2.0 * length(scenarios) + 1.4
  pdf(file.path(output_dir, paste0(name, ".pdf")), width = 13.5, height = height,
      useDingbats = FALSE)
  tryCatch(draw_figure(metric, scenarios), finally = dev.off())
  if (write_png) {
    png(file.path(output_dir, paste0(name, ".png")), width = 13.5, height = height,
        units = "in", res = 180)
    tryCatch(draw_figure(metric, scenarios), finally = dev.off())
  }
}

for (metric in c("coverage", "purity", "power", "roc")) {
  save_figure(metric, pure_rows, paste0(metric, "_pure"))
  save_figure(metric, mixed_rows, paste0(metric, "_mixed"))
}

if (write_roc_pages) {
  for (group in c("pure", "mixed")) {
    scenarios <- if (group == "pure") pure_rows else mixed_rows
    ks <- if (group == "pure") 1:5 else 2:5
    pdf(file.path(output_dir, paste0("roc_", group, "_by_L.pdf")),
        width = 13.5, height = 2.0 * length(scenarios) + 1.4, useDingbats = FALSE)
    tryCatch(for (k in ks) draw_figure("roc", scenarios, only_K = k), finally = dev.off())
  }
}

cat("\nFigures and summary tables saved in:", output_dir, "\n")
cat("Included", sum(audit$included), "unique simulations; skipped",
    sum(audit$duplicates), "duplicate copies; recorded", sum(audit$errors), "error entries.\n")
cat("Nonconverged replicate pairs:", sum(audit$nonconverged),
    "; excluded:", sum(audit$excluded_nonconverged), "\n")
cat("L in the ROC legend is the TRUE causal count; fitted SuSiE L =", fit_L, ".\n")
cat("Check configuration_counts.csv for incomplete or unbalanced conditions.\n")
