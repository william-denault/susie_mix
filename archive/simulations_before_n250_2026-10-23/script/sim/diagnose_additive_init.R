# Read-only analysis of saved simulation checkpoints. No fitting or GTEx access.
# source(.../diagnose_additive_init.R); diagnostics <- diagnose_additive_initialization()
.additive_diag_script_dir <- local({
  files <- vapply(sys.frames(), function(f) if (is.null(f$ofile)) "" else f$ofile, "")
  files <- files[!is.na(files) & nzchar(files)]
  if (length(files)) dirname(normalizePath(tail(files, 1L), winslash = "/")) else
    file.path(Sys.getenv("SUSIE_MIX_PROJECT_DIR", getwd()), "script/sim")
})

# Exact within-dataset AUC: probability a causal PIP exceeds a noncausal PIP,
# with half credit for ties. This is distinct from the pooled ROC across loci.
additive_diag_auc <- function(pip, truth) {
  k <- length(truth)
  (sum(rank(pip, ties.method = "average")[truth]) - k * (k + 1) / 2) /
    (k * (length(pip) - k))
}

# Area under the SAME threshold-grid ROC as plot_additive_init_roc.R, divided
# by max_fpr (mean TPR over that interval). Linear interpolation matches its
# plotted segments; this is not an exact all-distinct-PIP partial AUC.
additive_diag_pauc <- function(tp, fp, positives, negatives, max_fpr) {
  if (positives <= 0 || negatives <= 0) return(NA_real_)
  x <- rev(fp) / negatives; y <- rev(tp) / positives
  dx <- diff(x)
  keep <- which(dx > 0 & head(x, -1L) < max_fpr)
  width <- pmin(x[keep + 1L], max_fpr) - x[keep]
  end_y <- y[keep] + (y[keep + 1L] - y[keep]) * width / dx[keep]
  sum(width * (y[keep] + end_y) / 2) / max_fpr
}

diagnose_additive_initialization <- function(
    results_dir = file.path(Sys.getenv("SUSIE_MIX_PROJECT_DIR",
      "/project2/mstephens/wdenault/susie_mix"), "simulation results/additive_slide_init_v3"),
    output_dir = file.path(results_dir, "diagnostics"), pves = c(.025, .05),
    max_fpr = .25, top_examples = 3L, exclude_nonconverged = FALSE) {
  stopifnot(length(max_fpr) == 1L, is.finite(max_fpr), max_fpr > 0, max_fpr <= 1,
    length(top_examples) == 1L, is.finite(top_examples), top_examples >= 0,
    top_examples == floor(top_examples))
  helpers <- new.env(parent = environment())
  old <- options(susie.init.roc.autorun = FALSE)
  on.exit(options(old), add = TRUE)
  source(file.path(.additive_diag_script_dir, "plot_additive_init_roc.R"), local = helpers)
  source(file.path(.additive_diag_script_dir, "simulation_metric_helpers.R"), local = helpers)
  # Reuse the existing checkpoint validation, pairing, convergence policy and
  # threshold grid. These outputs go only to the new diagnostics directory.
  roc <- helpers$plot_additive_init_roc(results_dir, output_dir, pves = pves,
    zoom_max_fpr = max_fpr, exclude_nonconverged = exclude_nonconverged, write_png = FALSE)
  audit <- roc$replicates
  included <- audit[grepl("^included", audit$status), , drop = FALSE]
  methods <- helpers$additive_roc_methods
  colors <- helpers$additive_roc_colors
  rows <- causal_rows <- metric_rows <- blocks <- list()
  scalar <- function(x, fallback = NA) if (length(x)) x[1L] else fallback
  for (file in unique(included$file)) {
    saved <- readRDS(file)
    for (slot in included$slot[included$file == file]) {
      r <- saved$results[[slot]]; pve <- saved$settings$pve
      truth <- r$true_pos; k <- length(truth)
      if (r$p <= k) stop("Within-dataset AUC requires noncausal SNPs; seed ", r$seed)
      if (!all(vapply(r$fits[methods], function(f) "sets" %in% names(f), TRUE)))
        stop("Missing saved credible-set information; seed ", r$seed)
      pip <- lapply(r$fits[methods], `[[`, "pip")
      auc <- vapply(pip, additive_diag_auc, 0, truth = truth)
      ranks <- lapply(pip, function(x) rank(-x, ties.method = "min"))
      # Recompute CS metrics from saved membership, using the original helper.
      cs <- t(vapply(r$fits[methods], function(f) helpers$cs_summary(f$sets, f$sets$cs, truth),
        c(n_cs = 0, covered_cs = 0, purity_sum = 0, cs_size_sum = 0,
          cs_size_sum_sq = 0, recovered = 0, n_causal = 0)))
      metric_rows[[length(metric_rows) + 1L]] <- data.frame(scenario = "Additive only",
        pve = pve, K = k, configuration = "add2", seed = r$seed, method = methods,
        converged = vapply(r$fits[methods], function(f) isTRUE(f$converged), TRUE), cs, row.names = NULL)
      counts <- lapply(pip, helpers$additive_roc_counts, truth = truth)
      blocks[[length(blocks) + 1L]] <- list(pve = pve, seed = r$seed, k = k,
        nulls = r$p - k, counts = counts)
      for (m in 2:3) {
        d <- pip[[m]] - pip[[1L]]
        rows[[length(rows) + 1L]] <- data.frame(pve = pve, seed = r$seed,
          gene = scalar(r$gene, ""), p = r$p, causal_r = scalar(r$causal_r),
          method = methods[m], reference = methods[1L],
          mean_abs_pip_change = mean(abs(d)), q95_abs_pip_change = unname(quantile(abs(d), .95)),
          max_abs_pip_change = max(abs(d)), fraction_pip_changes_over_0.01 = mean(abs(d) > .01),
          mean_causal_pip_change = mean(d[truth]), mean_null_pip_change = mean(d[-truth]),
          reference_within_auc = auc[1L], method_within_auc = auc[m], within_auc_gain = auc[m] - auc[1L],
          mean_causal_rank_gain = mean(ranks[[1L]][truth] - ranks[[m]][truth]),
          reference_recovered = cs[1L, "recovered"], method_recovered = cs[m, "recovered"],
          extra_recovered = cs[m, "recovered"] - cs[1L, "recovered"],
          reference_false_cs = cs[1L, "n_cs"] - cs[1L, "covered_cs"],
          method_false_cs = cs[m, "n_cs"] - cs[m, "covered_cs"],
          all_converged = all(vapply(r$fits[methods], function(f) isTRUE(f$converged), TRUE)),
          file = file, slot = slot, row.names = NULL)
      }
      for (m in seq_along(methods)) {
        in_cs <- truth %in% unlist(r$fits[[methods[m]]]$sets$cs, use.names = FALSE)
        causal_rows[[length(causal_rows) + 1L]] <- data.frame(pve = pve, seed = r$seed,
          gene = scalar(r$gene, ""), snp_index = truth,
          snp = if (length(r$causal_snps) == k) r$causal_snps else as.character(truth),
          beta_standardized = if (length(r$beta_standardized) == k) r$beta_standardized else NA_real_,
          method = methods[m], pip = pip[[m]][truth], rank = ranks[[m]][truth],
          noncausal_fraction_above = vapply(pip[[m]][truth], function(x) mean(pip[[m]][-truth] > x), 0),
          recovered_in_cs = in_cs, row.names = NULL)
      }
    }
  }
  datasets <- do.call(rbind, rows); causal <- do.call(rbind, causal_rows)
  metrics <- do.call(rbind, metric_rows)
  cs_summary <- helpers$summarize_metrics(metrics, methods = methods)
  cs_summary$summary$false_cs <- cs_summary$summary$n_cs - cs_summary$summary$covered_cs
  cs_summary$summary$cs_fdr <- 1 - cs_summary$summary$coverage
  # Leave one dataset out of BOTH methods, with denominators recalculated.
  # Influence is the change in the pooled ROC gap, not an additive attribution.
  influence_rows <- sensitivity_rows <- overview_rows <- list()
  curves <- list()
  for (pve in sort(unique(datasets$pve))) {
    b <- blocks[vapply(blocks, function(x) x$pve == pve, TRUE)]
    n <- length(b); seeds <- vapply(b, `[[`, 0, "seed")
    positives <- sum(vapply(b, `[[`, 0, "k")); negatives <- sum(vapply(b, `[[`, 0, "nulls"))
    totals <- lapply(seq_along(methods), function(m) Reduce(`+`, lapply(b, function(x) x$counts[[m]])))
    area <- function(counts, pos, neg) vapply(counts, function(x)
      additive_diag_pauc(x[, "tp"], x[, "fp"], pos, neg, max_fpr), 0)
    full <- area(totals, positives, negatives)
    loo <- matrix(NA_real_, n, 3L)
    if (n > 1L) for (i in seq_len(n)) {
      left <- lapply(seq_along(methods), function(m) totals[[m]] - b[[i]]$counts[[m]])
      loo[i, ] <- area(left, positives - b[[i]]$k, negatives - b[[i]]$nulls)
    }
    for (m in 2:3) {
      gap <- full[m] - full[1L]
      influence_rows[[length(influence_rows) + 1L]] <- data.frame(pve = pve, seed = seeds,
        method = methods[m], max_fpr = max_fpr, full_gap = gap,
        gap_without_dataset = loo[, m] - loo[, 1L],
        influence = gap - (loo[, m] - loo[, 1L]))
      d <- datasets[datasets$pve == pve & datasets$method == methods[m], ]
      overview_rows[[length(overview_rows) + 1L]] <- data.frame(pve = pve, method = methods[m],
        n_replicates = n, max_fpr = max_fpr, pooled_pauc_gain = gap,
        mean_within_auc_gain = mean(d$within_auc_gain), median_within_auc_gain = median(d$within_auc_gain),
        n_within_auc_better = sum(d$within_auc_gain > 1e-12),
        n_within_auc_worse = sum(d$within_auc_gain < -1e-12),
        n_within_auc_tied = sum(abs(d$within_auc_gain) <= 1e-12),
        median_mean_abs_pip_change = median(d$mean_abs_pip_change),
        median_max_abs_pip_change = median(d$max_abs_pip_change),
        n_cs_recovery_better = sum(d$extra_recovered > 0), n_cs_recovery_worse = sum(d$extra_recovered < 0))
    }
    slide <- datasets[datasets$pve == pve & datasets$method == methods[2L], ]
    slide <- slide[match(seeds, slide$seed), ]
    order_influence <- order((full[2L] - full[1L]) - (loo[, 2L] - loo[, 1L]), decreasing = TRUE, na.last = TRUE)
    order_pip <- order(slide$max_abs_pip_change, decreasing = TRUE)
    subsets <- list(all = integer())
    for (count in c(1L, 5L, 10L)) if (count < n) {
      subsets[[paste0("remove_top_influence_", count)]] <- order_influence[seq_len(count)]
      subsets[[paste0("remove_largest_pip_change_", count)]] <- order_pip[seq_len(count)]
    }
    for (label in names(subsets)) {
      remove <- subsets[[label]]; keep <- setdiff(seq_len(n), remove)
      counts <- lapply(seq_along(methods), function(m) Reduce(`+`, lapply(b[keep], function(x) x$counts[[m]])))
      pos <- sum(vapply(b[keep], `[[`, 0, "k")); neg <- sum(vapply(b[keep], `[[`, 0, "nulls"))
      a <- area(counts, pos, neg)
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(pve = pve, subset = label,
        n_removed = length(remove), n_replicates = length(keep), max_fpr = max_fpr,
        method = methods, mean_tpr_over_fpr_interval = a, gain_over_susie = a - a[1L],
        removed_seeds = paste(seeds[remove], collapse = ";"))
      for (m in seq_along(methods)) curves[[length(curves) + 1L]] <- data.frame(pve = pve,
        subset = label, method = methods[m], threshold = c(helpers$additive_roc_thresholds, Inf),
        tp = counts[[m]][, "tp"], fp = counts[[m]][, "fp"], n_causal = pos, n_noncausal = neg,
        tpr = counts[[m]][, "tp"] / pos, fpr = counts[[m]][, "fp"] / neg)
    }
  }
  influence <- do.call(rbind, influence_rows)
  sensitivity <- do.call(rbind, sensitivity_rows); overview <- do.call(rbind, overview_rows)
  curves <- do.call(rbind, curves)
  key <- function(x) paste(x$pve, x$seed, x$method, sep = "|")
  datasets$pooled_roc_influence <- influence$influence[match(key(datasets), key(influence))]
  # Include clear CS wins, ranking wins, influential datasets, and reverse
  # examples. Selection is exploratory and its rule is recorded explicitly.
  selections <- list()
  for (pve in unique(datasets$pve)) {
    d <- datasets[datasets$pve == pve & datasets$method == methods[2L], ]
    candidates <- list(
      cs_win_no_extra_false_cs = which(d$extra_recovered > 0 & d$method_false_cs <= d$reference_false_cs),
      within_dataset_ranking_win = which(d$within_auc_gain > 1e-12),
      largest_positive_roc_influence = which(is.finite(d$pooled_roc_influence) & d$pooled_roc_influence > 0),
      within_dataset_ranking_loss = which(d$within_auc_gain < -1e-12))
    scores <- list(d$extra_recovered, d$within_auc_gain, d$pooled_roc_influence, -d$within_auc_gain)
    for (j in seq_along(candidates)) {
      ids <- candidates[[j]]; ids <- head(ids[order(scores[[j]][ids], decreasing = TRUE)], top_examples)
      if (length(ids)) selections[[length(selections) + 1L]] <- data.frame(selection = names(candidates)[j], d[ids, ], row.names = NULL)
    }
  }
  examples <- if (length(selections)) do.call(rbind, selections) else
    data.frame(selection = character(), datasets[FALSE, ], row.names = NULL)
  output <- list(overview = overview, per_dataset = datasets, causal_snps = causal,
    cs_summary = cs_summary$summary, cs_paired_differences = cs_summary$differences,
    roc_influence = influence, roc_sensitivity = sensitivity, roc_sensitivity_curves = curves,
    examples = examples)
  for (name in names(output)) write.csv(output[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  additive_diag_plots(output, methods, colors, output_dir, max_fpr)
  writeLines(c(
    "Saved-fit diagnostics: SuSiE versus continuous SuSiE-slide (and initialized additive SuSiE).",
    "No fits, package versions, input files, or checkpoint contents were changed.",
    paste0("ROC area is divided by max_fpr = ", max_fpr, "; it is mean TPR over that FPR interval."),
    "It uses trapezoids on the original ROC threshold grid, not an exact all-distinct-PIP curve.",
    "Within-dataset AUC uses all distinct PIPs with half credit for ties; pooled ROC also compares SNPs across datasets.",
    "Positive influence means omitting that dataset reduces the slide-minus-SuSiE pooled ROC gap.",
    "Influences do not sum to the total gap. Removal rankings are fixed from the full data, not iteratively reselected.",
    "Sensitivity removes the same datasets from every method and recomputes denominators and curves.",
    "CS power = distinct causal SNPs in any reported CS / all causal SNPs.",
    "CS coverage = reported CSs containing at least one causal SNP / all reported CSs; CS FDR = 1 - coverage.",
    "CS metrics use the saved 95% CSs and purity filtering; no-CS coverage/FDR is NA.",
    "ROC improvement concerns discrimination; it does not establish calibrated PIPs or nominal CS coverage.",
    "Examples are selected exploratory cases, including losses, not independent evidence of superiority.",
    "Example SNP axes are retained-column indices, not genomic distance or LD.",
    "Full X/y and component alpha are not saved by default; mechanistic LD/coding attribution needs a targeted rerun.",
    "Successful nonconverged fits follow the same inclusion policy as the existing ROC script; see roc_pip_replicates.csv."
  ), file.path(output_dir, "README.txt"))
  print(overview, row.names = FALSE)
  print(cs_summary$summary[c("pve", "method", "n_replicates", "power", "coverage", "cs_fdr", "cs_size")], row.names = FALSE)
  message("Diagnostics saved in: ", normalizePath(output_dir, winslash = "/"))
  invisible(output)
}

additive_diag_plots <- function(output, methods, colors, output_dir, max_fpr) {
  grDevices::pdf(file.path(output_dir, "diagnostics.pdf"), width = 10, height = 8, useDingbats = FALSE)
  tryCatch({
    for (pve in unique(output$per_dataset$pve)) {
      par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(1, 0, 2, 0))
      d <- output$per_dataset[output$per_dataset$pve == pve & output$per_dataset$method == methods[2L], ]
      hist(d$mean_abs_pip_change, main = "Magnitude of PIP changes", xlab = "Mean |slide PIP - SuSiE PIP| per dataset", col = colors[2L])
      hist(d$within_auc_gain, main = "Within-dataset ranking changes", xlab = "Slide - SuSiE AUC (exact; ties get half credit)", col = colors[2L])
      abline(v = 0, lty = 2)
      if (any(is.finite(d$pooled_roc_influence))) {
        plot(d$max_abs_pip_change, d$pooled_roc_influence, pch = 16, col = adjustcolor(colors[2L], .6),
          xlab = "Largest absolute PIP change in dataset", ylab = "Pooled ROC gap - gap without dataset",
          main = "Large PIP changes versus ROC influence")
        abline(h = 0, lty = 2)
      } else { plot.new(); title("Influence needs at least two datasets") }
      plot(NA, xlim = c(0, max_fpr), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate",
        main = "Pooled ROC: omit influential datasets")
      labels <- c("all", "remove_top_influence_5")
      if (!any(output$roc_sensitivity_curves$pve == pve & output$roc_sensitivity_curves$subset == labels[2L]))
        labels <- c("all", "remove_top_influence_1")
      labels <- labels[labels %in% output$roc_sensitivity_curves$subset[output$roc_sensitivity_curves$pve == pve]]
      for (j in seq_along(labels)) for (m in 1:2) {
        c <- output$roc_sensitivity_curves
        c <- c[c$pve == pve & c$subset == labels[j] & c$method == methods[m], ]
        if (nrow(c)) lines(rev(c$fpr), rev(c$tpr), col = colors[m], lty = j, lwd = 2)
      }
      subset_labels <- c("All datasets", if (length(labels) > 1L) gsub("_", " ", labels[2L]))
      legend("bottomright", legend = c(methods[1:2], subset_labels),
        col = c(colors[1:2], rep("black", length(labels))), lty = c(1, 1, seq_along(labels)), bty = "n", cex = .75)
      mtext(sprintf("PVE %g%% | %d datasets | continuous slide versus additive SuSiE", 100 * pve, nrow(d)), outer = TRUE, font = 2)
    }
  }, finally = grDevices::dev.off())
  examples <- output$examples
  examples <- examples[!duplicated(paste(examples$pve, examples$seed)), ]
  grDevices::pdf(file.path(output_dir, "examples.pdf"), width = 11, height = 8.5, useDingbats = FALSE)
  tryCatch({
    if (!nrow(examples)) { plot.new(); text(.5, .5, "No examples selected by the requested rules.") }
    cached_file <- NULL; saved <- NULL
    for (i in seq_len(nrow(examples))) {
      e <- examples[i, ]
      if (!identical(cached_file, e$file)) { saved <- readRDS(e$file); cached_file <- e$file }
      r <- saved$results[[e$slot]]; truth <- r$true_pos
      pips <- sapply(r$fits[methods], `[[`, "pip")
      par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1), oma = c(2, 0, 3, 0))
      plot(pips[, 1L], pips[, 2L], pch = 16, cex = .5, col = "#99999966", xlim = c(0, 1), ylim = c(0, 1),
        xlab = "SuSiE PIP", ylab = "SuSiE-slide PIP", main = "Causal SNPs highlighted")
      abline(0, 1, lty = 2)
      points(pips[truth, 1L], pips[truth, 2L], pch = 21, bg = "#E69F00", cex = 1.3)
      text(pips[truth, 1L], pips[truth, 2L], labels = seq_along(truth), pos = 3, cex = .8)
      matplot(seq_len(r$p), pips, type = "p", pch = c(1, 3, 4), cex = .5, col = colors, ylim = c(0, 1),
        xlab = "Retained SNP index (not genomic distance)", ylab = "PIP", main = "Same dataset, all three methods")
      abline(v = truth, lty = 3, col = "#E69F00")
      legend("topright", methods, col = colors, pch = c(1, 3, 4), bty = "n", cex = .7)
      plot.new(); title("Causal SNP PIPs and ranks")
      cr <- output$causal_snps
      cr <- cr[cr$pve == e$pve & cr$seed == e$seed, ]
      lines <- unlist(lapply(seq_along(truth), function(j) {
        z <- cr[cr$snp_index == truth[j], ]
        c(sprintf("Causal %d: column %d", j, truth[j]),
          sprintf("%s: PIP %.4g, rank %d, in CS %s", z$method, z$pip, z$rank, ifelse(z$recovered_in_cs, "yes", "no")), "")
      }))
      text(0, 1, paste(lines, collapse = "\n"), adj = c(0, 1), cex = .75)
      plot.new(); title("Credible sets and paired differences")
      lines <- vapply(methods, function(m) {
        sets <- r$fits[[m]]$sets$cs
        sprintf("%s: %d CSs; recovered %d/%d; false CSs %d\nsizes: %s", m, length(sets),
          sum(truth %in% unlist(sets)), length(truth), sum(!vapply(sets, function(s) any(s %in% truth), TRUE)),
          if (length(sets)) paste(lengths(sets), collapse = ", ") else "none")
      }, "")
      text(0, 1, paste(c(lines, "", sprintf("Within-dataset AUC gain: %+.4f", e$within_auc_gain),
        sprintf("Mean / max |PIP change|: %.4g / %.4g", e$mean_abs_pip_change, e$max_abs_pip_change),
        sprintf("Pooled ROC influence: %+.4g", e$pooled_roc_influence),
        sprintf("Correlation of causal genotypes: %.3f", e$causal_r)), collapse = "\n"), adj = c(0, 1), cex = .78)
      mtext(sprintf("PVE %g%% | seed %s | gene %s", 100 * e$pve, e$seed, e$gene), outer = TRUE, line = 1, font = 2)
      mtext(paste("Selected example:", gsub("_", " ", e$selection)), outer = TRUE, cex = .85)
      mtext("Exploratory selection; causal SNP identities and coefficients are in causal_snps.csv", side = 1, outer = TRUE, cex = .8)
    }
  }, finally = grDevices::dev.off())
}
