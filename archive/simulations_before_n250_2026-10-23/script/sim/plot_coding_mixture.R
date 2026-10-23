# Coding identification and mixture recovery from existing simulation saves.
# In R/RStudio: source(".../script/sim/plot_coding_mixture.R")
# No simulations or SuSiE fits are run. Only base R is required.

project_dir <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (.Platform$OS.type == "windows")
  "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix" else
  "/project2/mstephens/wdenault/susie_mix")
source(file.path(project_dir, "script/sim/coding_mixture_utils.R"))

coding_plot_settings <- list(
  chunk_dir = file.path(project_dir, "simulation results/chunks"),
  output_dir = file.path(project_dir, "simulation results/coding_figures"),
  file_pattern = "\\.RData$",
  n = 500, fit_L = 10, pve_values = c(.1, .2, .3, .4),
  snp_cutoff = .9,              # SNP detection for classification; also baseline call cutoff.
  thresholds = c(0, .1, .5, .8, .9, .95, .99, 1),
  n_bins = 10L,
  exclude_nonconverged = FALSE, # Only mixed-fit convergence matters here.
  max_reps_per_file = Inf,      # e.g. 5 for a quick preview; full checkpoints are still loaded.
  bootstrap_reps = 500L,
  bootstrap_seed = 20260911L,
  write_png = TRUE
)
cm_colors <- c(additive = "#2474B5", recessive = "#BB3E70", dominant = "#D77916")
cm_labels <- c(additive = "Additive", recessive = "Recessive", dominant = "Dominant")

cm_panel <- function(xlim = c(.7, 5.3), ylim = c(0, 1), xlab = "", ylab = "", title = "") {
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "",
       main = title, bty = "l", las = 1, xaxs = "i", yaxs = "i", cex.main = .95)
  if (nzchar(ylab)) mtext(ylab, side = 2, line = 2.8)
  abline(h = pretty(ylim), col = "#EEEEEE")
}

cm_start_layout <- function(nr, nc, title, subtitle = "") {
  par(mfrow = c(nr, nc), mar = c(3.4, 3.7, 2.4, .9),
      oma = c(2.5, .3, 4.2, .3), mgp = c(2.1, .65, 0), tcl = -.25,
      family = "sans", cex = .9)
  list(title = title, subtitle = subtitle)
}

cm_finish_layout <- function(heading) {
  mtext(heading$title, side = 3, outer = TRUE, line = 2.6, font = 2, cex = 1.25)
  mtext(heading$subtitle, side = 3, outer = TRUE, line = 1.0, cex = .85)
}

cm_legend <- function(classes = cm_classes, note = "") {
  # Explicit colored labels remain readable in both exported formats.
  xs <- seq(.25, .75, length.out = length(classes))
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  for (i in seq_along(classes)) text(xs[i], .018, cm_labels[classes[i]], col = cm_colors[classes[i]], font = 2, cex = .9)
  if (nzchar(note)) mtext(note, side = 1, line = -2, cex = .7)
}

cm_draw_baseline <- function(analysis, config) {
  d <- analysis$mixture[analysis$mixture$scenario == "Additive only", ]
  heading <- cm_start_layout(1, length(config$pve_values), "PIP allocation when every true effect is additive",
                  "Same pooled-PIP ratio used by the EM update; vertical bars: 95% seed-block bootstrap intervals")
  for (pve in config$pve_values) {
    cm_panel(xlab = "Number of true causal SNPs", ylab = "Share of all coding PIP mass",
             title = paste0("PVE = ", 100 * pve, "%"))
    abline(h = c(0, 1), lty = 2, col = "#777777")
    if (!any(d$pve == pve)) { text(3, .5, "No saved additive-only runs"); next }
    for (j in seq_along(cm_classes)) {
      x <- d[d$pve == pve & d$coding == cm_classes[j], ]
      x <- x[order(x$K), ]
      xpos <- x$K + (j - 2) * .06
      segments(xpos, x$lower, xpos, x$upper, col = cm_colors[j])
      # Separate controls with different causal-SNP selection designs in the CSV;
      # do not connect them as though they were one experimental series.
      for (control in unique(x$all_additive)) {
        ids <- which(x$all_additive == control)
        if (!anyDuplicated(x$K[ids])) lines(xpos[ids], x$estimated_share[ids], col = cm_colors[j])
      }
      points(xpos, x$estimated_share, pch = 16, col = cm_colors[j])
    }
  }
  cm_finish_layout(heading)
  cm_legend()
}

cm_draw_recovery <- function(analysis, config) {
  heading <- cm_start_layout(3, length(config$pve_values), "Recovery of the generating coding mixture",
                  "One point per causal allocation; pooled coding PIP share versus fraction of true effects; dashed line: equality")
  for (cl in cm_classes) for (pve in config$pve_values) {
    cm_panel(c(-.03, 1.03), c(0, 1), "True fraction of effects", paste(cm_labels[cl], "PIP share"),
             paste0(cm_labels[cl], "  |  PVE = ", 100 * pve, "%"))
    abline(0, 1, lty = 2, col = "#555555")
    d <- analysis$mixture[analysis$mixture$coding == cl & analysis$mixture$pve == pve, ]
    if (!nrow(d)) { text(.5, .5, "No saved results"); next }
    segments(d$true_share, d$lower, d$true_share, d$upper, col = adjustcolor(cm_colors[cl], .4))
    points(d$true_share, d$estimated_share, col = cm_colors[cl], pch = c(16, 17, 15, 18, 8)[pmin(d$K, 5)], cex = .8)
  }
  cm_finish_layout(heading)
  mtext("Symbols: K=1 circle; K=2 triangle; K=3 square; K=4 diamond; K=5 star. 95% seed-block intervals.",
        side = 1, outer = TRUE, line = .5, cex = .8)
}

cm_draw_confusion <- function(analysis, config) {
  heading <- cm_start_layout(1, length(config$pve_values), "Coding calls at the known causal SNP",
                  paste0("All causal SNPs in denominator; SNP PIP >= ", config$snp_cutoff,
                         " required; exact PIP ties are ambiguous; all K/configurations pooled"))
  for (pve in config$pve_values) {
    d <- analysis$confusion[analysis$confusion$pve == pve, ]
    counts <- matrix(0, 3, 5, dimnames = list(cm_classes, cm_calls))
    for (i in seq_len(nrow(d))) counts[d$true_coding[i], d$called_coding[i]] <-
      counts[d$true_coding[i], d$called_coding[i]] + d$count[i]
    rates <- counts / rowSums(counts)
    plot(NA, xlim = c(.5, 5.5), ylim = c(.5, 3.5), axes = FALSE,
         xlab = "Called coding", ylab = "", main = paste0("PVE = ", 100 * pve, "%"),
         xaxs = "i", yaxs = "i")
    palette <- colorRampPalette(c("#F5F8FC", "#135890"))(101)
    for (i in 1:3) for (j in 1:5) {
      v <- rates[i, j]
      fill <- if (is.finite(v)) palette[1 + round(v * 100)] else "#EEEEEE"
      rect(j - .49, 4 - i - .49, j + .49, 4 - i + .49, col = fill, border = "white")
      text(j, 4 - i, if (is.finite(v)) sprintf("%.0f%%", 100 * v) else "NA",
           col = if (is.finite(v) && v > .6) "white" else "#222222", cex = .9)
    }
    axis(1, at = 1:5, labels = c("Add", "Rec", "Dom", "Tie", "Not\ndetected"), tick = FALSE, cex.axis = .8)
    axis(2, at = 3:1, labels = c("Add", "Rec", "Dom"), las = 1, tick = FALSE)
    mtext("True coding", side = 2, line = 2.75)
  }
  cm_finish_layout(heading)
  mtext("Location is supplied by the simulation truth. This figure does not measure false calls at other SNPs.",
        side = 1, outer = TRUE, line = .6, cex = .85)
}

cm_classification_rates <- function(confusion) {
  groups <- split(seq_len(nrow(confusion)), interaction(confusion[c("pve", "K", "true_coding")], drop = TRUE))
  do.call(rbind, lapply(groups, function(ids) {
    d <- confusion[ids, ]
    n <- sum(d$count)
    detected <- sum(d$count[d$called_coding != "undetected"])
    correct <- sum(d$count[d$called_coding == d$true_coding])
    data.frame(d[1, c("pve", "K", "true_coding")], n_true = n, n_detected = detected,
               n_correct = correct, overall = cm_ratio(correct, n),
               conditional = cm_ratio(correct, detected), row.names = NULL)
  }))
}

cm_draw_classification <- function(analysis, config) {
  heading <- cm_start_layout(2, length(config$pve_values), "How often is the causal coding identified?",
                  paste0("At the known causal SNP, using the largest coding PIP; detection requires SNP PIP >= ", config$snp_cutoff))
  for (metric in c("overall", "conditional")) for (pve in config$pve_values) {
    cm_panel(xlab = "Number of true causal SNPs",
             ylab = if (metric == "overall") "Correct / all causal SNPs" else "Correct / detected causal SNPs",
             title = paste0("PVE = ", 100 * pve, "%"))
    d <- analysis$classification[analysis$classification$pve == pve, ]
    for (cl in cm_classes) {
      x <- d[d$true_coding == cl, ]; x <- x[order(x$K), ]
      lines(x$K, x[[metric]], type = "b", pch = 16, col = cm_colors[cl])
    }
  }
  cm_finish_layout(heading)
  cm_legend()
}

cm_draw_calibration <- function(analysis, config, additive_only = FALSE) {
  raw <- analysis$calibration
  if (additive_only) raw <- raw[raw$scenario == "Additive only", ]
  d <- if (nrow(raw)) cm_calibration_rates(raw, c("pve", "coding", "bin")) else raw
  heading <- cm_start_layout(1, length(config$pve_values),
    if (additive_only) "Coding-specific PIP calibration under additive truth" else "Coding-specific PIP calibration across all simulations",
    "A positive label requires the exact generating SNP AND coding; points pool predictors by PIP bin; dashed line: equality")
  for (pve in config$pve_values) {
    cm_panel(c(0, 1), c(0, 1), "Mean predicted PIP in bin", "Observed fraction of true predictors",
             paste0("PVE = ", 100 * pve, "%"))
    abline(0, 1, lty = 2, col = "#555555")
    for (cl in cm_classes) {
      x <- d[d$pve == pve & d$coding == cl & d$n > 0, ]; x <- x[order(x$bin), ]
      lines(x$mean_pip, x$observed_frequency, col = cm_colors[cl])
      points(x$mean_pip, x$observed_frequency, pch = 16, col = cm_colors[cl],
             cex = .55 + .6 * pmin(1, log10(pmax(x$n, 1)) / 5), xpd = NA)
    }
  }
  cm_finish_layout(heading)
  cm_legend()
}

cm_draw_discoveries <- function(analysis, config) {
  heading <- cm_start_layout(2, length(config$pve_values), paste0("Exact coding discovery at PIP >= ", config$snp_cutoff),
                  "All predictors tested; a wrong coding at a true SNP is a false discovery; all configurations pooled")
  d <- analysis$discovery[analysis$discovery$threshold == config$snp_cutoff, ]
  d <- aggregate(d[c("n_calls", "n_correct", "n_true")], d[c("pve", "K", "coding")], sum)
  d$precision <- cm_ratio(d$n_correct, d$n_calls)
  d$recall <- cm_ratio(d$n_correct, d$n_true)
  for (metric in c("precision", "recall")) for (pve in config$pve_values) {
    cm_panel(xlab = "Number of true causal SNPs", ylab = if (metric == "precision")
      "Correct / all calls of this coding" else "Recovered / true effects of this coding",
      title = paste0("PVE = ", 100 * pve, "%"))
    for (cl in cm_classes) {
      x <- d[d$pve == pve & d$coding == cl, ]; x <- x[order(x$K), ]
      lines(x$K, x[[metric]], type = "b", pch = 16, col = cm_colors[cl])
    }
  }
  cm_finish_layout(heading)
  cm_legend()
}

cm_draw_false_calls <- function(analysis, config) {
  d <- analysis$replicates[analysis$replicates$scenario == "Additive only", ]
  rates <- if (nrow(d)) aggregate(d[c("any_false_rec", "any_false_dom")], d[c("pve", "K")], mean) else NULL
  ymax <- if (is.null(rates)) .02 else max(.01, ceiling(110 * max(rates[c("any_false_rec", "any_false_dom")])) / 100)
  heading <- cm_start_layout(1, length(config$pve_values), "Confident nonadditive calls under purely additive truth",
                  paste0("Fraction of simulations with at least one false coding-specific call at PIP >= ", config$snp_cutoff,
                         "; every SNP in the locus is tested"))
  for (pve in config$pve_values) {
    cm_panel(ylim = c(0, ymax), xlab = "Number of true causal SNPs", ylab = "Simulations with >= 1 false call",
             title = paste0("PVE = ", 100 * pve, "%"))
    for (j in 2:3) {
      x <- d[d$pve == pve, ]; if (!nrow(x)) next
      points_by_K <- aggregate(x[[paste0("any_false_", cm_short[j])]], list(K = x$K), mean)
      lines(points_by_K$K, points_by_K$x, type = "b", pch = 16, col = cm_colors[j], xpd = NA)
    }
  }
  cm_finish_layout(heading)
  cm_legend(cm_classes[2:3])
}

cm_save_plot <- function(name, draw, analysis, config, rows = 1L, ...) {
  width <- max(6, 3.4 * length(config$pve_values))
  height <- 3.0 * rows + 1.35
  pdf(file.path(config$output_dir, paste0(name, ".pdf")), width = width, height = height, useDingbats = FALSE)
  tryCatch(draw(analysis, config, ...), finally = dev.off())
  if (config$write_png) {
    png(file.path(config$output_dir, paste0(name, ".png")), width = width, height = height, units = "in", res = 160)
    tryCatch(draw(analysis, config, ...), finally = dev.off())
  }
}

run_coding_plots <- function(config = coding_plot_settings) {
  stopifnot(config$snp_cutoff > 0, config$snp_cutoff <= 1,
            config$snp_cutoff %in% config$thresholds, config$thresholds[1] == 0,
            all(diff(config$thresholds) > 0), all(config$thresholds <= 1),
            config$n_bins >= 2L, config$bootstrap_reps >= 2L)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  analysis <- cm_read_simulations(config)
  message("Summarizing mixture recovery and seed-block uncertainty...")
  analysis$mixture <- cm_mixture_summary(analysis$replicates, config$bootstrap_reps, config$bootstrap_seed)
  analysis$confusion <- cm_confusion_summary(analysis$replicates)
  analysis$classification <- cm_classification_rates(analysis$confusion)
  analysis$calibration <- cm_calibration_rates(analysis$calibration)
  analysis$discovery <- cm_discovery_rates(analysis$discovery, analysis$replicates)
  baseline <- analysis$replicates[analysis$replicates$scenario == "Additive only", ]
  analysis$false_calls <- if (nrow(baseline)) aggregate(
    baseline[paste0("any_false_", cm_short)], baseline[c("pve", "K")], mean) else data.frame()
  tables <- c(mixture = "mixture_recovery.csv", confusion = "coding_confusion.csv",
              classification = "coding_classification.csv", calibration = "coding_pip_calibration.csv",
              discovery = "coding_discovery_accuracy.csv", false_calls = "additive_only_false_call_rates.csv")
  for (name in names(tables)) write.csv(analysis[[name]], file.path(config$output_dir, tables[[name]]), row.names = FALSE)
  configuration_counts <- aggregate(list(n_replicates = analysis$replicates$seed),
    analysis$replicates[c("scenario", "configuration", "design", "all_additive", "pve", "K")], length)
  write.csv(configuration_counts, file.path(config$output_dir, "configuration_counts.csv"), row.names = FALSE)
  saveRDS(analysis$replicates, file.path(config$output_dir, "replicate_coding_metrics.rds"))
  saveRDS(config, file.path(config$output_dir, "plot_settings.rds"))
  cm_render_plots(analysis, config)
  message("Saved coding figures and tables in: ", config$output_dir)
  message("Included ", sum(analysis$audit$included), "; duplicates ", sum(analysis$audit$duplicates),
          "; simulation errors ", sum(analysis$audit$errors), "; invalid/ambiguous maps or results ", sum(analysis$audit$invalid),
          "; nonconverged mixed fits ", sum(analysis$audit$nonconverged))
  invisible(analysis)
}

cm_render_plots <- function(analysis, config) {
  cm_save_plot("additive_only_pip_allocation", cm_draw_baseline, analysis, config)
  cm_save_plot("mixture_recovery", cm_draw_recovery, analysis, config, rows = 3L)
  cm_save_plot("coding_confusion", cm_draw_confusion, analysis, config)
  cm_save_plot("coding_classification", cm_draw_classification, analysis, config, rows = 2L)
  cm_save_plot("coding_pip_calibration", cm_draw_calibration, analysis, config)
  cm_save_plot("coding_pip_calibration_additive_only", cm_draw_calibration, analysis, config, additive_only = TRUE)
  cm_save_plot("coding_discovery_accuracy", cm_draw_discoveries, analysis, config, rows = 2L)
  cm_save_plot("additive_only_false_calls", cm_draw_false_calls, analysis, config)
}

# Set this option to FALSE before sourcing to load functions/settings only.
if (isTRUE(getOption("susie_mix.coding_plots.run", TRUE))) {
  coding_plot_results <- run_coding_plots()
}
