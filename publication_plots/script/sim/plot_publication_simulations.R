# Publication simulation panels, using existing summaries only.
# Run generate_publication_plots.R for these panels and Figure 1.
# The printed CI formula is omitted; coverage, power and CS-size bars remain.
# Purity has no error bars. Calibration retains its seed-block empirical bars.
make_publication_simulation_figures <- function(
    project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", getwd()),
    input_dir = file.path(dirname(project_dir), "simulation results/n250_r300_v1/figures"),
    output_dir = file.path(project_dir, "figures"),
    max_causal = 5L, write_png = TRUE) {
  source(file.path(project_dir, "script/sim/simulation_design.R"), local = TRUE)
  source(file.path(project_dir, "script/sim/simulation_plot_helpers.R"), local = TRUE)
  source(file.path(project_dir, "script/sim/publication_plot_helpers.R"), local = TRUE)
  stopifnot(max_causal %in% 3:5)
  method_names <- publication_methods
  method_labels <- publication_labels
  method_colors <- publication_colors
  method_symbols <- publication_shapes
  method_lines <- publication_lines
  pve_values <- publication_pve_values
  n_value <- publication_n
  target_coverage <- .95
  roc_max_fpr <- fdr_max <- .25
  write_roc_all_L <- TRUE
  write_roc_pages <- FALSE
  scenario_rows <- lapply(1:3, function(k) {
    vapply(combn(5L, k, simplify = FALSE), function(ids) {
      counts <- integer(5L)
      counts[ids] <- 1L
      sim_scenario_name(counts, display = TRUE)
    }, "")
  })
  scenario_groups <- list(pure = scenario_rows[[1]],
    pairs_1 = scenario_rows[[2]][1:5], pairs_2 = scenario_rows[[2]][6:10],
    triples_1 = scenario_rows[[3]][1:5], triples_2 = scenario_rows[[3]][6:10])
  message("Reading n=", n_value, " summaries; selecting the three publication methods and K = 1:", max_causal, ".")
  data <- read_publication_simulation_data(input_dir, max_causal, method_names)
  summary_table <- data$metrics
  roc_table <- data$curves
  pooled_roc_table <- data$pooled_curves
  calibration_table <- data$calibration
  pooled_calibration_table <- data$pooled_calibration
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Keep exact saved identifiers in data exports; use publication labels in plots.
  write.csv(summary_table, file.path(output_dir, "metric_summary.csv"), row.names = FALSE)
  write.csv(calibration_table, file.path(output_dir, "pip_calibration.csv"), row.names = FALSE)
  write.csv(pooled_calibration_table, file.path(output_dir, "pip_calibration_all_L.csv"), row.names = FALSE)
  saveRDS(roc_table, file.path(output_dir, "roc_counts.rds"))
  saveRDS(pooled_roc_table, file.path(output_dir, "roc_counts_all_L.rds"))
  saveRDS(list(methods = method_names, labels = method_labels, colors = method_colors,
    shapes = method_symbols, causal_counts = seq_len(max_causal), n = n_value,
    pve_values = pve_values, input_dir = normalizePath(input_dir, winslash = "/"),
    point_offsets = c(-.24, 0, .24),
    power_limits = c(0, 1), power_zoom_limits = c(.75, 1), purity_error_bars = FALSE,
    source_md5 = tools::md5sum(data$paths)), file.path(output_dir, "publication_settings.rds"))

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
  yticks <- if (metric == "power" && ylim[1] == .75) seq(.75, 1, .05) else if (metric %in% c("power", "roc", "pip_calibration")) seq(0, 1, .25) else panel_y_ticks(ylim)

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
        xlim <- c(min_K - 0.4, max_causal + 0.4)
        xticks <- seq.int(min_K, max_causal)
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
                                 pch = method_symbols[m], cex = .75,
                                 col = method_colors[m], xpd = NA)
        } else if (is_curve) {
          # The all-L table already pools the underlying TP/FP counts.
          # FDR need not increase monotonically: preserve threshold order,
          # rather than sorting FDR or reporting an optimized envelope.
          dm <- dm[order(dm$threshold, decreasing = TRUE), ]
          seg <- visible_curve_segments(if (is_fdr) dm$fdr else dm$fpr, dm$tpr, xlim)
          if (nrow(seg)) segments(seg[, "x0"], pmax(ylim[1], pmin(ylim[2], seg[, "y0"])),
                                 seg[, "x1"], pmax(ylim[1], pmin(ylim[2], seg[, "y1"])),
                                 col = method_colors[m], lty = method_lines[m], lwd = 1.5)
        } else {
          xpos <- dm$K + seq(-.24, .24, length.out = length(method_names))[m]
          lo <- dm[[paste0(metric, "_lo")]]
          hi <- dm[[paste0(metric, "_hi")]]
          if (metric != "purity" && !is.null(lo) && !is.null(hi)) {
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
                                 pch = method_symbols[m], cex = if (method_symbols[m] == 18) 1.2 else 0.95, col = method_colors[m], xpd = NA)
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
  if (is_pip && is.null(only_K)) title_text <- paste0(title_text, "  |  Causal counts up to ", max_causal, " pooled")
  if (length(methods) == 1L) title_text <- paste0(title_text, "  |  ", method_labels[methods])
  mtext(title_text, side = 3, outer = TRUE, line = 1.2, font = 2, cex = 1.1)
  mtext(if (is_calibration) "Mean PIP within bin" else if (is_fdr) "Empirical FDR" else if (is_roc) "False positive rate" else "Number of causal SNPs",
        side = 1, outer = TRUE, line = 0.5)
  ylab <- c(coverage = "Coverage", purity = "purity",
            power = "Power", cs_size = "Mean number of unique SNPs per CS",
            roc = "True positive rate (power)", power_fdr = "Power (true positive rate)",
            pip_calibration = "Fraction of SNPs in bin that are causal")
  mtext(ylab[metric], side = 2, outer = TRUE, line = 1.8)

  # Draw a common legend in the outer bottom margin.
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend(.5, .045, legend = unname(method_labels[methods]), col = method_colors[method_ids],
         pch = if (is_curve) NA else method_symbols[method_ids], lty = if (is_curve) method_lines[method_ids] else NA,
         lwd = 1.5, horiz = TRUE, xjust = .5, yjust = .5, bty = "n", cex = .95)
  if (is_calibration) {
    text(.5, .015, "Ten equal-width PIP bins; error bars: +/-2 empirical SE across seed blocks", cex = .85)

  }
}

save_figure <- function(metric, scenarios, name, only_K = NULL, methods = method_names, y_limits = NULL) {
  height <- 2.0 * length(scenarios) + 1.4
  pdf(file.path(output_dir, paste0(name, ".pdf")), width = 16, height = height,
      useDingbats = FALSE)
  tryCatch(draw_figure(metric, scenarios, only_K, methods = methods, y_limits = y_limits), finally = dev.off())
  if (write_png) {
    png(file.path(output_dir, paste0(name, ".png")), width = 16, height = height,
        units = "in", res = 180,
        type = if (capabilities("cairo")) "cairo-png" else getOption("bitmapType"))
    tryCatch(draw_figure(metric, scenarios, only_K, methods = methods, y_limits = y_limits), finally = dev.off())
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
    ks <- seq.int(min_K, max_causal)
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
for (group in names(scenario_groups))
  save_figure("power", scenario_groups[[group]], paste0("power_", group, "_zoom_075_1"),
              y_limits = c(.75, 1))
save_roc_figures()
save_roc_figures("power_fdr")

save_calibration_figures <- function() {
  # Calibration exports always pool all L and show one method per figure.
  for (group in names(scenario_groups)) {
    for (method in method_names) {
      method_suffix <- tolower(gsub("[^[:alnum:]]+", "_", method_labels[method]))
      save_figure("pip_calibration", scenario_groups[[group]],
        paste0("pip_calibration_", group, "_all_L_", method_suffix), methods = method)
    }
  }
}

save_calibration_figures()


  message("Publication simulation panels saved in: ", normalizePath(output_dir, winslash = "/"))
  invisible(list(data = data, output_dir = output_dir))
}
