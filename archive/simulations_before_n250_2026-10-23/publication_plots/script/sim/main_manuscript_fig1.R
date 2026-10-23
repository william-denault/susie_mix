# Main manuscript simulation figure: six columns, each with three rows.
# Run from the project root:
#   source("script/sim/main_manuscript_fig1.R")
#
# Required packages: ggplot2 and patchwork.
# Reads existing summaries only; no simulations or model fits are run.
#
# Rows: additive only, recessive only, partially recessive only.
# Columns:
#   A-C: SNP PIP calibration at PVE = 10%, one method per column.
#   D:   Power versus empirical FDR at PVE = 30%, all three methods.
#   E:   Mean credible-set purity at PVE = 30%, all three methods.
#   F:   Mean credible-set size at PVE = 30%, all three methods.
#
# Calibration and power/FDR pool available true causal counts K = 1 through 4.
# Calibration bars are the saved +/-2 seed-block empirical SE intervals.
# Power/FDR pools TP and FP counts at each PIP threshold, preserving threshold
# order; it does not average per-K rates or construct a monotone envelope.
# In panel D's recessive-only row, the SuSiE segment below 1% power is omitted
# from the display. Original and displayed curves are saved separately.
# Open/filled circles mark the saved 0.95/0.80 PIP thresholds. Saved counts use
# PIP >= cutoff; no PIPs equal these cutoffs in the selected PVE = 30% runs
# (checked against the checkpoints), so these also give the strict > counts.
# Purity = sum of CS minimum absolute correlations / number of reported CSs.
# Purity shows dodged points without intervals or connecting lines, with a
# shared lower y limit just below the lowest displayed mean.
# CS size counts unique biological SNPs per reported CS and retains the saved
# Gaussian intervals. Its shared y range includes all displayed intervals.
#
# To customize without running the defaults:
#   options(susie.manuscript.autorun = FALSE)
#   source("script/sim/main_manuscript_fig1.R")
#   result <- make_manuscript_simulation_figure(width = 14.4, height = 7.2)
#   result$figure  # patchwork object containing editable ggplots

make_manuscript_simulation_figure <- function(
    project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", getwd()),
    input_dir = file.path(project_dir, "input_data"),
    output_dir = file.path(project_dir, "figures/figure1"),
    calibration_pve = 0.10,
    performance_pve = 0.30,
    fdr_max = 0.25,
    recessive_susie_min_power = 0.01,
    purity_padding = 0.01,
    purity_dodge_width = 0.72,
    max_causal = 4L,
    width = 14.4,
    height = 7.2,
    base_size = 10,
    png_dpi = 300) {
  packages <- c("ggplot2", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the required R packages first: install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), "))")
  }
  stopifnot(calibration_pve > 0, calibration_pve < 1,
            performance_pve > 0, performance_pve < 1,
            fdr_max > 0, fdr_max <= 1,
            recessive_susie_min_power >= 0, recessive_susie_min_power < 1,
            purity_padding > 0,
            purity_dodge_width >= 0, purity_dodge_width < 1)
  # Reuse the established curve-clipping and axis helpers. These only define
  # functions; sourcing this file does not read checkpoints or run simulations.
  helper_path <- file.path(project_dir, "script/sim/simulation_plot_helpers.R")
  source(helper_path, local = TRUE)

  publication_helper <- file.path(project_dir, "script/sim/publication_plot_helpers.R")
  source(publication_helper, local = TRUE)
  methods <- publication_methods
  method_labels <- publication_labels
  scenarios <- c("Additive only", "Recessive only", "Partial recessive only")
  row_labels <- setNames(c("Additive\nonly", "Recessive\nonly", "Partially\nrecessive"), scenarios)
  method_colors <- publication_colors
  method_shapes <- publication_shapes
  method_lines <- setNames(c("solid", "dashed", "dotdash"), methods)
  method_point_sizes <- setNames(c(1.8, 1.8, 2.6), methods)
  saved <- read_publication_simulation_data(input_dir, max_causal, methods)
  paths <- saved$paths
  calibration <- saved$pooled_calibration
  curves <- saved$pooled_curves
  metrics <- saved$metrics

  select_data <- function(d, pve_value) {
    d <- d[d$scenario %in% scenarios & d$method %in% methods &
             abs(d$pve - pve_value) < 1e-8, , drop = FALSE]
    d$scenario <- factor(d$scenario, levels = scenarios)
    d$method <- factor(d$method, levels = methods)
    if (any(table(d$scenario, d$method) == 0)) {
      stop("A requested scenario/method is missing at PVE = ", pve_value)
    }
    d
  }
  calibration <- select_data(calibration, calibration_pve)
  curves <- select_data(curves, performance_pve)
  purity <- select_data(metrics, performance_pve)
  calibration_metrics <- select_data(metrics, calibration_pve)

  # Check that the pooled inputs describe the same available K values and
  # replicate counts as the metric summary; never silently label absent K = 5.
  k_values <- sort(unique(purity$K))
  for (s in scenarios) for (m in methods) {
    pm <- purity[purity$scenario == s & purity$method == m, ]
    cm <- calibration_metrics[calibration_metrics$scenario == s &
                                calibration_metrics$method == m, ]
    cb <- calibration[calibration$scenario == s & calibration$method == m, ]
    cr <- curves[curves$scenario == s & curves$method == m, ]
    baseline <- cr[cr$threshold == 0, ]
    if (!identical(sort(pm$K), k_values) ||
        !identical(sort(cm$K), k_values)) {
      stop("Available causal counts differ across the selected panels.")
    }
    if (nrow(cb) != 10L || anyDuplicated(cb$bin) || any(cb$K != 0) ||
        any(cb$n_replicates != sum(cm$n_replicates)) ||
        sum(cb$n_causal) != sum(cm$n_causal) ||
        nrow(baseline) != 1L || baseline$tp != sum(pm$n_causal) ||
        anyDuplicated(cr$threshold)) {
      stop("Inconsistent saved summaries for ", s, " / ", m,
           ". Refresh the simulation summaries before plotting.")
    }
  }
  if (any(!is.finite(purity$purity)) ||
      any(purity$purity < 0 | purity$purity > 1)) {
    stop("Missing or invalid purity estimates in the selected panels.")
  }
  purity_limits <- c(max(0, min(purity$purity) - purity_padding), 1)
  if (purity_limits[1] == min(purity$purity)) {
    purity_limits[1] <- min(purity$purity) - purity_padding
  }
  purity <- purity[order(purity$scenario, purity$method, purity$K), ]
  cs_size <- purity
  if (any(!is.finite(cs_size$cs_size)) || any(cs_size$cs_size < 1)) {
    stop("Missing or invalid credible-set sizes in the selected panels.")
  }
  cs_size_limits <- c(0, 1.04 * max(c(cs_size$cs_size, cs_size$cs_size_hi), na.rm = TRUE))
  calibration <- calibration[order(calibration$scenario, calibration$method,
                                   calibration$bin), ]
  # geom_path connects rows in their saved threshold order. geom_line would
  # sort by empirical FDR and incorrectly alter these potentially looping paths.
  curves <- curves[order(curves$scenario, curves$method, -curves$threshold), ]
  # Display-only removal of the annotated near-zero tail in the middle panel.
  # Do not connect the remaining curve to the origin or reorder it by FDR.
  # Set recessive_susie_min_power = 0 to show the full saved curve again.
  omit_curve_point <- curves$scenario == "Recessive only" &
    curves$method == "SuSiE" & curves$tpr < recessive_susie_min_power
  displayed_curves <- curves[!omit_curve_point, , drop = FALSE]
  # Include intersections of segments with FDR = fdr_max, not just the saved
  # vertices inside the window. The upper axis endpoint has no extra padding.
  power_limits <- c(0, visible_curve_max(displayed_curves, c(0, fdr_max)))
  if (power_limits[2] <= 0) stop("No positive power within the displayed FDR range.")
  pip_cutoffs <- c(0.95, 0.80)
  threshold_markers <- do.call(rbind, lapply(pip_cutoffs, function(cutoff) {
    d <- curves[abs(curves$threshold - cutoff) < 1e-12, ]
    if (nrow(d) != length(scenarios) * length(methods) ||
        any(table(d$scenario, d$method) != 1L)) {
      stop("Expected one saved curve point per scenario/method at PIP = ", cutoff)
    }
    d$pip_cutoff <- sprintf(">%.2f", cutoff)
    d
  }))
  threshold_markers$pip_cutoff <- factor(threshold_markers$pip_cutoff,
                                       levels = c(">0.95", ">0.80"))

  theme_manuscript <- ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.25),
      panel.border = ggplot2::element_rect(color = "#777777", linewidth = 0.35),
      panel.spacing.y = grid::unit(0.17, "in"),
      axis.ticks = ggplot2::element_line(color = "#777777", linewidth = 0.3),
      axis.text = ggplot2::element_text(color = "#222222", size = base_size - 1),
      axis.title = ggplot2::element_text(size = base_size),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(
        angle = 0, face = "bold", size = base_size, lineheight = 1.05,
        margin = ggplot2::margin(r = 7)),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.width = grid::unit(0.32, "in"))

  facet_rows <- function() {
    ggplot2::facet_grid(
      rows = ggplot2::vars(scenario), switch = "y", drop = FALSE,
      labeller = ggplot2::labeller(scenario = row_labels))
  }
  method_scales <- function(shape = TRUE, linetype = TRUE) {
    scales <- list(ggplot2::scale_color_manual(
      name = NULL, values = method_colors, breaks = methods, labels = method_labels, drop = FALSE))
    if (shape) scales <- c(scales, list(ggplot2::scale_shape_manual(
      name = NULL, values = method_shapes, breaks = methods, labels = method_labels, drop = FALSE)))
    if (linetype) scales <- c(scales, list(ggplot2::scale_linetype_manual(
      name = NULL, values = method_lines, breaks = methods, labels = method_labels, drop = FALSE)))
    scales
  }
  hide_row_labels <- ggplot2::theme(strip.text.y.left = ggplot2::element_blank())

  calibration_plots <- lapply(seq_along(methods), function(i) {
    d <- calibration[calibration$method == methods[i], ]
    p <- ggplot2::ggplot(d, ggplot2::aes(mean_pip, frequency, color = method,
                                       shape = method)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                          color = "#888888", linewidth = 0.4) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                            width = 0.025, linewidth = 0.3, alpha = 0.65,
                            na.rm = TRUE, show.legend = FALSE) +
      ggplot2::geom_point(size = unname(method_point_sizes[methods[i]]),
                          na.rm = TRUE, show.legend = FALSE) +
      facet_rows() + method_scales(linetype = FALSE) +
      ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1),
                                 expand = ggplot2::expansion(mult = 0.025)) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25),
                                 expand = ggplot2::expansion(mult = 0.025)) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste0(LETTERS[i], "  ", method_labels[methods[i]]),
        x = "Mean PIP", y = NULL) +
      theme_manuscript
    if (i > 1) p <- p + hide_row_labels
    p
  })

  power_plot <- ggplot2::ggplot(
    displayed_curves,
    ggplot2::aes(fdr, tpr, color = method, linetype = method, group = method)) +
    ggplot2::geom_path(linewidth = 0.55) +
    ggplot2::geom_point(
      data = threshold_markers,
      ggplot2::aes(fdr, tpr, color = method, shape = pip_cutoff),
      inherit.aes = FALSE, size = 2.6, stroke = 0.7, fill = "white",
      show.legend = c(color = FALSE, shape = TRUE)) +
    facet_rows() + method_scales(shape = FALSE) +
    ggplot2::scale_shape_manual(name = "PIP", values = c(">0.95" = 21, ">0.80" = 16)) +
    ggplot2::scale_x_continuous(
      breaks = pretty(c(0, fdr_max), n = 3),
      expand = ggplot2::expansion(mult = 0.025)) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, power_limits[2], by = 0.1),
      labels = function(x) formatC(x, format = "f", digits = 1)) +
    # Coordinate zoom keeps path segments crossing the right boundary intact.
    ggplot2::coord_cartesian(xlim = c(0, fdr_max), ylim = power_limits, expand = FALSE) +
    ggplot2::labs(title = "D  Power vs FDR",
                  x = "Empirical FDR", y = "Power") +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(order = 2, override.aes = list(color = "#333333"))) +
    theme_manuscript + hide_row_labels +
    ggplot2::theme(legend.title = ggplot2::element_text(size = base_size))

  purity_position <- ggplot2::position_dodge(width = purity_dodge_width)
  purity_plot <- ggplot2::ggplot(
    purity, ggplot2::aes(K, purity, color = method, shape = method,
                        group = method)) +
    ggplot2::geom_point(ggplot2::aes(size = method),
                        position = purity_position, show.legend = FALSE) +
    ggplot2::scale_size_manual(values = method_point_sizes, guide = "none") +
    facet_rows() + method_scales(linetype = FALSE) +
    ggplot2::scale_x_continuous(breaks = k_values,
                               expand = ggplot2::expansion(add = 0.2)) +
    ggplot2::scale_y_continuous(
      breaks = pretty(purity_limits, n = 4),
      labels = function(x) formatC(x, format = "f", digits = 2)) +
    ggplot2::coord_cartesian(ylim = purity_limits, expand = c(bottom = FALSE, top = FALSE)) +
    ggplot2::labs(title = "E  Credible-set purity",
                  x = "Number of causal SNPs", y = "purity") +
    theme_manuscript + hide_row_labels

  cs_size_plot <- ggplot2::ggplot(
    cs_size, ggplot2::aes(K, cs_size, color = method, shape = method,
                         group = method)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = cs_size_lo, ymax = cs_size_hi),
                            position = purity_position, width = 0.12,
                            linewidth = 0.3, alpha = 0.65,
                            na.rm = TRUE, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(size = method),
                        position = purity_position, show.legend = FALSE) +
    ggplot2::scale_size_manual(values = method_point_sizes, guide = "none") +
    facet_rows() + method_scales(linetype = FALSE) +
    ggplot2::scale_x_continuous(breaks = k_values,
                               expand = ggplot2::expansion(add = 0.2)) +
    ggplot2::scale_y_continuous(breaks = pretty(cs_size_limits, n = 4)) +
    ggplot2::coord_cartesian(ylim = cs_size_limits,
                             expand = c(bottom = FALSE, top = FALSE)) +
    ggplot2::labs(title = "F  CS size", x = "Number of causal SNPs",
                  y = "Mean SNPs per CS") +
    theme_manuscript + hide_row_labels

  figure <- patchwork::wrap_plots(
    c(calibration_plots, list(power_plot, purity_plot, cs_size_plot)),
    nrow = 1, widths = rep(1, 6), guides = "collect") &
    ggplot2::theme(legend.position = "bottom")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(output_dir, "main_simulations")
  ggplot2::ggsave(paste0(stem, ".pdf"), figure, width = width, height = height,
                  units = "in", device = grDevices::cairo_pdf, bg = "white")
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(paste0(stem, ".png"), figure, width = width, height = height,
                  units = "in", dpi = png_dpi, device = png_device, bg = "white")

  # Save exactly the plotted data and provenance for manuscript figure review.
  plot_data <- list(
    calibration = calibration, power_fdr = curves, threshold_markers = threshold_markers,
    power_fdr_display = displayed_curves,
    power_fdr_omitted = curves[omit_curve_point, , drop = FALSE],
    purity = purity, cs_size = cs_size,
    settings = list(methods = methods, labels = method_labels, colors = method_colors, max_causal = max_causal, scenarios = scenarios, causal_counts = k_values,
                    calibration_pve = calibration_pve, performance_pve = performance_pve,
                    fdr_max = fdr_max, power_limits = power_limits,
                    recessive_susie_min_power = recessive_susie_min_power,
                    purity_limits = purity_limits, purity_dodge_width = purity_dodge_width,
                    cs_size_limits = cs_size_limits),
    source_md5 = tools::md5sum(c(paths, helper_path, publication_helper)))
  saveRDS(plot_data, paste0(stem, "_data.rds"))
  message("Saved: ", normalizePath(output_dir, winslash = "/"))
  message("Pooled available causal counts: ", paste(k_values, collapse = ", "))
  message("Shared purity y range: ", paste(round(purity_limits, 4), collapse = " to "))
  message("Shared power y range: ", paste(round(power_limits, 4), collapse = " to "))
  message("Shared CS-size y range: ", paste(round(cs_size_limits, 4), collapse = " to "))
  invisible(list(figure = figure, plots = c(calibration_plots,
                 list(power_plot, purity_plot, cs_size_plot)), data = plot_data,
                 files = paste0(stem, c(".pdf", ".png", "_data.rds"))))
}

if (isTRUE(getOption("susie.manuscript.autorun", TRUE))) {
  manuscript_simulation_figure <- make_manuscript_simulation_figure()
}
