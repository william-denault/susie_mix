# Exact SNP-and-coding discovery by the largest PIP in each simulated locus.
# Run with Rscript, or source this file in R/RStudio. Base R only; no refits.

project_dir <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (.Platform$OS.type == "windows")
  "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix" else
  "/project2/mstephens/wdenault/susie_mix")
source(file.path(project_dir, "script/sim/coding_mixture_utils.R"))

lead_pip_settings <- list(
  chunk_dir = file.path(project_dir, "simulation results/chunks"),
  output_dir = file.path(project_dir, "simulation results/lead_pip_figures"),
  file_pattern = "\\.RData$",
  n = 500, fit_L = 10, pve_values = c(.1, .2, .3, .4),
  n_bins = 10L,
  tie_tolerance = 1e-12,       # Absolute PIP difference from the maximum.
  exclude_nonconverged = FALSE,
  max_reps_per_file = Inf,    # e.g. 5 for a preview; entire checkpoints are loaded.
  write_png = TRUE
)
lp_outcomes <- c("exact", "wrong_coding", "wrong_snp", "ambiguous", "no_support")
lp_outcome_labels <- c(exact = "Exact SNP + coding", wrong_coding = "Wrong coding at causal SNP",
  wrong_snp = "Noncausal SNP", ambiguous = "Tied leads", no_support = "All PIPs zero")
lp_outcome_colors <- c(exact = "#238B66", wrong_coding = "#D77916", wrong_snp = "#BB3E70",
  ambiguous = "#7F7F7F", no_support = "#D5D5D5")
lp_coding_colors <- c(additive = "#2474B5", recessive = "#BB3E70", dominant = "#D77916")
lp_coding_labels <- c(additive = "Additive", recessive = "Recessive", dominant = "Dominant")

lp_assess_replicate <- function(x, tie_tolerance = 1e-12) {
  prepared <- cm_prepare_replicate(x)
  pip <- prepared$pip
  peak <- max(pip)
  top <- if (peak > 0) which(peak - pip <= tie_tolerance) else integer()
  lead <- if (length(top) == 1L) top else NA_integer_
  # Select globally before consulting truth. Do not restrict to causal SNPs,
  # credible sets, or a PIP threshold; never break ties by predictor order.
  outcome <- if (!length(top)) "no_support" else if (length(top) > 1L) "ambiguous" else
    if (lead %in% x$true_pos_mix) "exact" else
      if (x$mix_to_add[lead] %in% x$true_pos) "wrong_coding" else "wrong_snp"
  runner_up <- if (length(pip) > 1L) max(pip[-which.max(pip)]) else 0
  n_true <- tabulate(match(x$causal_coding, cm_classes), nbins = 3L)
  metrics <- c(setNames(n_true, paste0("true_", cm_short)),
    lead_pip = peak, lead_index = lead, lead_snp_index = x$mix_to_add[lead],
    lead_coding_id = match(prepared$labels$coding[lead], cm_classes),
    outcome_id = match(outcome, lp_outcomes), n_top = length(top), pip_gap = peak - runner_up,
    n_top_exact = sum(top %in% x$true_pos_mix))
  list(metrics = metrics, map_origin = prepared$labels$origin, n_clamped = prepared$n_clamped)
}

lp_scopes <- function(replicates) {
  all <- replicates
  all$scope <- "All simulations"
  additive <- replicates[replicates$scenario == "Additive only", , drop = FALSE]
  additive$scope <- rep("Additive only", nrow(additive))
  rbind(all, additive)
}

lp_outcome_summary <- function(d, keys) {
  values <- data.frame(n_runs = rep(1L, nrow(d)), sum_lead_pip = d$lead_pip)
  for (outcome in lp_outcomes) values[[paste0("n_", outcome)]] <- as.integer(d$outcome == outcome)
  for (cl in cm_classes) values[[paste0("n_lead_", cl)]] <- as.integer(!is.na(d$called_coding) & d$called_coding == cl)
  out <- aggregate(values, d[keys], sum)
  out$n_unique <- out$n_exact + out$n_wrong_coding + out$n_wrong_snp
  out$mean_lead_pip <- cm_ratio(out$sum_lead_pip, out$n_runs)
  for (outcome in lp_outcomes) out[[paste0("fraction_", outcome)]] <- cm_ratio(out[[paste0("n_", outcome)]], out$n_runs)
  out$exact_accuracy_unique <- cm_ratio(out$n_exact, out$n_unique)
  for (cl in cm_classes) out[[paste0("share_unique_", cl)]] <- cm_ratio(out[[paste0("n_lead_", cl)]], out$n_unique)
  out
}

lp_summarize <- function(replicates, config) {
  d <- lp_scopes(replicates)
  d$bin <- pmin(config$n_bins, floor(d$lead_pip * config$n_bins) + 1L)
  bins <- lp_outcome_summary(d, c("scope", "pve", "K", "bin"))
  bins$bin_lower <- (bins$bin - 1L) / config$n_bins
  bins$bin_upper <- bins$bin / config$n_bins
  unique <- d[!is.na(d$called_coding), , drop = FALSE]
  if (nrow(unique)) {
    accuracy <- aggregate(data.frame(n_leads = rep(1L, nrow(unique)), sum_lead_pip = unique$lead_pip,
      n_exact = as.integer(unique$outcome == "exact")), unique[c("scope", "pve", "K", "called_coding", "bin")], sum)
    accuracy$mean_lead_pip <- accuracy$sum_lead_pip / accuracy$n_leads
    accuracy$exact_accuracy <- accuracy$n_exact / accuracy$n_leads
  } else {
    accuracy <- data.frame(scope = character(), pve = numeric(), K = integer(), called_coding = character(),
      bin = integer(), n_leads = integer(), sum_lead_pip = numeric(), n_exact = integer(),
      mean_lead_pip = numeric(), exact_accuracy = numeric())
  }
  summary <- lp_outcome_summary(d, c("scope", "pve", "K"))
  conditions <- lp_outcome_summary(replicates, c("scenario", "configuration", "design", "all_additive", "pve", "K"))
  list(bins = bins, accuracy = accuracy, summary = summary, conditions = conditions,
    additive = summary[summary$scope == "Additive only", , drop = FALSE])
}

lp_layout <- function(nr, nc) {
  par(mfrow = c(nr, nc), mar = c(3.1, 3.7, 2.2, .8), oma = c(3.3, .2, 4.4, .2),
    mgp = c(1.9, .6, 0), tcl = -.25, family = "sans", cex = .85)
}

lp_panel <- function(title, ylab, xlab = "Lead PIP bin", xlim = c(0, 1), counts = FALSE) {
  plot(NA, xlim = xlim, ylim = if (counts) c(0, 1.13) else c(0, 1), axes = FALSE,
    xlab = xlab, ylab = "", main = title, xaxs = "i", yaxs = "i", cex.main = .95)
  axis(1, at = if (identical(xlim, c(0, 1))) c(0, .5, 1) else pretty(xlim), las = 1)
  axis(2, at = c(0, .25, .5, .75, 1), las = 1)
  box(bty = "l")
  abline(h = c(0, .25, .5, .75, 1), col = "#EEEEEE")
  mtext(ylab, side = 2, line = 2.7)
}

lp_heading <- function(title, subtitle) {
  width <- par("din")[1]
  mtext(paste(strwrap(title, width = floor(width * 9)), collapse = "\n"),
    side = 3, outer = TRUE, line = 2.6, font = 2, cex = 1.2)
  mtext(paste(strwrap(subtitle, width = floor(width * 14)), collapse = "\n"),
    side = 3, outer = TRUE, line = 1, cex = .8)
}

lp_legend <- function(labels, colors, fill = FALSE) {
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  ncol <- if (length(labels) > 3 && par("din")[1] < 10) 3L else length(labels)
  legend("bottom", legend = unname(labels), col = if (fill) NA else unname(colors),
    fill = if (fill) unname(colors) else NULL, border = NA,
    pch = if (fill) NULL else 16, ncol = ncol, bty = "n", cex = .85, inset = .006)
}

lp_draw_outcomes <- function(analysis, config, additive_only = FALSE) {
  scope <- if (additive_only) "Additive only" else "All simulations"
  d <- analysis$bins[analysis$bins$scope == scope, , drop = FALSE]
  ks <- sort(unique(analysis$replicates$K))
  lp_layout(length(ks), length(config$pve_values))
  for (K in ks) for (pve in config$pve_values) {
    lp_panel(paste0("K = ", K, "  |  PVE = ", 100 * pve, "%"), "Fraction of simulations", counts = TRUE)
    x <- d[d$K == K & d$pve == pve, , drop = FALSE]
    if (!nrow(x)) { text(.5, .5, "No saved simulations"); next }
    for (i in seq_len(nrow(x))) {
      heights <- as.numeric(x[i, paste0("fraction_", lp_outcomes)])
      edges <- c(0, cumsum(heights))
      padding <- .05 / config$n_bins
      rect(x$bin_lower[i] + padding, head(edges, -1), x$bin_upper[i] - padding, tail(edges, -1),
        col = lp_outcome_colors, border = "white", lwd = .3)
      text((x$bin_lower[i] + x$bin_upper[i]) / 2, 1.065, x$n_runs[i], cex = .68, col = "#555555")
    }
  }
  lp_heading(if (additive_only) "Lead-PIP discoveries under purely additive truth" else "Does the largest PIP identify the exact SNP and coding?",
    "One global lead per locus; ties and zero support remain in the denominator. Counts above bars; empty bins are blank.")
  lp_legend(lp_outcome_labels, lp_outcome_colors, fill = TRUE)
}

lp_draw_accuracy <- function(analysis, config) {
  d <- analysis$accuracy[analysis$accuracy$scope == "All simulations", , drop = FALSE]
  ks <- sort(unique(analysis$replicates$K))
  lp_layout(length(ks), length(config$pve_values))
  for (K in ks) for (pve in config$pve_values) {
    lp_panel(paste0("K = ", K, "  |  PVE = ", 100 * pve, "%"), "Exact-match fraction", "Mean lead PIP in bin")
    abline(0, 1, lty = 2, col = "#777777")
    panel <- d[d$K == K & d$pve == pve, , drop = FALSE]
    if (!nrow(panel)) { text(.5, .5, "No unambiguous leads"); next }
    for (cl in cm_classes) {
      x <- panel[panel$called_coding == cl, , drop = FALSE]
      x <- x[order(x$bin), , drop = FALSE]
      # Empty bins break curves instead of implying observations in gaps.
      xx <- yy <- rep(NA_real_, config$n_bins)
      xx[x$bin] <- x$mean_lead_pip; yy[x$bin] <- x$exact_accuracy
      lines(xx, yy, col = lp_coding_colors[cl])
      points(x$mean_lead_pip, x$exact_accuracy, col = lp_coding_colors[cl], pch = 16,
        cex = .55 + .5 * pmin(1, log10(pmax(1, x$n_leads)) / 3), xpd = NA)
    }
  }
  lp_heading("Exact discovery accuracy by the lead's reported coding",
    "Positive, unambiguous leads only; colors indicate the called coding. Larger points have more leads; dashed line: equality.")
  lp_legend(lp_coding_labels, lp_coding_colors)
}

lp_draw_additive <- function(analysis, config) {
  d <- analysis$additive
  ks <- sort(unique(analysis$replicates$K))
  lp_layout(1, length(config$pve_values))
  par(mar = c(3.1, 3.7, 3.1, .8))
  for (pve in config$pve_values) {
    lp_panel("", "Share of unambiguous leads",
      "Number of true causal SNPs", range(ks) + c(-.35, .35))
    mtext(paste0("PVE = ", 100 * pve, "%"), side = 3, line = 1.55, font = 2, cex = .95)
    x <- d[d$pve == pve, , drop = FALSE]; x <- x[order(x$K), , drop = FALSE]
    if (!nrow(x)) { text(mean(ks), .5, "No additive-only simulations"); next }
    for (cl in cm_classes) lines(x$K, x[[paste0("share_unique_", cl)]], type = "b", pch = 16,
      col = lp_coding_colors[cl], xpd = NA)
    # Panel-level exclusions are visible; the table also reports them for each K.
    mtext(paste0("Unique leads: ", sum(x$n_unique), " / ", sum(x$n_runs), " runs"), side = 3, line = .3, cex = .7)
  }
  lp_heading("Which coding wins when every true effect is additive?",
    "Count of leads using each coding / all positive, unambiguous leads; every SNP competes. This is a ranking frequency, not a PIP-mass share.")
  lp_legend(lp_coding_labels, lp_coding_colors)
}

lp_save_plot <- function(name, draw, analysis, config, rows = 1L, ...) {
  width <- max(7, 3.4 * length(config$pve_values))
  height <- 2.45 * rows + 1.55
  pdf(file.path(config$output_dir, paste0(name, ".pdf")), width = width, height = height, useDingbats = FALSE)
  tryCatch(draw(analysis, config, ...), finally = dev.off())
  if (config$write_png) {
    png(file.path(config$output_dir, paste0(name, ".png")), width = width, height = height, units = "in", res = 160)
    tryCatch(draw(analysis, config, ...), finally = dev.off())
  }
}

lp_render_plots <- function(analysis, config) {
  nr <- length(unique(analysis$replicates$K))
  lp_save_plot("lead_pip_outcomes", lp_draw_outcomes, analysis, config, rows = nr)
  lp_save_plot("lead_pip_outcomes_additive_only", lp_draw_outcomes, analysis, config, rows = nr, additive_only = TRUE)
  lp_save_plot("lead_pip_exact_accuracy", lp_draw_accuracy, analysis, config, rows = nr)
  lp_save_plot("additive_only_lead_coding", lp_draw_additive, analysis, config)
}

run_lead_pip_plots <- function(config = lead_pip_settings) {
  stopifnot(config$n_bins >= 2L, config$n_bins == as.integer(config$n_bins),
    is.finite(config$tie_tolerance), config$tie_tolerance >= 0, config$tie_tolerance < 1,
    length(config$pve_values) > 0L, config$max_reps_per_file > 0)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  stream <- cm_read_simulations(config,
    assess_replicate = function(x, ...) lp_assess_replicate(x, config$tie_tolerance))
  replicates <- stream$replicates
  replicates$called_coding <- cm_classes[replicates$lead_coding_id]
  replicates$outcome <- lp_outcomes[replicates$outcome_id]
  analysis <- c(list(replicates = replicates, audit = stream$audit, issues = stream$issues), lp_summarize(replicates, config))
  tables <- c(bins = "lead_pip_outcomes_by_bin.csv", accuracy = "lead_pip_accuracy_by_coding.csv",
    summary = "lead_pip_summary.csv", conditions = "lead_pip_outcomes_by_condition.csv", additive = "additive_only_lead_coding.csv")
  for (name in names(tables)) write.csv(analysis[[name]], file.path(config$output_dir, tables[[name]]), row.names = FALSE)
  saveRDS(replicates, file.path(config$output_dir, "lead_pip_replicates.rds"))
  saveRDS(config, file.path(config$output_dir, "plot_settings.rds"))
  saveRDS(analysis, file.path(config$output_dir, "lead_pip_analysis.rds"))
  lp_render_plots(analysis, config)
  message("Saved lead-PIP figures and tables in: ", config$output_dir)
  message("Included ", nrow(replicates), "; exact ", sum(replicates$outcome == "exact"),
    "; wrong coding ", sum(replicates$outcome == "wrong_coding"), "; noncausal SNP ", sum(replicates$outcome == "wrong_snp"),
    "; tied ", sum(replicates$outcome == "ambiguous"), "; zero support ", sum(replicates$outcome == "no_support"))
  message("Duplicates ", sum(stream$audit$duplicates), "; simulation errors ", sum(stream$audit$errors),
    "; invalid ", sum(stream$audit$invalid), "; nonconverged ", sum(stream$audit$nonconverged))
  invisible(analysis)
}

# Set this option to FALSE before sourcing to load functions/settings only.
if (isTRUE(getOption("susie_mix.lead_pip_plots.run", TRUE))) {
  lead_pip_results <- run_lead_pip_plots()
}
