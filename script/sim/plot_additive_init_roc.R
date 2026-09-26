# SNP-level ROC curves for the additive initialization experiment; base R only.
# In RCC RStudio: set SUSIE_MIX_PROJECT_DIR, then source this whole file.
# To load functions without plotting: options(susie.init.roc.autorun = FALSE).

additive_roc_methods <- c("SuSiE", "SuSiE-slide", "SuSiE-init-slide")
additive_roc_colors <- c("#2489FF", "#009E73", "#D81B60")

# Match the main simulation plots: exact TP/FP counts on a dense common PIP
# threshold grid, including both endpoints. No SNP-level vectors are pooled
# across checkpoints in memory. A discovery means PIP >= threshold.
additive_roc_thresholds <- sort(unique(c(
  0, 10^seq(-12, -2, length.out = 301), seq(.01, .99, by = .001),
  1 - 10^seq(-2, -12, length.out = 301), 1
)))

additive_roc_counts <- function(pip, truth, thresholds = additive_roc_thresholds) {
  if (!is.numeric(pip) || !length(pip) || any(!is.finite(pip)) ||
      any(pip < -1e-10 | pip > 1 + 1e-10)) stop("Invalid saved SNP PIPs.")
  if (!is.numeric(truth) || !length(truth) || any(!is.finite(truth)) ||
      any(truth != floor(truth) | truth < 1 | truth > length(pip)) || anyDuplicated(truth))
    stop("Invalid causal SNP indices for the saved PIPs.")
  if (!is.numeric(thresholds) || length(thresholds) < 2L || any(!is.finite(thresholds)) ||
      thresholds[1L] != 0 || tail(thresholds, 1L) != 1 || any(diff(thresholds) <= 0))
    stop("Thresholds must increase strictly from 0 to 1.")
  bins <- findInterval(pmin(1, pmax(0, pip)), thresholds)
  all_bins <- tabulate(bins, length(thresholds))
  causal_bins <- tabulate(bins[truth], length(thresholds))
  cbind(tp = c(rev(cumsum(rev(causal_bins))), 0),
        fp = c(rev(cumsum(rev(all_bins - causal_bins))), 0))
}

plot_additive_init_roc <- function(
    results_dir = file.path(Sys.getenv("SUSIE_MIX_PROJECT_DIR",
      "/project2/mstephens/wdenault/susie_mix"), "simulation results/additive_slide_init_v3"),
    output_dir = file.path(results_dir, "figures"),
    pves = c(.025, .05), zoom_max_fpr = .25,
    exclude_nonconverged = FALSE, write_png = TRUE) {
  stopifnot(length(pves) > 0, all(is.finite(pves)), all(pves > 0 & pves < 1),
    !anyDuplicated(pves), length(zoom_max_fpr) == 1L,
    is.finite(zoom_max_fpr), zoom_max_fpr > 0, zoom_max_fpr <= 1)
  files <- list.files(results_dir, pattern = "^additive_init_chunk[0-9]+\\.rds$",
                      full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("No initialization checkpoints found in: ", results_dir,
    "\nExpected pve_0.025/chunks and pve_0.05/chunks subfolders.")
  thresholds <- c(additive_roc_thresholds, Inf)
  totals <- array(0, dim = c(length(thresholds), 2L, length(additive_roc_methods), length(pves)))
  positives <- negatives <- included <- numeric(length(pves))
  audit <- list()
  reference <- NULL
  seen <- new.env(hash = TRUE, parent = emptyenv())
  for (file in files) {
    saved <- readRDS(file)
    if (!is.list(saved$settings) || !is.list(saved$results)) stop("Invalid checkpoint: ", file)
    s <- saved$settings
    j <- match(s$pve, pves)
    if (length(j) != 1L || is.na(j)) next
    design <- s
    design$chunk <- design$pve <- NULL
    # Different PLINK workspace sizes do not change the simulation design.
    if (identical(design$genotype_inputs$mode, "gtex_fallback"))
      design$genotype_inputs$plink_memory <- NULL
    if (is.null(reference)) reference <- design
    if (!identical(design, reference)) stop("Incompatible simulation settings in: ", file)
    message("Reading PVE ", 100 * s$pve, "%: ", basename(file))
    for (i in seq_along(saved$results)) {
      r <- saved$results[[i]]
      seed <- if (is.null(r)) NA_real_ else r$seed
      status <- "pending"
      detail <- ""
      if (!is.null(r)) {
        if (length(seed) != 1L || !is.finite(seed)) stop("Invalid replicate seed in: ", file)
        key <- paste(s$pve, seed, sep = "|")
        if (exists(key, envir = seen, inherits = FALSE)) stop("Duplicate PVE/seed across checkpoints: ", key)
        assign(key, TRUE, envir = seen)
        if (!is.null(r$error)) {
          status <- "failed"
          detail <- r$error
        } else {
          if (length(r$p) != 1L || !is.finite(r$p) || r$p != floor(r$p) || r$p < 2L ||
              !all(additive_roc_methods %in% names(r$fits))) stop("Invalid saved fits in: ", file, "; seed ", seed)
          if (length(r$true_pos) != s$K) stop("Causal count disagrees with checkpoint settings: ", file)
          counts <- lapply(additive_roc_methods, function(method) {
            pip <- r$fits[[method]]$pip
            if (length(pip) != r$p) stop("PIP length disagrees with SNP count: ", method, "; seed ", seed)
            additive_roc_counts(pip, r$true_pos)
          })
          converged <- vapply(r$fits[additive_roc_methods], function(f) isTRUE(f$converged), logical(1))
          detail <- paste(additive_roc_methods[!converged], collapse = "; ")
          if (exclude_nonconverged && !all(converged)) {
            status <- "excluded_nonconverged"
          } else {
            status <- if (all(converged)) "included" else "included_nonconverged"
            for (m in seq_along(additive_roc_methods)) totals[, , m, j] <- totals[, , m, j] + counts[[m]]
            positives[j] <- positives[j] + length(r$true_pos)
            negatives[j] <- negatives[j] + r$p - length(r$true_pos)
            included[j] <- included[j] + 1
          }
        }
      }
      audit[[length(audit) + 1L]] <- data.frame(pve = s$pve, chunk = s$chunk,
        slot = i, seed = seed, status = status, detail = detail, file = file)
    }
    rm(saved)
  }
  if (!length(audit)) stop("No checkpoints for the requested PVE values in: ", results_dir)
  audit <- do.call(rbind, audit)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(audit, file.path(output_dir, "roc_pip_replicates.csv"), row.names = FALSE)
  if (!any(included > 0)) stop("No successful replicates available for ROC plotting. See roc_pip_replicates.csv in ", output_dir)
  if (any(included > 0 & (positives <= 0 | negatives <= 0))) stop("ROC requires both causal and noncausal SNPs.")
  rows <- list()
  for (j in which(included > 0)) for (m in seq_along(additive_roc_methods)) {
    rows[[length(rows) + 1L]] <- data.frame(pve = pves[j], method = additive_roc_methods[m],
      threshold = thresholds, tp = totals[, 1, m, j], fp = totals[, 2, m, j],
      n_causal = positives[j], n_noncausal = negatives[j], n_replicates = included[j],
      tpr = totals[, 1, m, j] / positives[j], fpr = totals[, 2, m, j] / negatives[j])
  }
  curves <- do.call(rbind, rows)
  write.csv(curves, file.path(output_dir, "roc_pip_counts.csv"), row.names = FALSE)

  draw <- function(max_fpr) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(1, length(pves)), mar = c(4.2, 4.4, 2.4, 1),
        oma = c(3.8, .3, 3.4, .3), mgp = c(2.6, .7, 0), tcl = -.25, family = "sans")
    for (j in seq_along(pves)) {
      plot(NA, xlim = c(0, max_fpr), ylim = c(0, 1), xaxs = "i", yaxs = "i",
        axes = FALSE, xlab = "False positive rate", ylab = "True positive rate (power)")
      xticks <- pretty(c(0, max_fpr), 5)
      xticks <- xticks[xticks >= 0 & xticks <= max_fpr]
      yticks <- seq(0, 1, .2)
      abline(v = xticks, h = yticks, col = "#E2E2E2", lwd = .7)
      abline(0, 1, col = "#888888", lty = 3)
      axis(1, at = xticks, col = "#777777")
      axis(2, at = yticks, las = 1, col = "#777777")
      box(col = "#777777")
      title(main = sprintf("PVE = %g%%", 100 * pves[j]), line = 1.1, cex.main = 1.05)
      mtext(sprintf("%d successful replicates", included[j]), side = 3, line = .1, cex = .8)
      if (included[j] == 0) text(max_fpr / 2, .5, "No successful replicates", col = "#777777")
      for (m in seq_along(additive_roc_methods)) {
        d <- curves[curves$pve == pves[j] & curves$method == additive_roc_methods[m], ]
        if (!nrow(d)) next
        d <- d[order(d$threshold, decreasing = TRUE), ]
        lines(d$fpr, d$tpr, col = additive_roc_colors[m], lty = m, lwd = 2)
      }
    }
    mtext("Detection of causal SNPs using PIPs", side = 3, outer = TRUE, line = 1.8, font = 2, cex = 1.15)
    mtext(sprintf("Additive effects | %d causal SNPs | n = %d | fitted L = %d",
                  reference$K, reference$n, reference$fit_L),
          side = 3, outer = TRUE, line = .4, cex = .85)
    par(fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4), new = TRUE)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    legend(.5, .075, legend = additive_roc_methods, col = additive_roc_colors,
      lty = seq_along(additive_roc_methods), lwd = 2, horiz = TRUE,
      xjust = .5, yjust = .5, bty = "n", cex = if (length(pves) > 1) .95 else .68)
    text(.5, .025, "Pooled SNPs; PIP >= threshold; same replicates for all methods", cex = .75)
  }
  outputs <- character()
  save_plot <- function(name, max_fpr) {
    path <- file.path(output_dir, paste0(name, ".pdf"))
    grDevices::pdf(path, width = 5.5 * length(pves), height = 5.7, useDingbats = FALSE)
    tryCatch(draw(max_fpr), finally = grDevices::dev.off())
    outputs <<- c(outputs, path)
    if (write_png) {
      path <- file.path(output_dir, paste0(name, ".png"))
      png_args <- list(filename = path, width = 5.5 * length(pves), height = 5.7, units = "in", res = 200)
      if (capabilities("cairo")) png_args$type <- "cairo"
      do.call(grDevices::png, png_args)
      tryCatch(draw(max_fpr), finally = grDevices::dev.off())
      outputs <<- c(outputs, path)
    }
  }
  save_plot("roc_pip", 1)
  save_plot("roc_pip_zoom", zoom_max_fpr)
  message("ROC plots and CSVs saved in: ", normalizePath(output_dir, winslash = "/"))
  print(table(PVE = audit$pve, status = audit$status))
  invisible(list(curves = curves, replicates = audit, files = outputs))
}

if (isTRUE(getOption("susie.init.roc.autorun", TRUE)))
  additive_init_roc <- plot_additive_init_roc()
