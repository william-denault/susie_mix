# Full posterior CDFs of slider delta. Run with Source in RStudio or Rscript.
# Base R only. Reads one gene RDS at a time; no refitting is performed.

# Settings ----
project <- Sys.getenv("SUSIE_MIX_PROJECT_DIR",
                      "/project2/mstephens/wdenault/susie_mix")
root <- file.path(project, "results_slide_prior_em")
iteration <- NULL                 # Latest iteration with saved gene RDS files.
# iteration <- "iteration_024"    # Or select an iteration explicitly.
unforced_leads_only <- TRUE       # Same filter as your existing delta plots.
min_entries <- 3L                 # Same length(out) > 2 filter; use 0 for all.

# Extract the full probability distributions ----
read_slider_posteriors <- function(root, iteration = NULL,
                                   unforced_leads_only = TRUE, min_entries = 3L) {
  result_files <- function(i) list.files(file.path(root, i, "results"),
    pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  if (is.null(iteration)) {
    candidates <- list.dirs(root, recursive = FALSE, full.names = FALSE)
    candidates <- candidates[grepl("^iteration_[0-9]+$", candidates)]
    candidates <- candidates[order(as.numeric(sub("iteration_", "", candidates)),
                                   decreasing = TRUE)]
    for (candidate in candidates) {
      if (length(result_files(candidate))) {
        iteration <- candidate
        break
      }
    }
    if (is.null(iteration)) stop("No iteration contains gene RDS files under ", root)
  }
  files <- result_files(iteration)
  if (!length(files)) stop("No gene RDS files in ", file.path(root, iteration, "results"))
  message("Reading ", iteration, ": ", length(files), " saved gene files.")
  message("This is a snapshot of saved fits; iteration completeness is not checked.")
  # In online EM, these fits can use different batch priors. They are the saved
  # posteriors, not a new E-step under the last prior_history.csv weights.
  rows <- lead_prob <- cs_prob <- list()
  grid <- NULL
  n <- filtered_files <- forced_leads <- 0L
  for (k in seq_along(files)) {
    tryCatch({
      out <- readRDS(files[k])
      if (!is.list(out)) stop("Expected a list of tissue results.")
      if (!is.null(out$error)) stop("Saved gene error: ", paste(out$error, collapse = "; "))
      if (length(out) < min_entries) {
        filtered_files <- filtered_files + 1L
        next
      }
      for (o in seq_along(out)) {
        if (!is.list(out[[o]])) next
        fit <- out[[o]]$fit_slide_prior
        if (!length(fit$sets$cs)) next
        tissue <- names(out)[o]
        if (!length(tissue) || is.na(tissue) || !nzchar(tissue)) tissue <- as.character(o)
        if (isFALSE(fit$converged)) stop("Unconverged fit in ", tissue)
        g <- fit$delta_grid
        if (!is.matrix(fit$alpha) || !is.numeric(g) || !length(g) ||
            any(!is.finite(g)) || anyDuplicated(g))
          stop("Missing/invalid alpha or delta_grid in ", tissue)
        if (!is.numeric(fit$alpha_delta) ||
            !identical(dim(fit$alpha_delta), c(dim(fit$alpha), length(g))))
          stop("Missing/incompatible alpha_delta in ", tissue,
               ". The full CDF requires the saved joint posterior; ",
               "delta_summary.rds cannot reconstruct it.")
        ord <- order(g)
        g <- g[ord]
        if (is.null(grid)) grid <- g
        if (!identical(as.numeric(g), as.numeric(grid)))
          stop("Different delta grids across saved fits.")
        p <- ncol(fit$alpha)
        if (!is.logical(fit$delta_forced) || length(fit$delta_forced) != p ||
            anyNA(fit$delta_forced)) stop("Invalid delta_forced in ", tissue)
        components <- fit$sets$cs_index
        if (length(components) != length(fit$sets$cs) || anyNA(components) ||
            any(components < 1 | components > nrow(fit$alpha) |
                  components != floor(components))) stop("Invalid cs_index in ", tissue)
        for (s in seq_along(fit$sets$cs)) {
          cs <- fit$sets$cs[[s]]
          ell <- components[s]
          if (!length(cs) || anyNA(cs) || anyDuplicated(cs) ||
              any(cs < 1 | cs > p | cs != floor(cs))) stop("Invalid CS members in ", tissue)
          a <- fit$alpha[ell, cs]
          if (any(!is.finite(a)) || any(a < 0) || sum(a) <= 0)
            stop("Invalid CS alpha in ", tissue)
          lead <- cs[which.max(a)]
          if (unforced_leads_only && fit$delta_forced[lead]) {
            forced_leads <- forced_leads + 1L
            next
          }
          # alpha_delta[ell, j, d] = P(J_ell = j, delta_ell = grid[d] | data).
          # Sum over CS SNPs, then normalize: condition on J_ell being in this CS.
          # For the lead curve, condition on J_ell being the lead SNP instead.
          joint <- matrix(fit$alpha_delta[ell, cs, ord, drop = FALSE],
                          nrow = length(cs), ncol = length(grid))
          if (any(!is.finite(joint)) || any(joint < 0) ||
              any(abs(rowSums(joint) - a) > 1e-7))
            stop("Joint posterior does not agree with alpha in ", tissue)
          w_cs <- colSums(joint)
          w_lead <- joint[match(lead, cs), ]
          if (sum(w_cs) <= 0 || sum(w_lead) <= 0) stop("Zero posterior mass in ", tissue)
          cs_name <- names(fit$sets$cs)[s]
          if (!length(cs_name) || is.na(cs_name) || !nzchar(cs_name)) cs_name <- paste0("L", ell)
          n <- n + 1L
          rows[[n]] <- data.frame(file = basename(files[k]), tissue = tissue,
            cs = cs_name, component = ell, n_snps = length(cs),
            cs_coverage = sum(a), lead_index = lead, lead_alpha = max(a),
            lead_delta_forced = fit$delta_forced[lead],
            cs_forced_mass = sum(a * fit$delta_forced[cs]) / sum(a))
          cs_prob[[n]] <- w_cs / sum(w_cs)
          lead_prob[[n]] <- w_lead / sum(w_lead)
        }
      }
    }, error = function(e) stop(basename(files[k]), ": ", conditionMessage(e), call. = FALSE))
    if (k %% 100L == 0L || k == length(files)) message(k, "/", length(files), " files; ", n, " retained CSs")
  }
  if (!n) stop("No credible sets remain after filtering.")
  signals <- do.call(rbind, rows)
  lead_pmf <- do.call(rbind, lead_prob)
  cs_pmf <- do.call(rbind, cs_prob)
  colnames(lead_pmf) <- colnames(cs_pmf) <- as.character(grid)
  to_cdf <- function(m) {
    ans <- m
    if (ncol(m) > 1L) for (j in 2:ncol(m)) ans[, j] <- ans[, j - 1L] + m[, j]
    ans[, ncol(m)] <- 1
    ans
  }
  # Rows of these matrices align exactly with the signals data frame.
  lead_cdf <- to_cdf(lead_pmf)
  cs_cdf <- to_cdf(cs_pmf)
  signals$lead_delta_mean <- as.numeric(lead_pmf %*% grid)
  signals$cs_delta_mean <- as.numeric(cs_pmf %*% grid)
  # Each retained CS has equal weight, including in the all-tissue curve.
  # This describes selected signals; it is not a population CDF or an ECDF of means.
  groups <- c(list(seq_len(n)), unname(split(seq_len(n), signals$tissue)))
  labels <- c("All tissues", names(split(seq_len(n), signals$tissue)))
  tissue_cdf <- do.call(rbind, lapply(seq_along(groups), function(i) {
    ix <- groups[[i]]
    data.frame(tissue = labels[i], pooled = i == 1L, n_cs = length(ix),
      n_genes = length(unique(signals$file[ix])), delta = grid,
      lead_pmf = colMeans(lead_pmf[ix, , drop = FALSE]),
      cs_pmf = colMeans(cs_pmf[ix, , drop = FALSE]),
      lead_cdf = colMeans(lead_cdf[ix, , drop = FALSE]),
      cs_cdf = colMeans(cs_cdf[ix, , drop = FALSE]))
  }))
  message("Excluded ", filtered_files, " gene files by min_entries and ",
          forced_leads, " CSs with a forced-additive lead.")
  list(iteration = iteration, root = root, n_files = length(files),
       unforced_leads_only = unforced_leads_only, min_entries = min_entries,
       filtered_files = filtered_files, forced_leads_excluded = forced_leads,
       grid = grid, signals = signals, lead_pmf = lead_pmf, cs_pmf = cs_pmf,
       lead_cdf = lead_cdf, cs_cdf = cs_cdf, tissue_cdf = tissue_cdf)
}

# Plot tissue averages of the full posterior CDFs ----
plot_slider_posterior_cdf <- function(x, which = c("lead", "cs")) {
  which <- match.arg(which)
  d <- x$tissue_cdf
  tissues <- unique(d$tissue[!d$pooled])
  nr <- ceiling(length(tissues) / min(5L, length(tissues)))
  old <- par(mfrow = c(nr, min(5L, length(tissues))),
             mar = c(3.5, 3.8, 2.8, 0.7), oma = c(0, 0, 4, 0), mgp = c(2.3, .7, 0))
  on.exit(par(old))
  pooled <- d[d$pooled, ]
  # Pad the support to show the jumps at its endpoints, including delta = -1.
  pad <- max(diff(range(x$grid)) * .025, .025)
  xx <- c(min(x$grid) - pad, x$grid, max(x$grid) + pad)
  for (tissue in tissues) {
    a <- d[!d$pooled & d$tissue == tissue, ]
    yy <- a[[paste0(which, "_cdf")]]
    plot(xx, c(0, yy, 1), type = "s", ylim = c(0, 1), lwd = 2,
         main = paste0(tissue, "\n", a$n_cs[1], " CSs; ", a$n_genes[1], " genes"),
         cex.main = .85, xlab = expression(delta),
         ylab = expression("Mean posterior " * P(delta <= x)))
    abline(v = 0, col = "grey80")
    lines(xx, c(0, pooled[[paste0(which, "_cdf")]], 1), type = "s",
          col = "#C33B32", lty = 2, lwd = 1.8)
  }
  title <- if (which == "lead") "Lead-SNP posterior CDF" else
    "Credible-set posterior CDF"
  mtext(paste(x$iteration, title, sep = " - "), outer = TRUE, line = 2.4, cex = 1)
  mtext(paste0("Black: tissue; red dashed: all tissues pooled; equal weight per CS",
               if (x$unforced_leads_only) "; unforced leads only" else "; all leads"),
        outer = TRUE, line = .9, cex = .85)
  invisible(x)
}

# Generate and save ----
posterior_cdf <- read_slider_posteriors(root, iteration,
                                       unforced_leads_only, min_entries)
output_dir <- file.path(root, posterior_cdf$iteration, "posterior_cdf")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(posterior_cdf, file.path(output_dir, "posterior_cdf.rds"))
write.csv(posterior_cdf$tissue_cdf, file.path(output_dir, "posterior_cdf_by_tissue.csv"),
          row.names = FALSE)
write.csv(posterior_cdf$signals, file.path(output_dir, "posterior_cdf_signals.csv"),
          row.names = FALSE)
save_cdf_plot <- function(which) {
  nt <- length(unique(posterior_cdf$signals$tissue))
  width <- max(6, min(5L, nt) * 3)
  height <- ceiling(nt / min(5L, nt)) * 2.7 + .8
  for (ext in c("pdf", "png")) {
    filename <- file.path(output_dir, paste0(which, "_posterior_cdf.", ext))
    if (ext == "pdf") pdf(filename, width = width, height = height) else
      png(filename, width = width, height = height, units = "in", res = 180)
    tryCatch(plot_slider_posterior_cdf(posterior_cdf, which), finally = dev.off())
  }
}
save_cdf_plot("lead")
save_cdf_plot("cs")
if (interactive()) plot_slider_posterior_cdf(posterior_cdf, "lead")
message("Saved posterior CDFs to ", output_dir)

# RStudio examples after running this script:
# plot_slider_posterior_cdf(posterior_cdf, "cs")
# posterior_cdf$signals[1, ]             # Gene/tissue/CS for matrix row 1.
# plot(posterior_cdf$grid, posterior_cdf$lead_cdf[1, ], type = "s",
#      ylim = c(0, 1), xlab = "delta", ylab = "Posterior CDF")
# At arbitrary x, F(x) = sum(posterior_cdf$lead_pmf[1, posterior_cdf$grid <= x]).
