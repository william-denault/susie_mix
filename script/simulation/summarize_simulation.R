#!/usr/bin/env Rscript
# Base-R aggregation: all main comparisons use paired, converged replicates.

simulation_summarize <- function(output, allow_incomplete = FALSE) {
  manifests <- list.files(output, "^manifest_shard_[0-9]+\\.rds$", full.names = TRUE)
  if (!length(manifests)) stop("No run manifests found in ", output)
  manifest <- readRDS(manifests[1])
  if (any(vapply(manifests, function(p) readRDS(p)$signature != manifest$signature,
                 logical(1)))) stop("Incompatible manifests in this output directory.")
  config <- manifest$config
  scenarios <- manifest$scenarios
  files <- list.files(file.path(output, "replicates"), "^rep_[0-9]+\\.rds$",
                       full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("No replication checkpoints found.")
  summary_dir <- file.path(output, "summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  metrics <- sets <- curves <- calibration <- statuses <- variances <- rejected <- list()
  method_status <- list()
  for (i in seq_along(files)) {
    record <- readRDS(files[i])
    if (!identical(record$signature, manifest$signature)) stop("Incompatible checkpoint: ", files[i])
    meta <- cbind(record$scenario, replication = record$replication)
    rownames(meta) <- NULL
    statuses[[i]] <- cbind(meta, status = record$status,
      gene = if (is.null(record$gene)) NA_character_ else record$gene$gene_name,
      error = if (is.null(record$error)) "" else record$error)
    if (!is.null(record$phenotype))
      variances[[i]] <- cbind(meta, as.data.frame(as.list(record$phenotype$variance)))
    if (length(record$rejected_loci))
      rejected[[i]] <- cbind(meta, do.call(rbind, record$rejected_loci))
    for (method in names(record$methods)) {
      result <- record$methods[[method]]
      m <- cbind(meta, method = method)
      rownames(m) <- NULL
      method_status[[length(method_status) + 1L]] <- cbind(m,
        converged = isTRUE(result$converged), seconds = result$seconds,
        niter = if (is.null(result$niter)) NA_integer_ else result$niter,
        error = if (is.null(result$error)) "" else result$error,
        warnings = paste(result$warnings, collapse = " | "))
      if (record$status != "ok") next
      key <- paste(record$scenario$scenario_id, record$replication, method, sep = "_")
      metrics[[key]] <- cbind(m, result$metrics)
      if (nrow(result$cs)) sets[[key]] <- cbind(m, result$cs)
      curves[[key]] <- cbind(m, result$curve)
      calibration[[key]] <- cbind(m, result$calibration)
    }
  }
  bind <- function(x) if (length(Filter(Negate(is.null), x))) do.call(rbind, x) else data.frame()
  write_table <- function(x, name) {
    rownames(x) <- NULL
    write.csv(x, file.path(summary_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
  }
  statuses <- bind(statuses)
  if (anyDuplicated(paste(statuses$scenario_id, statuses$replication)))
    stop("Duplicate scenario/replication checkpoints.")
  completion <- do.call(rbind, lapply(scenarios$scenario_id, function(id) {
    s <- statuses[statuses$scenario_id == id, , drop = FALSE]
    cbind(scenarios[scenarios$scenario_id == id, , drop = FALSE],
      target = config$replications, attempted = nrow(s),
      paired_success = sum(s$status == "ok"),
      failed = sum(s$status != "ok"), missing = config$replications - nrow(s))
  }))
  complete <- all(completion$paired_success >= config$replications) &&
    all(completion$failed == 0L) && all(completion$missing == 0L)
  write_table(completion, "completion")
  write_table(statuses, "replicate_status")
  write_table(bind(method_status), "method_status")
  write_table(bind(variances), "phenotype_variance")
  write_table(bind(rejected), "rejected_loci")
  writeLines(c(
    paste("Mode:", if (config$pilot) "PILOT - not scientific benchmark results" else "PRODUCTION"),
    paste("All scenarios complete:", complete),
    paste("Paired converged replicates:", sum(completion$paired_success)),
    paste("Required datasets:", nrow(scenarios) * config$replications),
    "Main metrics exclude both methods whenever either fit fails or does not converge.",
    "Inspect completion.csv and method_status.csv before interpreting the figures."),
    file.path(summary_dir, "STATUS.txt"))
  if (!length(metrics)) stop("No paired converged fits. Status tables were written.")
  metrics <- bind(metrics)
  sets <- bind(sets)
  curves <- bind(curves)
  calibration <- bind(calibration)
  write_table(metrics, "replicate_metrics")
  write_table(sets, "credible_sets")
  # Standard errors are across independent replicates (loci), not across correlated SNPs.
  mcse <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_
  }
  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  group_rows <- function(d, columns) split(seq_len(nrow(d)),
    interaction(d[, columns, drop = FALSE], drop = TRUE, lex.order = TRUE))
  group_cols <- c("scenario_id", "architecture", "k", "pve", "method")
  metric_names <- c("power", "fdp", "cs_power", "cs_coverage", "cs_coding_coverage",
    "cs_size_snp", "cs_size_predictor", "cs_purity_snp", "cs_purity_predictor", "brier_score", "n_cs")
  metric_summary <- do.call(rbind, lapply(group_rows(metrics, c("scenario_id", "method")), function(ix) {
    d <- metrics[ix, , drop = FALSE]
    do.call(rbind, lapply(metric_names, function(nm) cbind(d[1, group_cols, drop = FALSE],
      metric = nm, n_replicates = nrow(d), n_defined = sum(is.finite(d[[nm]])),
      mean = finite_mean(d[[nm]]), mcse = mcse(d[[nm]]))))
  }))
  write_table(metric_summary, "metrics_summary")
  # Paired differences quantify the comparison with its Monte Carlo uncertainty.
  paired <- do.call(rbind, lapply(split(metrics, metrics$scenario_id), function(d) {
    a <- d[d$method == "susie", ]
    b <- d[d$method == "susie_mix", ]
    b <- b[match(a$replication, b$replication), ]
    do.call(rbind, lapply(metric_names, function(nm) {
      delta <- b[[nm]] - a[[nm]]
      cbind(a[1, c("scenario_id", "architecture", "k", "pve"), drop = FALSE], metric = nm,
        n_defined = sum(is.finite(delta)), mix_minus_susie = finite_mean(delta), mcse = mcse(delta))
    }))
  }))
  write_table(paired, "paired_differences")
  # Pooled CS coverage/size also reported; per-replicate versions above give equal
  # weight to each locus with >=1 reported CS. Neither assigns coverage=0 to no-CS fits.
  if (nrow(sets)) {
    cs_summary <- do.call(rbind, lapply(group_rows(sets, c("scenario_id", "method")), function(ix) {
      d <- sets[ix, ]
      cbind(d[1, group_cols, drop = FALSE], n_sets = nrow(d),
        coverage = mean(d$covered_snp), coding_coverage = mean(d$covered_coding),
        mean_size_snp = mean(d$size_snp), median_size_snp = median(d$size_snp),
        mean_size_predictor = mean(d$size_predictor),
        mean_purity_snp = mean(d$purity_snp),
        mean_purity_predictor = mean(d$purity_predictor),
        fraction_purity_snp_approximate = mean(d$purity_snp_approximate))
    }))
    write_table(cs_summary, "cs_summary")
  }
  curve_summary <- do.call(rbind, lapply(group_rows(curves, c("scenario_id", "method", "threshold")), function(ix) {
    d <- curves[ix, ]
    cbind(d[1, c(group_cols, "threshold"), drop = FALSE], n_replicates = nrow(d),
      power = mean(d$power), power_mcse = mcse(d$power),
      fdr = mean(d$fdp), fdr_mcse = mcse(d$fdp),
      pooled_fdp = sum(d$fp) / max(1, sum(d$discoveries)),
      tp = sum(d$tp), fp = sum(d$fp), discoveries = sum(d$discoveries))
  }))
  write_table(curve_summary, "power_fdr")
  # Ratios use pooled counts; uncertainty uses independent replicate clusters.
  ratio_mcse <- function(numerator, denominator) {
    if (length(denominator) < 2L || sum(denominator) == 0) return(NA_real_)
    ratio <- sum(numerator) / sum(denominator)
    sd(numerator - ratio * denominator) / sqrt(length(denominator)) / mean(denominator)
  }
  cal_summary <- do.call(rbind, lapply(group_rows(calibration, c("scenario_id", "method", "bin")), function(ix) {
    d <- calibration[ix, ]
    n <- sum(d$n)
    cbind(d[1, c(group_cols, "bin", "lower", "upper"), drop = FALSE],
      n = n, n_nonempty_replicates = sum(d$n > 0), n_causal = sum(d$n_causal),
      mean_pip = if (n) sum(d$sum_pip) / n else NA_real_,
      causal_fraction = if (n) sum(d$n_causal) / n else NA_real_,
      causal_fraction_mcse = ratio_mcse(d$n_causal, d$n),
      calibration_gap = if (n) sum(d$n_causal - d$sum_pip) / n else NA_real_,
      calibration_gap_mcse = ratio_mcse(d$n_causal - d$sum_pip, d$n))
  }))
  write_table(cal_summary, "pip_calibration")
  simulation_plots(metric_summary, curve_summary, cal_summary, summary_dir,
    pilot = config$pilot, complete = complete, pip_threshold = config$pip_threshold,
    coverage = config$coverage)
  message("Summaries and PNG figures written to ", summary_dir)
  if (!complete && !allow_incomplete)
    stop("Benchmark incomplete. See completion.csv. Use --allow-incomplete only to inspect partial results.")
  invisible(list(completion = completion, metrics = metric_summary,
                  curve = curve_summary, calibration = cal_summary))
}

simulation_plots <- function(metrics, curves, calibration, summary_dir, pilot, complete,
                              pip_threshold, coverage) {
  figure_dir <- file.path(summary_dir, "figures")
  dir.create(figure_dir, showWarnings = FALSE)
  architectures <- unique(metrics$architecture)
  colors <- c(susie = "#2463A5", susie_mix = "#D45F22")
  label <- paste(if (pilot) "PILOT |" else "", if (!complete) "INCOMPLETE |" else "")
  panels <- function(path, draw) {
    png(path, width = 2100, height = 1050, res = 140)
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(2, 4), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2, 0))
    for (a in architectures) draw(a)
    plot.new()
    legend("center", c("SuSiE", "SuSiE-mix"), col = colors, lty = 1, lwd = 2, bty = "n")
    mtext(label, outer = TRUE, line = 0, cex = 0.9)
  }
  metric_labels <- c(power = paste0("Power (SNP PIP >= ", pip_threshold, ")"), cs_power = "Causal SNPs in any CS",
    cs_coverage = "CS coverage (biological SNP)", cs_size_snp = "Mean CS size (unique SNPs)",
    cs_purity_snp = "Mean CS purity (additive LD)", cs_purity_predictor = "Mean CS purity (fitted coding)")
  for (h in sort(unique(metrics$pve))) {
    tag <- sprintf("pve_%02d", round(100 * h))
    for (nm in names(metric_labels)) {
      d <- metrics[metrics$pve == h & metrics$metric == nm, ]
      panels(file.path(figure_dir, paste0(nm, "_", tag, ".png")), function(a) {
        z <- d[d$architecture == a, ]
        ymax <- if (nm == "cs_size_snp")
          max(c(1, z$mean, z$mean + 1.96 * z$mcse), na.rm = TRUE) else 1
        plot(NA, xlim = c(1, 5), ylim = c(0, ymax), xlab = "Number of causal SNPs",
             ylab = metric_labels[[nm]], main = paste(gsub("_", " + ", a), "|", 100 * h, "%"))
        if (nm == "cs_coverage") abline(h = coverage, lty = 3, col = "grey55")
        for (method in names(colors)) {
          q <- z[z$method == method, ]
          q <- q[order(q$k), ]
          lines(q$k, q$mean, type = "b", col = colors[[method]], pch = 16, lwd = 2)
          use <- is.finite(q$mean) & is.finite(q$mcse) & q$mcse > 0
          if (any(use)) arrows(q$k[use], pmax(0, q$mean[use] - 1.96 * q$mcse[use]),
            q$k[use], pmin(ymax, q$mean[use] + 1.96 * q$mcse[use]),
            angle = 90, code = 3, length = 0.035, col = colors[[method]])
        }
        if (!any(is.finite(z$mean)))
          text(3, ymax / 2, if (nrow(z)) "No defined estimates (no reported CS)" else
            "No completed replicates", cex = 0.8)
      })
    }
    for (k in 1:5) {
      panels(file.path(figure_dir, paste0("power_fdr_", tag, "_k", k, ".png")), function(a) {
        z <- curves[curves$pve == h & curves$k == k & curves$architecture == a, ]
        plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Empirical FDR (mean FDP)",
          ylab = "Power", main = paste(gsub("_", " + ", a), "| K =", k))
        for (method in names(colors)) {
          q <- z[z$method == method, ]
          q <- q[order(q$threshold, decreasing = TRUE), ]
          lines(q$fdr, q$power, col = colors[[method]], lwd = 2)
        }
        if (!nrow(z)) text(0.5, 0.5, "No scenario / completed replicates", cex = 0.8)
      })
      panels(file.path(figure_dir, paste0("pip_calibration_", tag, "_k", k, ".png")), function(a) {
        z <- calibration[calibration$pve == h & calibration$k == k & calibration$architecture == a, ]
        plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Mean SNP PIP in bin",
          ylab = "Observed causal fraction", main = paste(gsub("_", " + ", a), "| K =", k))
        abline(0, 1, lty = 3, col = "grey55")
        for (method in names(colors)) {
          q <- z[z$method == method & z$n > 0, ]
          q <- q[order(q$bin), ]
          lines(q$mean_pip, q$causal_fraction, type = "b", pch = 16,
                col = colors[[method]], lwd = 2)
        }
        if (!nrow(z)) text(0.5, 0.5, "No scenario / completed replicates", cex = 0.8)
      })
    }
  }
}

if (sys.nframe() == 0L) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "../.."),
                        winslash = "/", mustWork = TRUE)
  args <- commandArgs(TRUE)
  allow <- "--allow-incomplete" %in% args
  args <- setdiff(args, "--allow-incomplete")
  if (!length(args)) {
    output <- file.path(root, "simulation results", "production")
  } else if (length(args) == 2L && args[1] == "--output") {
    output <- args[2]
  } else stop("Usage: summarize_simulation.R [--output PATH] [--allow-incomplete]")
  simulation_summarize(output, allow)
}
