# Compact sufficient statistics and analytic intervals for simulation metrics.
cs_summary <- function(sets, cs_snps, truth) {
  n_cs <- length(cs_snps)
  if (length(sets$cs) != n_cs) stop("Saved CS mapping has the wrong length.")
  hit <- vapply(cs_snps, function(cs) any(cs %in% truth), logical(1))
  purity <- sets$purity$min.abs.corr
  if (n_cs && (length(purity) != n_cs || any(!is.finite(purity)))) {
    stop("Missing or invalid minimum-correlation purity for a reported CS.")
  }
  sizes <- vapply(cs_snps, function(cs) length(unique(cs)), integer(1))
  c(n_cs = n_cs, covered_cs = sum(hit), purity_sum = sum(purity),
    cs_size_sum = sum(sizes), cs_size_sum_sq = sum(as.double(sizes)^2),
    recovered = sum(unique(truth) %in% unlist(cs_snps, use.names = FALSE)),
    n_causal = length(unique(truth)))
}

summarize_metrics <- function(replicates, level = .95,
                              methods = unique(replicates$method),
                              proportion_n = c("denominator", "replicates")) {
  proportion_n <- match.arg(proportion_n)
  stopifnot(level > 0, level < 1)
  if (!"cs_size_sum_sq" %in% names(replicates)) {
    stop("Saved CS-size sums of squares are missing; rerun plot_simulations.R to rebuild the compact summaries.")
  }
  z <- qnorm((1 + level) / 2)
  group_vars <- c("scenario", "pve", "K", "method")
  cell_vars <- c("scenario", "pve", "K")
  numerators <- c(coverage = "covered_cs", purity = "purity_sum",
                  power = "recovered", cs_size = "cs_size_sum")
  denominators <- c(coverage = "n_cs", purity = "n_cs", power = "n_causal", cs_size = "n_cs")
  count_vars <- unique(c(unname(numerators), unname(denominators), "cs_size_sum_sq"))
  groups <- split(seq_len(nrow(replicates)),
                  interaction(replicates[cell_vars], drop = TRUE, sep = "|"))
  summary_rows <- difference_rows <- list()
  for (ids in groups) {
    cell <- replicates[ids, , drop = FALSE]
    seeds <- sort(unique(cell$seed))
    n_blocks <- length(seeds)
    estimates <- influences <- list()
    reference_keys <- NULL
    for (method in methods) {
      d <- cell[cell$method == method, , drop = FALSE]
      if (!nrow(d)) stop("All selected methods are required for paired comparisons.")
      keys <- sort(paste(d$configuration, d$seed, sep = "|"))
      if (anyDuplicated(keys)) stop("Duplicate configuration/seed within a method.")
      if (is.null(reference_keys)) reference_keys <- keys
      if (!identical(keys, reference_keys)) stop("Methods must use the same replicates for paired comparisons.")
      totals_by_seed <- rowsum(as.matrix(d[count_vars]), group = d$seed)
      totals <- totals_by_seed[match(as.character(seeds), rownames(totals_by_seed)), , drop = FALSE]
      point <- colSums(totals)
      row <- data.frame(d[1, group_vars, drop = FALSE], as.list(point),
        n_replicates = nrow(d), n_configurations = length(unique(d$configuration)),
        n_seed_blocks = n_blocks, interval_level = level,
        proportion_n_definition = proportion_n, row.names = NULL)
      estimates[[method]] <- numeric()
      influences[[method]] <- list()
      for (metric in names(numerators)) {
        a <- numerators[[metric]]; b <- denominators[[metric]]
        denominator <- unname(point[b])
        estimate <- if (denominator > 0) unname(point[a] / denominator) else NA_real_
        n <- if (metric != "cs_size" && proportion_n == "replicates") nrow(d) else denominator
        se <- NA_real_
        if (is.finite(estimate) && n > 0) {
          if (metric == "cs_size") {
            # Sample variance of individual CS sizes, not of replicate means.
            if (n > 1) {
              variance <- max(0, (point["cs_size_sum_sq"] - point[a]^2 / n) / (n - 1))
              se <- unname(sqrt(variance / n))
            }
          } else {
            # Requested approximation also applies to continuous mean purity.
            se <- sqrt(max(0, estimate * (1 - estimate)) / n)
          }
        }
        row[[metric]] <- estimate
        row[[paste0(metric, "_n")]] <- n
        row[[paste0(metric, "_se")]] <- se
        row[[paste0(metric, "_lo")]] <- max(0, estimate - z * se)
        row[[paste0(metric, "_hi")]] <- if (metric == "cs_size") estimate + z * se else
          min(1, estimate + z * se)
        estimates[[method]][metric] <- estimate
        influences[[method]][[metric]] <- if (denominator > 0)
          (totals[, a] - estimate * totals[, b]) / denominator else rep(NA_real_, n_blocks)
      }
      summary_rows[[length(summary_rows) + 1L]] <- row
    }
    # Analytic paired ratio differences: keep repeated uses of each seed
    # together. No random resampling is used for these separate comparisons.
    pairs <- if (length(methods) > 1) combn(methods, 2L, simplify = FALSE) else list()
    for (pair in pairs) for (metric in names(numerators)) {
      influence <- influences[[pair[2]]][[metric]] - influences[[pair[1]]][[metric]]
      se <- if (n_blocks > 1 && all(is.finite(influence)))
        sqrt(n_blocks / (n_blocks - 1) * sum(influence^2)) else NA_real_
      delta <- unname(estimates[[pair[2]]][metric] - estimates[[pair[1]]][metric])
      difference_rows[[length(difference_rows) + 1L]] <- data.frame(
        cell[1, cell_vars, drop = FALSE], metric = metric,
        reference = pair[1], method = pair[2], difference = delta,
        standard_error = se, lower = delta - z * se, upper = delta + z * se,
        n_seed_blocks = n_blocks, interval_level = level,
        interval_method = "paired_seed_normal", row.names = NULL)
    }
  }
  list(summary = do.call(rbind, summary_rows), differences = do.call(rbind, difference_rows))
}
