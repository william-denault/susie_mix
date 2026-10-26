# Small synthetic summaries check that K=5 contributes to every pooled panel.
local({
  source("publication_plots/script/sim/simulation_plot_helpers.R", local = TRUE)
  source("publication_plots/script/sim/publication_plot_helpers.R", local = TRUE)
  stopifnot(identical(figure_y_limits("coverage", data.frame(coverage = c(.64, .95))), c(.6, 1)),
    identical(figure_y_limits("coverage", data.frame(coverage = c(.8, .95))), c(.7, 1)))
  input <- tempfile("publication-counts-")
  dir.create(input)
  metrics <- expand.grid(scenario = c("Additive only", "Recessive only", "Partial recessive only"),
    pve = publication_pve_values, K = 1:5, method = publication_methods,
    stringsAsFactors = FALSE)
  metrics$n_replicates <- 2L
  metrics$n_causal <- 2L * metrics$K
  metrics$purity <- .9
  metrics$coverage <- metrics$power <- metrics$cs_size <- 1
  metrics$cs_size_lo <- metrics$cs_size_hi <- 1
  curves <- lapply(seq_len(nrow(metrics)), function(i) {
    d <- metrics[i, c("scenario", "pve", "K", "method")]
    transform(d[rep(1, 4), ], threshold = c(0, .8, .95, Inf),
      tp = c(2*d$K, 2*d$K, 0, 0), fp = c(2*(10-d$K), 0, 0, 0),
      tpr = c(1, 1, 0, 0), fpr = c(1, 0, 0, 0))
  })
  groups <- lapply(seq_len(nrow(metrics)), function(i) {
    d <- metrics[i, c("scenario", "pve", "K", "method")]
    v <- calibration_bin_totals(c(rep(.9, d$K), rep(.1, 10-d$K)), seq_len(d$K))
    values <- rbind(v, v)
    rownames(values) <- c("1", "2")
    list(meta = d, totals = values)
  })
  cache <- list(signature = data.frame(file = "fixture_n250_L10_results.RData"), groups = groups)
  write.csv(metrics, file.path(input, "metric_summary.csv"), row.names = FALSE)
  saveRDS(do.call(rbind, curves), file.path(input, "roc_counts.rds"))
  saveRDS(cache, file.path(input, "pip_calibration_seed_counts.rds"))
  five <- read_publication_simulation_data(input)
  four <- read_publication_simulation_data(input, max_causal = 4L)
  stopifnot(setequal(five$metrics$K, 1:5), all(five$pooled_curves$tp[five$pooled_curves$threshold == 0] == 30L),
    all(four$pooled_curves$tp[four$pooled_curves$threshold == 0] == 20L),
    all(five$pooled_calibration$n_replicates == 10L),
    all(four$pooled_calibration$n_replicates == 8L))
  totals <- aggregate(five$pooled_calibration$n_causal,
    five$pooled_calibration[c("scenario", "pve", "method")], sum)
  stopifnot(all(totals$x == 30L))
  write.csv(metrics[metrics$K < 5, ], file.path(input, "metric_summary.csv"), row.names = FALSE)
  stopifnot(inherits(try(read_publication_simulation_data(input), silent = TRUE), "try-error"))
  write.csv(metrics, file.path(input, "metric_summary.csv"), row.names = FALSE)
  cache$signature$file <- "fixture_n500_L10_results.RData"
  saveRDS(cache, file.path(input, "pip_calibration_seed_counts.rds"))
  stopifnot(inherits(try(read_publication_simulation_data(input), silent = TRUE), "try-error"))
  cat("PASS: K=5 included in metrics, pooled ROC counts and calibration; incomplete K range and old n=500 inputs rejected.\n")
})
