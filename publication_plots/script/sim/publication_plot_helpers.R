# Publication labels are separate from the saved method identifiers.
# In particular, SuSiE-slide below ALWAYS means the finite-prior method.
publication_methods <- c("SuSiE", "SuSiE-mix", "SuSiE-slide-prior")
publication_labels <- setNames(c("SuSiE", "SuSiE-MiX", "SuSiE-slide"), publication_methods)
publication_colors <- setNames(c("#2489FF", "#D81B60", "#009E73"), publication_methods)
publication_shapes <- setNames(c(16, 17, 18), publication_methods)
publication_lines <- setNames(c(1, 2, 4), publication_methods)

# Requires simulation_plot_helpers.R (for seed-block calibration intervals).
# Filter BEFORE pooling: an all-K summary cannot be trimmed after pooling.
read_publication_simulation_data <- function(input_dir, max_causal = 4L,
                                             methods = publication_methods) {
  stopifnot(length(max_causal) == 1L, is.finite(max_causal),
            max_causal >= 1L, max_causal == as.integer(max_causal),
            all(methods %in% publication_methods), !anyDuplicated(methods))
  paths <- file.path(input_dir, c("metric_summary.csv", "roc_counts.rds",
                                 "pip_calibration_seed_counts.rds"))
  if (any(!file.exists(paths))) stop("Missing saved summaries:\n",
    paste(paths[!file.exists(paths)], collapse = "\n"),
    "\nRun the full simulation summary script first.")
  keep <- function(d) d$method %in% methods & d$K %in% seq_len(max_causal)
  metrics <- read.csv(paths[1], stringsAsFactors = FALSE)
  metrics <- metrics[keep(metrics), , drop = FALSE]
  if (!nrow(metrics) || !setequal(unique(metrics$method), methods))
    stop("No complete set of requested methods in the selected causal counts.")
  cells <- split(metrics, interaction(metrics[c("scenario", "pve", "K")], drop = TRUE))
  for (d in cells) {
    if (nrow(d) != length(methods) || !setequal(d$method, methods) ||
        length(unique(d$n_replicates)) != 1L || length(unique(d$n_causal)) != 1L)
      stop("Selected methods must have matching simulation counts in every cell.")
  }
  curves <- readRDS(paths[2])
  curves <- curves[keep(curves), , drop = FALSE]
  curves$fdr <- with(curves, ifelse(tp + fp > 0, fp / (tp + fp), 0))
  baseline <- curves[curves$threshold == 0, ]
  key <- function(d) paste(d$scenario, d$pve, d$K, d$method, sep = "|")
  mt <- match(key(metrics), key(baseline))
  if (anyNA(mt) || anyDuplicated(key(baseline)) ||
      !setequal(key(metrics), key(baseline)) ||
      any(baseline$tp[mt] != metrics$n_causal))
    stop("ROC counts and metric summaries disagree; refresh the saved summaries.")

  pooled_curves <- aggregate(curves[c("tp", "fp")],
    curves[c("scenario", "pve", "method", "threshold")], sum)
  group_key <- function(d) paste(d$scenario, d$pve, d$method, sep = "|")
  pooled_baseline <- pooled_curves[pooled_curves$threshold == 0, ]
  mt <- match(group_key(pooled_curves), group_key(pooled_baseline))
  pooled_curves$tpr <- pooled_curves$tp / pooled_baseline$tp[mt]
  pooled_curves$fpr <- pooled_curves$fp / pooled_baseline$fp[mt]
  pooled_curves$fdr <- with(pooled_curves, ifelse(tp + fp > 0, fp / (tp + fp), 0))

  cache <- readRDS(paths[3])
  groups <- Filter(function(g) isTRUE(keep(g$meta)), cache$groups)
  calibration <- calibration_summary(groups)
  totals <- aggregate(calibration$n_causal, calibration[c("scenario", "pve", "K", "method")], sum)
  mt <- match(key(metrics), key(totals))
  rep_counts <- unique(calibration[c("scenario", "pve", "K", "method", "n_replicates")])
  mr <- match(key(metrics), key(rep_counts))
  if (anyNA(mt) || anyNA(mr) || !setequal(key(metrics), key(totals)) ||
      any(totals$x[mt] != metrics$n_causal) ||
      any(rep_counts$n_replicates[mr] != metrics$n_replicates))
    stop("Calibration cache and metric summaries disagree; refresh the saved summaries.")
  pooled_calibration <- calibration_summary(groups, pool_K = TRUE)
  list(metrics = metrics, curves = curves, pooled_curves = pooled_curves,
       calibration = calibration, pooled_calibration = pooled_calibration,
       paths = paths, max_causal = max_causal)
}
