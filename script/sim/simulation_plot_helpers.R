# Exact bin totals and seed-level uncertainty for SNP PIP calibration.
calibration_breaks <- (0:10) / 10

calibration_bin_totals <- function(pip, truth) {
  if (!length(pip) || any(!is.finite(pip)) || any(pip < -1e-10 | pip > 1 + 1e-10) ||
      !length(truth) || anyNA(truth) || any(truth != as.integer(truth)) ||
      any(truth < 1 | truth > length(pip))) stop("Invalid calibration PIPs or causal indices.")
  pip <- pmin(1, pmax(0, pip))
  # [0,.1), ... , [.9,1]; exact zero and one are retained.
  bin <- findInterval(pip, calibration_breaks, rightmost.closed = TRUE)
  sums <- numeric(10)
  z <- rowsum(pip, bin, reorder = FALSE)
  sums[as.integer(rownames(z))] <- z[, 1]
  c(tabulate(bin, 10), tabulate(bin[unique(truth)], 10), sums, 1)
}

calibration_merge <- function(state, meta, seeds, values) {
  if (!length(seeds)) return(invisible(NULL))
  key <- paste(meta$scenario, meta$pve, meta$K, meta$method, sep = "|")
  block <- rowsum(values, as.character(seeds), reorder = FALSE)
  old <- state[[key]]
  if (is.null(old)) {
    state[[key]] <- list(meta = meta, totals = block)
  } else {
    ids <- union(rownames(old$totals), rownames(block))
    combined <- matrix(0, length(ids), 31, dimnames = list(ids, NULL))
    combined[rownames(old$totals), ] <- old$totals
    combined[rownames(block), ] <- combined[rownames(block), , drop = FALSE] + block
    state[[key]] <- list(meta = meta, totals = combined)
  }
  invisible(NULL)
}

calibration_summary <- function(groups, pool_K = FALSE) {
  state <- new.env(hash = TRUE, parent = emptyenv())
  for (g in groups) {
    meta <- g$meta
    if (pool_K) meta$K <- 0L
    calibration_merge(state, meta, rownames(g$totals), g$totals)
  }
  rows <- lapply(as.list(state), function(g) {
    d <- g$totals
    count <- colSums(d[, 1:10, drop = FALSE])
    causal <- colSums(d[, 11:20, drop = FALSE])
    pip_sum <- colSums(d[, 21:30, drop = FALSE])
    frequency <- ifelse(count > 0, causal / count, NA_real_)
    n_seeds <- nrow(d)
    # Delta-method empirical SE of a pooled ratio, clustered by simulation
    # seed. Empty-bin seeds contribute zero; SNPs in LD are not independent.
    residual <- d[, 11:20, drop = FALSE] - sweep(d[, 1:10, drop = FALSE], 2, frequency, "*")
    se <- if (n_seeds > 1) sqrt(n_seeds / (n_seeds - 1) * colSums(residual^2)) / count else rep(NA_real_, 10)
    data.frame(g$meta, bin = 1:10, lower_pip = head(calibration_breaks, -1),
      upper_pip = tail(calibration_breaks, -1), n_snps = count, n_causal = causal,
      pip_sum = pip_sum, mean_pip = ifelse(count > 0, pip_sum / count, NA_real_),
      frequency = frequency, empirical_se = se,
      lower = pmax(0, frequency - 2 * se), upper = pmin(1, frequency + 2 * se),
      n_seed_blocks = n_seeds, n_nonempty_seed_blocks = colSums(d[, 1:10, drop = FALSE] > 0),
      n_replicates = sum(d[, 31]), row.names = NULL)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$scenario, out$pve, out$K, out$method, out$bin), ]
}

calibration_selection <- function(replicates, methods) {
  fields <- c("scenario", "pve", "K", "configuration", "seed", "method")
  d <- replicates[replicates$method %in% methods, fields]
  d <- d[do.call(order, unname(d)), ]
  rownames(d) <- NULL
  d
}

calibration_source_signature <- function(files) {
  info <- file.info(files)
  data.frame(file = basename(files), bytes = info$size, modified = as.numeric(info$mtime))
}

# Used when refreshing figures from existing summaries that predate calibration.
# Only the exact replicate/method selection used for those summaries is read.
collect_calibration <- function(files, replicates, methods, cache_file = NULL) {
  selection <- calibration_selection(replicates, methods)
  signature <- calibration_source_signature(files)
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$selection, selection) && identical(cached$signature, signature)) {
      cat("Reusing exact PIP calibration bin totals.\n")
      return(cached$groups)
    }
  }
  wanted <- selection[selection$method == methods[1], ]
  key <- function(d) paste(d$configuration, d$pve, d$seed, sep = "|")
  wanted_keys <- key(wanted)
  if (anyDuplicated(wanted_keys)) stop("Duplicate calibration replicate identities.")
  for (method in methods) if (!setequal(key(selection[selection$method == method, ]), wanted_keys))
    stop("Calibration methods must use identical replicates.")
  seen <- logical(nrow(wanted))
  state <- new.env(hash = TRUE, parent = emptyenv())
  for (f in seq_along(files)) {
    cat("Calibration: reading", f, "of", length(files), ":", basename(files[f]), "\n")
    info <- sim_parse_checkpoints(files[f])
    configuration <- sim_configuration(as.numeric(info[sim_count_columns]))
    saved <- new.env(parent = emptyenv())
    load(files[f], envir = saved)
    results <- saved$results
    seeds <- vapply(results, function(x) as.numeric(x$seed), 0)
    lookup <- match(paste(configuration, info$pve, seeds, sep = "|"), wanted_keys)
    eligible <- which(!is.na(lookup) & vapply(results, function(x) is.null(x$error), TRUE))
    eligible <- eligible[!seen[lookup[eligible]]]
    if (!length(eligible)) next
    meta <- wanted[lookup[eligible[1]], c("scenario", "pve", "K"), drop = FALSE]
    for (method in methods) {
      field <- c(SuSiE = "susie_pip", `SuSiE-mix` = "susie_mix_pip_snp",
                 `SuSiE-slide` = "susie_slide_pip")[[method]]
      values <- t(vapply(results[eligible], function(x) calibration_bin_totals(x[[field]], x$true_pos), numeric(31)))
      calibration_merge(state, cbind(meta, method = method), seeds[eligible], values)
    }
    seen[lookup[eligible]] <- TRUE
    rm(saved, results)
  }
  if (!all(seen)) stop("Current checkpoints are missing ", sum(!seen), " replicates used by the saved summaries; rerun plot_simulations.R.")
  groups <- as.list(state)
  if (!is.null(cache_file)) saveRDS(list(selection = selection, signature = signature, groups = groups), cache_file)
  groups
}

# Clip actual line segments before drawing, retaining their threshold order.
# Explicit clipping also avoids label loss on Windows PNG devices when curves
# extend beyond the displayed window.
visible_curve_segments <- function(x, y, xlim) {
  empty <- matrix(numeric(), ncol = 4,
                  dimnames = list(NULL, c("x0", "y0", "x1", "y1")))
  if (length(x) < 2) return(empty)
  a <- head(x, -1); b <- tail(x, -1)
  u <- head(y, -1); v <- tail(y, -1)
  ok <- is.finite(a) & is.finite(b) & is.finite(u) & is.finite(v) &
    pmax(a, b) >= xlim[1] & pmin(a, b) <= xlim[2]
  a <- a[ok]; b <- b[ok]; u <- u[ok]; v <- v[ok]
  if (!length(a)) return(empty)
  start <- rep(0, length(a)); end <- rep(1, length(a))
  moving <- a != b
  t0 <- (xlim[1] - a[moving]) / (b[moving] - a[moving])
  t1 <- (xlim[2] - a[moving]) / (b[moving] - a[moving])
  start[moving] <- pmax(0, pmin(t0, t1))
  end[moving] <- pmin(1, pmax(t0, t1))
  cbind(x0 = pmax(xlim[1], pmin(xlim[2], a + start * (b - a))),
        y0 = u + start * (v - u),
        x1 = pmax(xlim[1], pmin(xlim[2], a + end * (b - a))),
        y1 = u + end * (v - u))
}

# Maximal visible line height, including intersections with the x-window edges.
visible_curve_max <- function(d, xlim) {
  values <- numeric()
  # Keep distinct panel/method curves separate when scanning a whole figure.
  # Joining them would invent line segments and incorrect boundary heights.
  keys <- intersect(c("scenario", "pve", "K", "method"), names(d))
  for (dm in split(d, d[keys], drop = TRUE)) {
    dm <- dm[order(dm$threshold, decreasing = TRUE), ]
    x <- dm$fdr; y <- dm$tpr
    ok <- is.finite(x) & is.finite(y)
    values <- c(values, y[ok & x >= xlim[1] & x <= xlim[2]])
    if (length(x) > 1) for (edge in xlim) {
      a <- head(x, -1); b <- tail(x, -1)
      crossing <- which(head(ok, -1) & tail(ok, -1) & a != b &
                          pmin(a, b) <= edge & pmax(a, b) >= edge)
      values <- c(values, y[crossing] + (edge - a[crossing]) /
                    (b[crossing] - a[crossing]) * (y[crossing + 1] - y[crossing]))
    }
  }
  if (length(values)) max(values) else 0
}

figure_y_limits <- function(metric, d, xlim = c(0, .25)) {
  if (metric == "coverage") return(c(.7, 1))
  values <- if (metric %in% c("coverage", "purity", "cs_size")) d[[metric]] else numeric()
  values <- values[is.finite(values)]
  if (metric == "purity") {
    lower <- if (length(values)) min(values) else 0
    # An all-one figure needs a nonzero range to draw axes.
    return(c(if (lower >= 1) .95 else lower, 1))
  }
  if (metric %in% c("cs_size", "power_fdr")) {
    upper <- if (metric == "power_fdr") visible_curve_max(d, xlim) else
      if (length(values)) max(values) else 0
    return(c(0, if (upper > 0) upper else 1))
  }
  c(0, 1)
}

panel_y_ticks <- function(limits) {
  ticks <- pretty(limits, n = 4)
  ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
  # Include the actual endpoint values, not a rounded expansion of the limits.
  ticks <- ticks[ticks > limits[1] + diff(limits) * .12 &
                   ticks < limits[2] - diff(limits) * .12]
  sort(unique(c(limits, ticks)))
}
