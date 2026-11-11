# Plot saved priors without changing any history or fitted results.
# Accepts a CSV filename or an already-loaded data.frame/tibble.
plot_slider_prior_history <- function(history, tissues = NULL, ncol = 6L,
                                      batches_per_iteration = 9L) {
  if (is.character(history) && length(history) == 1L) {
    history <- read.csv(history, check.names = FALSE, stringsAsFactors = FALSE)
  }
  h <- as.data.frame(history)
  if (!"batches_in_iteration" %in% names(h)) h$batches_in_iteration <- batches_per_iteration
  weights <- sprintf("w_%02d", seq_len(17L))
  required <- c("iteration", "tissue", "phase", "batch", "prior_role", "stage", weights)
  missing <- setdiff(required, names(h))
  if (length(missing)) stop("Missing history columns: ", paste(missing, collapse = ", "))

  # Retain the original full-EM history, the initial online prior, and all
  # completed updates. Later input snapshots repeat the preceding update.
  keep <- is.na(h$stage) | h$stage == 0L | h$prior_role == "after_update"
  h <- h[keep, , drop = FALSE]
  if (is.null(tissues)) tissues <- unique(h$tissue)
  if (!length(tissues) || any(!tissues %in% h$tissue)) stop("No data for requested tissues")
  h <- h[h$tissue %in% tissues, , drop = FALSE]
  h <- h[order(h$iteration, h$batch, h$tissue), , drop = FALSE]
  if (any(!is.finite(as.matrix(h[, weights])))) stop("Prior weights contain missing/nonfinite values")

  grid_columns <- sprintf("slider_%02d", seq_len(17L))
  slider <- if (all(grid_columns %in% names(h))) {
    as.numeric(unlist(h[1L, grid_columns], use.names = FALSE))
  } else seq(-1, 1, length.out = 17L)

  # Each saved update has the same colour across tissues.
  key <- paste(h$iteration, h$batch, sep = "/")
  update_keys <- unique(key)
  colours <- setNames(grDevices::hcl.colors(length(update_keys), "Dark 3"), update_keys)
  online <- h$phase == "minibatch" & h$prior_role == "after_update"
  completed <- sort(unique(h$iteration[online & h$batch == h$batches_in_iteration]))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  ncol <- min(as.integer(ncol), length(tissues))
  par(mfrow = c(ceiling(length(tissues) / ncol), ncol),
      mar = c(3.5, 4.0, 2.2, 0.8), oma = c(3.5, 0, 0, 0),
      mgp = c(2.3, 0.65, 0))

  for (tissue in tissues) {
    d <- h[h$tissue == tissue, , drop = FALSE]
    W <- as.matrix(d[, weights])
    is_online <- d$phase == "minibatch" & d$prior_role == "after_update"
    is_end <- is_online & d$batch == d$batches_in_iteration
    line_colours <- colours[paste(d$iteration, d$batch, sep = "/")]

    plot(slider, W[1L, ], type = "n", ylim = range(W), main = tissue,
         xlab = "Slider value", ylab = "Slider prior probability",
         cex.main = 1, cex.lab = 0.85, cex.axis = 0.8)
    abline(v = 0, col = "grey80", lwd = 0.8)
    for (j in seq_len(nrow(d))) {
      lines(slider, W[j, ], col = line_colours[j],
            lty = if (is_online[j]) 2 else 1, lwd = 1.2)
    }
    # Draw completed online iterations last so their endpoints stay visible.
    # The latest batch of an unfinished iteration is never marked as an end.
    for (j in which(is_end)) {
      lines(slider, W[j, ], col = "black", lty = 1, lwd = 3)
    }
  }
  mtext("Solid: full EM / initial prior     Dashed: online batch updates     Thick black: completed online pass",
        side = 1, outer = TRUE, line = 0.8, cex = 0.8)
  mtext(paste("Completed online iterations:",
              if (length(completed)) paste(sprintf("%03d", completed), collapse = ", ") else "none yet"),
        side = 1, outer = TRUE, line = 2.1, cex = 0.8)
  invisible(list(history = h, completed_online_iterations = completed))
}
