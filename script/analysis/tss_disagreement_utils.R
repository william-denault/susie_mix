# Select individual CSs whose biological lead SNP is absent from every CS of
# the other model in the SAME gene/tissue. Coding differences alone are not
# TSS disagreement. Shared leads are removed even in partially shared regions.
select_tss_disagreement <- function(cs, additive_model = "SuSiE", mixed_model = "SuSiE-mix",
                                    plot_limit_kb = 200) {
  required <- c("gene", "tissue", "model", "lead_snp", "distance_to_tss_kb")
  if (!all(required %in% names(cs))) stop("CS summary is missing TSS agreement fields: ",
                                         paste(setdiff(required, names(cs)), collapse = ", "))
  x <- as.data.frame(cs)
  models <- c(additive_model, mixed_model)
  if (length(models) != 2L || anyNA(models) || anyDuplicated(models)) stop("Specify two distinct models.")
  x <- x[x$model %in% models, , drop = FALSE]
  if (anyNA(x[c("gene", "tissue")]) || any(!nzchar(x$gene) | !nzchar(x$tissue))) {
    stop("TSS comparisons need nonempty gene and tissue names.")
  }
  part <- function(v) paste0(nchar(v), ":", v)
  region <- paste0(part(x$gene), part(x$tissue))
  lead <- trimws(as.character(x$lead_snp))
  valid_lead <- !is.na(lead) & nzchar(lead)
  lead_key <- paste0(region, part(lead))
  x$tss_other_model_cs_count <- integer(nrow(x))
  x$tss_shared_lead <- rep(FALSE, nrow(x))
  x$tss_agreement_status <- rep("unknown_lead", nrow(x))
  for (model_name in models) {
    own <- which(x$model == model_name)
    other <- which(x$model != model_name)
    counts <- table(region[other])
    other_n <- as.integer(counts[region[own]])
    other_n[is.na(other_n)] <- 0L
    shared <- valid_lead[own] & lead_key[own] %in% lead_key[other[valid_lead[other]]]
    unknown_other <- region[own] %in% region[other[!valid_lead[other]]]
    status <- ifelse(shared, "agreed_lead", ifelse(other_n == 0L, "other_model_no_cs",
              ifelse(unknown_other, "unknown_other_model_lead", "different_lead")))
    status[!valid_lead[own]] <- "unknown_lead"
    x$tss_other_model_cs_count[own] <- other_n
    x$tss_shared_lead[own] <- shared
    x$tss_agreement_status[own] <- status
  }
  x$tss_include <- x$tss_agreement_status %in% c("different_lead", "other_model_no_cs")
  # This flag describes the displayed window; KDE fitting uses all finite distances.
  x$tss_finite_distance <- is.finite(x$distance_to_tss_kb)
  x$tss_in_plot_window <- x$tss_finite_distance &
    abs(x$distance_to_tss_kb) <= plot_limit_kb
  x$tss_agreement_rule <- rep("same biological lead SNP within gene/tissue; coding ignored", nrow(x))
  summaries <- list()
  for (tissue_name in c("All tissues", sort(unique(as.character(x$tissue))))) {
    for (model_name in models) {
      z <- x[x$model == model_name & (tissue_name == "All tissues" | x$tissue == tissue_name), , drop = FALSE]
      summaries[[length(summaries) + 1L]] <- data.frame(
        tissue = tissue_name, model = model_name, n_cs_input = nrow(z),
        n_cs_agreed_excluded = sum(z$tss_shared_lead),
        n_cs_different_lead = sum(z$tss_agreement_status == "different_lead"),
        n_cs_other_model_no_cs = sum(z$tss_agreement_status == "other_model_no_cs"),
        n_cs_unknown_excluded = sum(grepl("^unknown", z$tss_agreement_status)),
        n_cs_selected = sum(z$tss_include),
        n_selected_missing_distance = sum(z$tss_include & !z$tss_finite_distance),
        n_selected_outside_window = sum(z$tss_include & z$tss_finite_distance & !z$tss_in_plot_window),
        n_cs_in_window = sum(z$tss_include & z$tss_in_plot_window))
    }
  }
  audit_columns <- unique(c(intersect(c("gene", "tissue", "model", "model_key", "cs_id", "cs_number",
                                         "lead_snp", "lead_coding", "distance_to_tss_kb", "em_iteration"), names(x)),
                            grep("^tss_", names(x), value = TRUE)))
  list(selected = x[x$tss_include, , drop = FALSE], audit = x[audit_columns],
       summary = do.call(rbind, summaries))
}

# Gaussian KDE of raw signed distances. A fixed bandwidth gives both models
# and all tissues the same smoothing. The display limits never truncate the
# observations or renormalize the curve inside the displayed window.
summarize_tss_kde <- function(x, model_name, tissue_name = "All tissues",
                              plot_limit_kb = 200, bandwidth_kb = 10, n = 1025L) {
  if (length(bandwidth_kb) != 1L || !is.finite(bandwidth_kb) || bandwidth_kb <= 0)
    stop("TSS bandwidth_kb must be a finite positive number.")
  if (length(plot_limit_kb) != 1L || !is.finite(plot_limit_kb) || plot_limit_kb <= 0)
    stop("TSS plot_limit_kb must be a finite positive number.")
  if (length(n) != 1L || !is.finite(n) || n < 3 || n != as.integer(n))
    stop("TSS KDE grid size must be an integer of at least 3.")
  if (!"distance_to_tss_kb" %in% names(x)) stop("Missing distance_to_tss_kb.")
  d <- x$distance_to_tss_kb[is.finite(x$distance_to_tss_kb)]
  grid <- seq(-plot_limit_kb, plot_limit_kb, length.out = n)
  estimated <- rep(NA_real_, n)
  if (length(d) == 1L) {
    estimated <- dnorm(grid, mean = d, sd = bandwidth_kb)
  } else if (length(d) > 1L) {
    kde <- density(d, bw = bandwidth_kb, kernel = "gaussian", n = n,
                   from = -plot_limit_kb, to = plot_limit_kb)
    grid <- kde$x
    estimated <- kde$y
  }
  data.frame(
    tissue = tissue_name, model = model_name,
    analysis_set = "CSs without a shared biological lead SNP",
    distance_to_tss_kb = grid, density = estimated,
    bandwidth_kb = bandwidth_kb, kernel = "gaussian",
    normalization = "all finite selected distances; unit area over full domain",
    n_cs_selected = nrow(x), n_cs_total = length(d),
    n_cs_in_window = sum(abs(d) <= plot_limit_kb),
    n_cs_missing_distance = nrow(x) - length(d),
    n_cs_outside_window = sum(abs(d) > plot_limit_kb),
    plot_limit_kb = plot_limit_kb,
    stringsAsFactors = FALSE
  )
}

plot_tss_kde <- function(x, model_colors, plot_limit_kb = 200,
                         main_title = "TSS distance for disagreeing CS leads",
                         compact = FALSE, set_margins = FALSE) {
  x <- as.data.frame(x)
  valid <- is.finite(x$density)
  ymax <- if (any(valid)) max(1e-6, 1.12 * max(x$density[valid])) else 1
  if (set_margins) par(
    mar = if (compact) c(3.8, 4.2, 3.2, 1) else c(5, 5, 4, 1),
    mgp = if (compact) c(2.4, .7, 0) else c(2.8, .8, 0)
  )
  plot(NA_real_, NA_real_, xlim = c(-plot_limit_kb, plot_limit_kb), ylim = c(0, ymax),
       xaxs = "i", yaxs = "i", xaxt = "n", bty = "l",
       xlab = "Distance to TSS (kb)", ylab = "Density (per kb)",
       main = main_title, cex.main = if (compact) .85 else 1,
       cex.lab = if (compact) .8 else 1, cex.axis = if (compact) .75 else 1)
  axis(1, at = seq(-plot_limit_kb, plot_limit_kb, length.out = 5),
       cex.axis = if (compact) .75 else 1)
  abline(v = 0, col = "gray65", lty = 3)
  for (model_name in names(model_colors)) {
    z <- x[x$model == model_name, , drop = FALSE]
    if (any(is.finite(z$density))) lines(z$distance_to_tss_kb, z$density,
                                        col = model_colors[[model_name]], lwd = 2.5)
  }
  if (!any(valid)) text(0, ymax / 2, "No disagreeing CS leads with finite TSS distance",
                        cex = if (compact) .7 else .9)
  legend("topright", legend = vapply(names(model_colors), function(model_name) {
    total <- unique(x$n_cs_total[x$model == model_name])
    sprintf("%s (n = %d)", model_name, if (length(total)) total[1] else 0L)
  }, character(1)), col = unname(model_colors), lwd = 2.5, bty = "n",
    cex = if (compact) .65 else .8)
  bandwidth <- unique(x$bandwidth_kb)
  if (length(bandwidth) == 1L) mtext(
    sprintf("Gaussian KDE; bandwidth = %g kb; n = all finite selected CSs", bandwidth),
    side = 3, line = .15, cex = if (compact) .55 else .65)
  invisible(NULL)
}
