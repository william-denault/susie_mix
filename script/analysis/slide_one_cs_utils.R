# Source slide_analysis_utils.R before this file.
slide_one_cs_comparisons <- function(res, results_dir, expected_coding = "dominant",
    association_threshold = 5e-8, minimum_mean_reads = 100, both_one_cs = FALSE,
    coding_selection = c("endpoints", "direction"), delta_tolerance = 1e-6) {
  expected_coding <- match.arg(expected_coding, c("dominant", "recessive"))
  coding_selection <- match.arg(coding_selection)
  cases <- slide_primary(res, association_threshold, minimum_mean_reads)
  count <- cases[[paste0("n_", expected_coding, "_slide")]]
  if (coding_selection == "direction") count <- count + cases[[paste0("n_partial_", expected_coding, "_slide")]]
  cases <- cases[which(cases$ncs_slide == 1L & count == 1L & (!both_one_cs | cases$ncs_susie == 1L)), , drop = FALSE]
  cases <- cases[order(cases$gene, cases$tissue), , drop = FALSE]
  empty <- data.frame(gene = character(), tissue = character(), status = character(), message = character(),
    additive_lead_snp = character(), lead_snp = character(), additive_n_cs = integer(), n_cs = integer(),
    lead_coding = character(), slide_coding_class = character(), lead_delta = numeric(),
    lead_pip = numeric(), distance_bp = numeric(), distance_kb = numeric())
  cached_gene <- NULL
  rows <- list()
  for (i in seq_len(nrow(cases))) {
    gene <- cases$gene[i]
    tissue <- cases$tissue[i]
    if (!identical(cached_gene, gene)) {
      out <- tryCatch(readRDS(file.path(results_dir, cases$result_file[i])), error = identity)
      cached_gene <- gene
    }
    rows[[i]] <- tryCatch({
      if (inherits(out, "error")) stop(conditionMessage(out))
      x <- out[[tissue]]
      map <- x$add_predictor_map
      slide_validate_fit(x[["susie_add"]], map)
      slide_validate_fit(x[["fit_slide"]], map, TRUE)
      add_fit <- x$susie_add
      fit <- x$fit_slide
      add_leads <- slide_leads(add_fit, map)
      lead <- slide_leads(fit, map, TRUE, delta_tolerance)
      allowed <- c(expected_coding, if (coding_selection == "direction") paste0("partial_", expected_coding))
      if (nrow(lead) != 1L || !lead$lead_coding %in% allowed || nrow(add_leads) != cases$ncs_susie[i])
        stop("Saved fits differ from selected summary; regenerate the slide summary.")
      if (nrow(add_leads)) {
        selected <- which.max(add_leads$lead_pip)
        add_index <- add_leads$index[selected]
        add_snps <- map$snp[get_cs(add_fit)[[selected]]]
        selection <- if (nrow(add_leads) == 1L) "Single credible-set lead" else "Highest-PIP credible-set lead"
      } else {
        add_index <- which.max(add_fit$pip)
        add_snps <- character()
        selection <- "Highest-PIP SNP (no credible set)"
      }
      slide_snps <- map$snp[get_cs(fit)[[1]]]
      coordinates <- get_predictor_coordinates(map)
      indices <- c(add_index, lead$index)
      pos <- coordinates$position[indices]
      chr <- tolower(sub("^chr", "", coordinates$chromosome[indices], ignore.case = TRUE))
      valid <- all(is.finite(pos) & pos > 0) && !anyNA(chr) && chr[1] == chr[2]
      distance <- if (valid) abs(diff(pos)) else NA_real_
      shared <- length(intersect(add_snps, slide_snps))
      data.frame(gene = gene, tissue = tissue, status = "compared", message = "Saved-fit comparison",
        min_pv = cases$min_pv[i], mean_count = cases$mean_count[i], additive_n_cs = nrow(add_leads), n_cs = 1L,
        additive_lead_snp = map$snp[add_index], lead_snp = lead$lead_snp,
        lead_coding = expected_coding, slide_coding_class = lead$lead_coding,
        lead_delta = lead$lead_delta, lead_delta_forced = lead$lead_delta_forced,
        lead_component = lead$component, additive_lead_pip = add_fit$pip[add_index], lead_pip = lead$lead_pip,
        additive_lead_selection = selection, same_lead_snp = map$snp[add_index] == lead$lead_snp,
        distance_bp = distance, distance_kb = distance / 1000, valid_coordinates = valid,
        selected_cs_overlap_snps = shared,
        selected_cs_jaccard = if (length(add_snps)) shared / length(union(add_snps, slide_snps)) else NA_real_,
        dif_elbo_slide_vs_add = max_elbo(fit) - max_elbo(add_fit),
        coding_selection = coding_selection)
    }, error = function(e) data.frame(gene = gene, tissue = tissue, status = "error", message = conditionMessage(e)))
  }
  result <- slide_bind_rows(rows, empty)
  for (nm in setdiff(names(empty), names(result))) result[[nm]] <- rep(NA, nrow(result))
  result[order(-result$distance_bp, result$gene, result$tissue, na.last = TRUE), , drop = FALSE]
}

plot_one_cs_slide <- function(project_dir, expected_coding = "dominant",
    results_dir = file.path(project_dir, "results"), summary_dir = file.path(project_dir, "results_slide/summary"),
    output_dir = file.path(project_dir, "results_slide/plot", paste0("one_cs_", expected_coding)),
    association_threshold = 5e-8, minimum_mean_reads = 100, both_one_cs = FALSE,
    coding_selection = "endpoints", top_n = 20L) {
  summaries <- slide_load_summary(summary_dir)
  tolerance <- unique(summaries$res$delta_tolerance)
  if (length(tolerance) != 1L) stop("Summary contains multiple delta tolerances.")
  comparisons <- slide_one_cs_comparisons(summaries$res, results_dir, expected_coding,
    association_threshold, minimum_mean_reads, both_one_cs, coding_selection, tolerance)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(comparisons, file.path(output_dir, "lead_snp_comparisons.csv"), row.names = FALSE)
  completed <- comparisons[comparisons$status == "compared", , drop = FALSE]
  plot_input <- completed
  plot_input$status <- rep("plotted", nrow(plot_input))
  overview <- plot_one_cs_lead_distance_overview(plot_input, expected_coding,
    file.path(output_dir, "lead_snp_distance_overview.png"), top_n = top_n,
    title_prefix = paste0("SuSiE-slide (", coding_selection, ") | "), annotation_mode = "largest_shifts",
    lead_pip_label = "SuSiE-slide lead SNP PIP",
    pip_note = "Distance uses base-pair coordinates. Slide PIP is SNP-level; coding is the fitted lead delta.")
  write.csv(overview$top_cases, file.path(output_dir, "largest_lead_shifts.csv"), row.names = FALSE)
  distances <- completed$distance_kb[is.finite(completed$distance_kb)]
  stats <- data.frame(coding = expected_coding, coding_selection = coding_selection, delta_tolerance = tolerance,
    association_threshold = association_threshold, minimum_mean_reads = minimum_mean_reads,
    both_one_cs = both_one_cs, n_candidates = nrow(comparisons), n_compared = nrow(completed),
    n_errors = sum(comparisons$status == "error"), n_valid_distances = length(distances),
    n_changed_leads = overview$n_changed, n_over_100kb = sum(distances > 100),
    median_distance_kb = if (length(distances)) median(distances) else NA_real_,
    max_distance_kb = if (length(distances)) max(distances) else NA_real_)
  write.csv(stats, file.path(output_dir, "lead_distance_summary.csv"), row.names = FALSE)
  print(stats)
  invisible(list(comparisons = comparisons, summary = stats, overview = overview))
}
