# Saved-fit comparisons require only base R. Expression plotting is optional.
em_one_cs_lead <- function(fit, predictor_map) {
  validate_predictor_map(fit, predictor_map, "lead predictor map")
  if (length(fit$pip) != ncol(fit$alpha) || !any(is.finite(fit$pip))) stop("No usable lead PIPs.")
  cs <- get_cs(fit)
  indices <- if (length(cs)) vapply(seq_along(cs), function(i) {
    get_cs_lead_predictor(fit, cs, i)
  }, integer(1)) else seq_along(fit$pip)
  scores <- fit$pip[indices]
  scores[!is.finite(scores)] <- -Inf
  which_lead <- which.max(scores)
  index <- indices[which_lead]
  list(snp = as.character(predictor_map$snp[index]),
       coding = as.character(predictor_map$coding[index]), pip = fit$pip[index],
       snps = if (length(cs)) unique(predictor_map$snp[cs[[which_lead]]]) else character(),
       selection = if (!length(cs)) "Highest-PIP SNP (no credible set)" else
         if (length(cs) == 1L) "Single credible-set lead" else "Highest-PIP credible-set lead")
}

em_one_cs_distances <- function(context, res, expected_coding = "dominant",
                                association_threshold = 5e-8, minimum_mean_reads = 100,
                                both_one_cs = FALSE) {
  expected_coding <- match.arg(expected_coding, c("dominant", "recessive"))
  count_column <- if (expected_coding == "dominant") "n_dom_weighted" else "n_rec_weighted"
  keep <- is.finite(res$min_pv) & res$min_pv < association_threshold &
    is.finite(res$mean_count) & res$mean_count >= minimum_mean_reads &
    res$ncs_weighted_fit_mix == 1L & res[[count_column]] == 1L
  if (both_one_cs) keep <- keep & res$ncs_susie == 1L
  cases <- res[which(keep), , drop = FALSE]
  cases <- cases[order(cases$gene, cases$tissue), , drop = FALSE]
  rows <- list()
  # Typed empty table also replaces stale output when no cases qualify.
  empty <- data.frame(gene = character(), tissue = character(), status = character(),
                      message = character(), em_iteration = integer(),
                      additive_lead_snp = character(), lead_snp = character(),
                      additive_n_cs = integer(), n_cs = integer(), lead_coding = character(),
                      lead_pip = numeric(), distance_bp = numeric(), distance_kb = numeric())
  cached_gene <- NULL
  for (i in seq_len(nrow(cases))) {
    gene <- cases$gene[i]
    tissue <- cases$tissue[i]
    if (!identical(cached_gene, gene)) {
      joined <- tryCatch(em_join_gene(file.path(context$iteration_dir, "results", paste0(gene, ".rds")),
                                     context), error = identity)
      cached_gene <- gene
    }
    rows[[i]] <- tryCatch({
      if (inherits(joined, "error")) stop(conditionMessage(joined))
      x <- joined[[tissue]]
      if (is.null(x) || length(get_cs(x$weighted_fit_mix)) != 1L) stop("Selected EM tissue/CS is unavailable.")
      add <- em_one_cs_lead(x$susie_add, x$add_predictor_map)
      mix <- em_one_cs_lead(x$weighted_fit_mix, x$mix_predictor_map)
      if (mix$coding != expected_coding) stop("Selected coding differs from saved EM fit.")
      chr <- parse_variant_chromosome(c(add$snp, mix$snp))
      pos <- parse_variant_position(c(add$snp, mix$snp))
      valid <- all(is.finite(pos) & pos > 0) && !anyNA(chr) && chr[1] == chr[2]
      distance <- if (valid) abs(diff(pos)) else NA_real_
      shared <- length(intersect(add$snps, mix$snps))
      data.frame(
        gene = gene, tissue = tissue, status = "compared", message = "Saved-fit comparison",
        em_iteration = context$iteration, min_pv = cases$min_pv[i], mean_count = cases$mean_count[i],
        additive_n_cs = length(get_cs(x$susie_add)), n_cs = 1L,
        additive_lead_snp = add$snp, lead_snp = mix$snp, lead_coding = mix$coding,
        additive_lead_pip = add$pip, lead_pip = mix$pip,
        additive_lead_selection = add$selection, same_lead_snp = add$snp == mix$snp,
        distance_bp = distance, distance_kb = distance / 1000,
        valid_coordinates = valid, selected_cs_overlap_snps = shared,
        selected_cs_jaccard = if (length(add$snps)) shared / length(union(add$snps, mix$snps)) else NA_real_,
        likelihood_statistic = 2 * (get_log_lik_metric(x$weighted_fit_mix) - get_log_lik_metric(x$susie_add)),
        stringsAsFactors = FALSE)
    }, error = function(e) data.frame(gene = gene, tissue = tissue, status = "error",
                                      message = conditionMessage(e), em_iteration = context$iteration))
  }
  if (!length(rows)) return(empty)
  columns <- unique(c(names(empty), unlist(lapply(rows, names))))
  rows <- lapply(rows, function(row) {
    for (column in setdiff(columns, names(row))) row[[column]] <- NA
    row[columns]
  })
  out <- do.call(rbind, rows)
  out[order(-out$distance_bp, out$gene, out$tissue, na.last = TRUE), , drop = FALSE]
}

em_plot_one_cs <- function(project_dir, iteration = "latest", expected_coding = "dominant",
                           association_threshold = 5e-8, minimum_mean_reads = 100,
                           both_one_cs = FALSE, expression_plots = FALSE,
                           datadir = "/project2/mstephens/gtex", top_n = 20L) {
  context <- em_analysis_context(project_dir, iteration)
  source(file.path(context$project_dir, "script/analysis/generate_summary_results.R"), local = TRUE)
  source(file.path(context$project_dir, "script/analysis/one_cs_lead_distance_plot_utils.R"), local = TRUE)
  # Bind helpers to the shared summary definitions without polluting the caller.
  em_join_gene <- em_join_gene
  environment(em_join_gene) <- environment()
  em_one_cs_lead <- em_one_cs_lead
  environment(em_one_cs_lead) <- environment()
  compare <- em_one_cs_distances
  environment(compare) <- environment()
  summaries <- em_load_summary(context)
  comparisons <- compare(context, summaries$res, expected_coding,
                         association_threshold, minimum_mean_reads, both_one_cs)
  output_dir <- file.path(context$iteration_dir, "plot", paste0("one_cs_", expected_coding))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(comparisons, file.path(output_dir, "lead_snp_comparisons.csv"), row.names = FALSE)
  completed <- comparisons[comparisons$status == "compared", , drop = FALSE]
  # The shared overview consumes the historical plotting status label. These
  # rows come from all saved fits, independent of expression-plot success.
  overview_input <- completed
  overview_input$status <- rep("plotted", nrow(overview_input))
  overview <- plot_one_cs_lead_distance_overview(
    overview_input, expected_coding, file.path(output_dir, "lead_snp_distance_overview.png"),
    top_n = top_n, title_prefix = paste0("EM ", context$iteration, " | "), annotate_examples = FALSE)
  write.csv(overview$top_cases, file.path(output_dir, "largest_lead_shifts.csv"), row.names = FALSE)
  distance <- completed$distance_bp[is.finite(completed$distance_bp)]
  stats <- data.frame(
    em_iteration = context$iteration, coding = expected_coding,
    association_threshold = association_threshold, minimum_mean_reads = minimum_mean_reads,
    both_one_cs = both_one_cs, n_candidates = nrow(comparisons), n_compared = nrow(completed),
    n_errors = sum(comparisons$status == "error"), n_valid_distances = length(distance),
    n_changed_leads = overview$n_changed, n_over_100kb = sum(distance > 1e5),
    median_distance_kb = if (length(distance)) median(distance) / 1000 else NA_real_,
    max_distance_kb = if (length(distance)) max(distance) / 1000 else NA_real_)
  write.csv(stats, file.path(output_dir, "lead_distance_summary.csv"), row.names = FALSE)
  print(stats)
  if (expression_plots) {
    source(file.path(context$project_dir, "script/analysis/fit_mix_expression_plot_utils.R"), local = TRUE)
    run_one_cs_fit_mix_plots(
      completed, expected_coding, file.path(output_dir, "expression"), "plot_summary.csv",
      project_dir = context$project_dir, datadir = datadir,
      result_reader = function(file) em_join_gene(file.path(context$iteration_dir, "results", basename(file)), context),
      mixed_fit_name = "weighted_fit_mix", mixed_label = paste0("EM ", context$iteration, " SuSiE-mix"))
  }
  invisible(list(comparisons = comparisons, summary = stats, overview = overview))
}
