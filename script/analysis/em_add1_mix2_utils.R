# Compare the original additive fit with two-CS fits from a completed EM iteration.
em_plot_add1_mix2 <- function(
    project_dir, iteration = "latest", datadir = "/project2/mstephens/gtex",
    association_threshold = 1e-8, minimum_mean_reads = 100) {
  context <- em_analysis_context(project_dir, iteration)
  source(file.path(context$project_dir, "script/analysis/generate_summary_results.R"), local = TRUE)
  source(file.path(context$project_dir, "script/analysis/fit_mix_expression_plot_utils.R"), local = TRUE)
  source(file.path(context$project_dir, "script/analysis/add1_mix2_plot_utils.R"), local = TRUE)
  library(susieR)
  # Resolve join validation against the locally sourced summary helpers.
  em_join_gene <- em_join_gene
  environment(em_join_gene) <- environment()
  summaries <- em_load_summary(context)
  cases <- select_add1_mix2_cases(
    summaries$res, association_threshold, minimum_mean_reads,
    mixed_count_column = "ncs_weighted_fit_mix",
    overlap_column = "overlap_snp_add_weighted"
  )
  plot_dir <- file.path(context$iteration_dir, "plot", "additive_1cs_vs_mix_2cs")
  message("Using EM iteration ", context$iteration, "; plots: ", plot_dir)
  plot_summary <- run_add1_mix2_plots(
    cases = cases, plot_dir = plot_dir, project_dir = context$project_dir,
    datadir = datadir, results_dir = context$baseline_dir,
    temp_dir = file.path(context$iteration_dir, "temp_plink_add1_mix2"),
    result_reader = function(file) em_join_gene(
      file.path(context$iteration_dir, "results", basename(file)), context
    ),
    mixed_fit_name = "weighted_fit_mix",
    mixed_label = paste0("EM ", context$iteration, " SuSiE-mix"),
    em_iteration = context$iteration
  )
  invisible(list(context = context, cases = cases, plot_summary = plot_summary))
}
