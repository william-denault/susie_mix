# Reproduce the matched-CS distance/PIP figure for additive SuSiE versus the
# latest COMPLETE EM SuSiE-mix iteration. Base R only; no refitting or GTEx load.
#
# RSTUDIO: open this file, edit settings if needed, and click Source.
# First generate that iteration's summaries with generate_summary_results_em.R.
# CLI (from the project root):
# Rscript --vanilla script/analysis/plot_matched_cs_lead_distances_em.R [PROJECT_DIR] [latest|ITERATION] [TISSUE]
#
# Outputs: selected iteration/plot/matched_cs_lead_distances/ (PNG, PDF, CSVs).
# All tissues, equal positive CS counts, no extra P-value/read/coding filters.
# The mixed-model source is weighted_fit_mix, NOT the original susie_mix fit.
em_matched_cs_settings <- list(
  project_dir = Sys.getenv("SUSIE_MIX_PROJECT_DIR", if (dir.exists("script/analysis")) getwd() else
                            "/project2/mstephens/wdenault/susie_mix"),
  iteration = "latest", # Or an explicit completed iteration, e.g. 18.
  tissue = NULL,        # NULL = all tissues; e.g. "Blood" for a single tissue.
  output_dir = NULL,    # NULL = the selected iteration's default plot folder.
  show_plot = interactive()
)

em_matched_cs_load_helpers <- function(project_dir) {
  helpers <- new.env(parent = globalenv())
  source(file.path(project_dir, "script/analysis/em_analysis_utils.R"), local = helpers)
  old <- options(susie_mix.matched_cs.autorun = FALSE)
  on.exit(options(old), add = TRUE)
  source(file.path(project_dir, "script/analysis/plot_matched_cs_lead_distances.R"), local = helpers)
  helpers
}

run_em_matched_cs_lead_distances <- function(settings = em_matched_cs_settings) {
  helpers <- em_matched_cs_load_helpers(settings$project_dir)
  context <- helpers$em_analysis_context(settings$project_dir, settings$iteration)
  summaries <- helpers$em_load_summary(context)
  tissue <- settings$tissue
  if (!is.null(tissue) && (length(tissue) != 1L || is.na(tissue) ||
                           !tissue %in% summaries$res$tissue)) {
    stop("Choose a tissue in the EM summary. Available: ",
         paste(sort(unique(summaries$res$tissue)), collapse = ", "))
  }
  message("Matching additive SuSiE against EM weighted_fit_mix (iteration ", context$iteration, ")...")
  result <- helpers$match_equal_cs_regions(
    summaries$res, summaries$cs,
    mixed_model_key = "weighted_fit_mix", mixed_count_column = "ncs_weighted_fit_mix")
  # Keep source identities explicit even in empty tables. mixed_* fields are EM.
  for (nm in c("pairs", "region_audit", "orphan_cs")) {
    result[[nm]]$em_iteration <- rep(context$iteration, nrow(result[[nm]]))
    result[[nm]]$comparison_mixed_model <- rep("weighted_fit_mix", nrow(result[[nm]]))
  }
  pairs <- result$pairs
  selected <- if (is.null(tissue)) pairs else pairs[pairs$tissue == tissue, , drop = FALSE]
  distance_summary <- helpers$matched_cs_distance_summary(
    selected, if (is.null(tissue)) "All tissues" else tissue)
  per_tissue <- lapply(sort(unique(selected$tissue)), function(t) {
    helpers$matched_cs_distance_summary(selected[selected$tissue == t, , drop = FALSE], t)
  })
  per_tissue <- if (length(per_tissue)) do.call(rbind, per_tissue) else distance_summary[FALSE, ]
  distance_summary$em_iteration <- context$iteration
  per_tissue$em_iteration <- rep(context$iteration, nrow(per_tissue))
  out <- settings$output_dir
  if (is.null(out)) {
    out <- file.path(context$iteration_dir, "plot", "matched_cs_lead_distances")
    if (!is.null(tissue)) out <- file.path(out, "by_tissue", gsub("[^[:alnum:]_-]", "_", tissue))
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tables <- list(matched_cs_pairs = pairs, region_matching_audit = result$region_audit,
                 orphan_cs_rows = result$orphan_cs, distance_summary = distance_summary,
                 distance_summary_by_tissue = per_tissue)
  for (nm in names(tables)) write.csv(tables[[nm]], file.path(out, paste0(nm, ".csv")),
                                     row.names = FALSE, na = "")
  draw <- function(output_file = NULL) helpers$plot_matched_cs_distances(
    pairs, output_file, tissue, mixed_label = "EM SuSiE-mix",
    analysis_label = paste("EM iteration", context$iteration))
  for (ext in c("png", "pdf")) draw(file.path(out, paste0("matched_cs_lead_distances.", ext)))
  metadata <- list(
    iteration = context$iteration, iteration_dir = context$iteration_dir,
    mixed_model_key = "weighted_fit_mix", mixed_count_column = "ncs_weighted_fit_mix",
    region_summary = file.path(context$summary_dir, "res_summary.RData"),
    cs_summary = file.path(context$summary_dir, "res_cs_summary.RData"),
    summary_metadata = summaries$metadata, settings = settings,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  saveRDS(metadata, file.path(out, "analysis_metadata.rds"))
  writeLines(c(
    paste("Additive SuSiE / EM SuSiE-mix lead distances: iteration", context$iteration),
    paste("Region summary:", metadata$region_summary),
    paste("CS summary:", metadata$cs_summary), paste("Generated:", metadata$generated_at),
    paste("Plot tissue:", if (is.null(tissue)) "All tissues" else tissue),
    "Latest selects the latest complete iteration using manifests, completion markers, result files and priors.",
    "Completion does not itself establish convergence of the EM prior updates. Summary provenance is checked.",
    "Mixed-model source: weighted_fit_mix; counts: ncs_weighted_fit_mix. Original susie_mix is not substituted.",
    "Cohort: all available summary regions with equal positive additive/EM CS counts; no further association, read-count, convergence or coding filters.",
    "Summary-generation omissions remain documented in summary/res_errors.csv and summary/res_cs_errors.csv.",
    "Matching is identical to the original figure: sort GRCh38 lead positions and pair ranks, minimizing total one-to-one absolute distance.",
    "Ties minimize sum of squared distances, then use SNP ID / CS number. This does not maximize unchanged SNPs or establish signal equivalence.",
    "One row per matched CS pair; a k-CS region contributes k observations. region_weight = 1/k is supplied but not used in this plot.",
    "Coordinates are parsed from the actual lead_snp. Invalid or incomplete regions are excluded in full and audited.",
    "Agreement compares SNP IDs; distance compares positions. Different variants may have 0 bp distance.",
    "Y is the EM lead predictor's coding-specific PIP, not its component alpha or PIP summed across codings.",
    "Zero-distance pairs with valid PIPs form the left marginal density (Gaussian reflection at 0/1, nrd0 bandwidth, unit area).",
    "Positive distances with valid PIPs form the right scatter (logarithmic kb); color is EM lead coding, shape is CS count.",
    "Missing/invalid PIPs are retained in distance and agreement counts. All-pair percentages include zeros; cutoffs are strictly >10 kb and >100 kb.",
    "matched_cs_pairs.csv and region_matching_audit.csv retain all tissues; the optional tissue setting selects the figure and distance summaries.",
    "All mixed_* pair columns refer to the selected EM fit. comparison_mixed_model and em_iteration record this explicitly."
  ), file.path(out, "matching_notes.txt"))
  print(table(result$region_audit$status))
  print(distance_summary[c("em_iteration", "tissue", "n_regions", "n_pairs", "n_same_lead_snp",
                           "n_different_lead_snp", "n_over_10kb", "n_over_100kb")], row.names = FALSE)
  if (nrow(result$orphan_cs)) warning(nrow(result$orphan_cs), " orphan CS rows; see orphan_cs_rows.csv.")
  message("Saved EM distance/support figure and tables to ", normalizePath(out, winslash = "/"))
  if (interactive() && isTRUE(settings$show_plot)) tryCatch(draw(), error = function(e) {
    message("Files saved; plot preview unavailable: ", conditionMessage(e))
  })
  invisible(c(result, list(distance_summary = distance_summary,
                           distance_summary_by_tissue = per_tissue,
                           metadata = metadata, output_dir = out)))
}

# Helper-only imports/tests can set options(susie_mix.em_matched_cs.autorun = FALSE).
if (isTRUE(getOption("susie_mix.em_matched_cs.autorun", TRUE))) {
  em_matched_cs_args <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
  if ("--help" %in% em_matched_cs_args) {
    cat("Usage: Rscript --vanilla script/analysis/plot_matched_cs_lead_distances_em.R [PROJECT_DIR] [latest|ITERATION] [TISSUE]\n",
        "First run generate_summary_results_em.R for the same iteration. Source in RStudio uses the settings above.\n",
        "Default output: results_em/iteration_NNN/plot/matched_cs_lead_distances/\n", sep = "")
  } else {
    if (length(em_matched_cs_args) > 3L) stop("Expected [PROJECT_DIR] [latest|ITERATION] [TISSUE].")
    if (length(em_matched_cs_args) >= 1L) em_matched_cs_settings$project_dir <- em_matched_cs_args[1]
    if (length(em_matched_cs_args) >= 2L) em_matched_cs_settings$iteration <- em_matched_cs_args[2]
    if (length(em_matched_cs_args) >= 3L) em_matched_cs_settings$tissue <- em_matched_cs_args[3]
    em_matched_cs_results <- run_em_matched_cs_lead_distances(em_matched_cs_settings)
    em_matched_cs_pairs <- em_matched_cs_results$pairs
  }
}
