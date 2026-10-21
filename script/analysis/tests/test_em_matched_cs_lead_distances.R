# Rscript --vanilla script/analysis/tests/test_em_matched_cs_lead_distances.R
# Synthetic saved EM summaries: no GTEx data, packages or model fitting.
options(susie_mix.em_matched_cs.autorun = FALSE)
source("script/analysis/plot_matched_cs_lead_distances_em.R")
helpers <- em_matched_cs_load_helpers(getwd())
expect_error <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
variant <- function(p) if (length(p)) paste0("chr1_", p, "_A_G_b38_A") else character()
regions <- sets <- list()
add_case <- function(gene, add, em, original = add, pip = .7, coding = "dominant", tissue = "Blood") {
  regions[[length(regions) + 1L]] <<- data.frame(
    gene = gene, tissue = tissue, ncs_susie = length(add), ncs_susie_mix = length(original),
    ncs_weighted_fit_mix = length(em), has_weighted_fit_mix = TRUE, em_iteration = 18L,
    min_pv = 1, mean_count = 1, converged_mix = FALSE, converged_weighted_fit_mix = TRUE,
    log_lik_mix = -1000, log_lik_weighted_mix = -10, result_file = paste0(gene, ".rds"))
  for (model in c("susie_add", "susie_mix", "weighted_fit_mix")) {
    pos <- switch(model, susie_add = add, susie_mix = original, weighted_fit_mix = em)
    if (!length(pos)) next
    sets[[length(sets) + 1L]] <<- data.frame(
      gene = gene, tissue = tissue, model_key = model, cs_number = seq_along(pos),
      lead_snp = variant(pos), lead_coding = if (model == "weighted_fit_mix") coding else "additive",
      lead_pip = if (model == "weighted_fit_mix") pip else .99,
      em_iteration = 18L, result_file = paste0(gene, ".rds"))
  }
}
add_case("SAME", 1000, 1000, pip = .95, coding = "additive")
add_case("SHIFT", 1000, 21000, pip = .35)
add_case("MULTI", c(100, 200), c(200, 300), pip = c(.2, .8))
add_case("EM_EQUAL", 500, 700, original = c(500, 800), pip = .6, tissue = "Liver")
add_case("BASE_EQUAL", 1000, c(1000, 2000))
add_case("MISSINGPIP", 10, 200010, pip = NA_real_, coding = "recessive")
add_case("INVALID", 10, 0)
add_case("ZERO", numeric(), numeric())
res_summary <- do.call(rbind, regions)
res_cs_summary <- do.call(rbind, sets)
res_errors <- cs_errors <- data.frame(error = character())

dir.create("tmp", showWarnings = FALSE)
root <- tempfile("em-matched-cs-", tmpdir = normalizePath("tmp", winslash = "/"))
dir.create(file.path(root, "script/analysis"), recursive = TRUE)
dir.create(file.path(root, "script/scan_tissue_attempt"), recursive = TRUE)
stopifnot(all(file.copy(file.path("script/analysis", c("em_analysis_utils.R",
  "plot_matched_cs_lead_distances.R", "plot_matched_cs_lead_distances_em.R")),
  file.path(root, "script/analysis"))))
stopifnot(file.copy("script/scan_tissue_attempt/em_utils.R", file.path(root, "script/scan_tissue_attempt")))
history <- NULL
manifest <- data.frame(gene = res_summary$gene, chunk = 1L)
for (iteration in 18:19) {
  d <- file.path(root, "results_em", sprintf("iteration_%03d", iteration))
  dir.create(file.path(d, "results"), recursive = TRUE)
  dir.create(file.path(d, "completed"))
  priors <- data.frame(iteration = iteration, tissue = c("Blood", "Liver"),
                       pi_add = .85, pi_rec = .05, pi_dom = .1)
  history <- rbind(history, priors)
  write.csv(manifest, file.path(d, "manifest.csv"), row.names = FALSE)
  write.csv(priors, file.path(d, "priors.csv"), row.names = FALSE)
  if (iteration == 19L) next # latest must skip this incomplete iteration.
  # Context checks the raw-file inventory; plotting itself reads only summaries.
  for (gene in manifest$gene) saveRDS(list(), file.path(d, "results", paste0(gene, ".rds")))
  saveRDS(list(iteration = iteration, chunk = 1L, n_genes = nrow(manifest)),
          file.path(d, "completed/chunk_001.done"))
  dir.create(file.path(d, "summary"))
  summary_dir <- file.path(d, "summary")
  save(res_summary, res_errors, file = file.path(summary_dir, "res_summary.RData"))
  save(res_cs_summary, cs_errors, file = file.path(summary_dir, "res_cs_summary.RData"))
  saveRDS(list(iteration = iteration, priors = priors), file.path(summary_dir, "em_summary_metadata.rds"))
}
write.csv(history, file.path(root, "results_em/prior_history.csv"), row.names = FALSE)
settings <- modifyList(em_matched_cs_settings, list(project_dir = root, show_plot = FALSE))
result <- suppressWarnings(run_em_matched_cs_lead_distances(settings))
p <- result$pairs; s <- result$distance_summary
status <- setNames(result$region_audit$status, result$region_audit$gene)
stopifnot(result$metadata$iteration == 18L, nrow(p) == 6L, s$n_regions == 5L,
          s$n_same_lead_snp == 1L, s$n_different_lead_snp == 5L,
          s$n_over_10kb == 2L, s$n_over_100kb == 1L, s$n_missing_or_invalid_mixed_pip == 1L,
          all(p$em_iteration == 18L), all(p$comparison_mixed_model == "weighted_fit_mix"),
          all(grepl("|weighted_fit_mix|", p$mixed_cs_id, fixed = TRUE)),
          all(p$converged_weighted_fit_mix), all(p$log_lik_weighted_mix == -10),
          !"converged_mix" %in% names(p), !"log_lik_mix" %in% names(p),
          p$distance_bp[p$gene == "SHIFT"] == 20000, p$mixed_lead_pip[p$gene == "SHIFT"] == .35,
          all(p$distance_bp[p$gene == "MULTI"] == 100), # Not identical-lead-first matching.
          status["EM_EQUAL"] == "matched", status["BASE_EQUAL"] == "unequal_cs_counts",
          status["INVALID"] == "invalid_lead_coordinates", status["ZERO"] == "equal_zero_cs",
          all(p$min_pv == 1), all(p$mean_count == 1))
expected <- c("matched_cs_pairs.csv", "region_matching_audit.csv", "orphan_cs_rows.csv",
              "distance_summary.csv", "distance_summary_by_tissue.csv", "analysis_metadata.rds",
              "matching_notes.txt", "matched_cs_lead_distances.png", "matched_cs_lead_distances.pdf")
stopifnot(all(file.info(file.path(result$output_dir, expected))$size > 0),
          readRDS(file.path(result$output_dir, "analysis_metadata.rds"))$iteration == 18L)

# Tissue selection affects the figure/summaries, while the pair audit stays pooled.
liver <- suppressWarnings(run_em_matched_cs_lead_distances(modifyList(settings, list(tissue = "Liver"))))
stopifnot(liver$distance_summary$n_pairs == 1L, nrow(liver$pairs) == 6L,
          grepl("by_tissue/Liver$", liver$output_dir))
expect_error(suppressWarnings(run_em_matched_cs_lead_distances(
  modifyList(settings, list(tissue = "Unknown")))), "Choose a tissue")
expect_error(run_em_matched_cs_lead_distances(modifyList(settings, list(iteration = 19))), "unfinished")

# Detect mixed/stale iterations instead of silently plotting the original fit.
res_cs_summary$em_iteration[1] <- 17L
save(res_cs_summary, cs_errors, file = file.path(summary_dir, "res_cs_summary.RData"))
expect_error(suppressWarnings(run_em_matched_cs_lead_distances(settings)), "mixed-iteration EM CS")
res_cs_summary$em_iteration[] <- 18L
save(res_cs_summary, cs_errors, file = file.path(summary_dir, "res_cs_summary.RData"))

# Source follows the same RStudio entry path without triggering the old script.
local({
  old_options <- options(susie_mix.em_matched_cs.autorun = TRUE)
  on.exit(options(old_options))
  saved <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
  on.exit(if (is.na(saved)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else
    Sys.setenv(SUSIE_MIX_PROJECT_DIR = saved), add = TRUE)
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
  env <- new.env(parent = globalenv())
  suppressWarnings(source("script/analysis/plot_matched_cs_lead_distances_em.R", local = env))
  stopifnot(env$em_matched_cs_results$metadata$iteration == 18L,
            identical(env$em_matched_cs_pairs, p),
            !file.exists(file.path(root, "plot/matched_cs_lead_distances")))
})
cat("PASS: EM source/count/PIP selection, latest complete iteration, provenance, cohort, matching, exports and RStudio Source.\n",
    "Synthetic fixture outputs: ", result$output_dir, "\n", sep = "")
