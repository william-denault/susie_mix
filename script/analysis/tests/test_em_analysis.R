# Rscript --vanilla script/analysis/tests/test_em_analysis.R
# Synthetic iteration 12; no GTEx or refitting. Optional data.table checks.
source("script/analysis/em_analysis_utils.R")
source("script/analysis/em_one_cs_utils.R")
source("script/analysis/generate_summary_results.R")
source("script/analysis/descriptive_results_weighted.R")
root <- tempfile("em-analysis-")
dir.create(file.path(root, "script/scan_tissue_attempt"), recursive = TRUE)
dir.create(file.path(root, "script/analysis"), recursive = TRUE)
file.copy("script/scan_tissue_attempt/em_utils.R", file.path(root, "script/scan_tissue_attempt"))
file.copy(list.files("script/analysis", "\\.R$", full.names = TRUE), file.path(root, "script/analysis"))
dir.create(file.path(root, "results"))
make_fit <- function(p, cs, lead = 1L) {
  alpha <- matrix(rep(.02 / max(1, p - 1), p), nrow = 1)
  alpha[1, lead] <- .98
  list(alpha = alpha, pip = as.numeric(alpha), elbo = c(-20, -10), KL = 2,
       converged = TRUE, sets = if (cs) list(cs = list(L1 = as.integer(lead)), cs_index = 1L) else NULL)
}
snps <- c("chr1_100000_A_G_b38_A", "chr1_350000_C_T_b38_C")
add_map <- data.frame(predictor_index = 1:2, predictor_name = snps, snp = snps, coding = "additive")
coding <- rep(c("additive", "recessive", "dominant"), each = 2)
mix_map <- data.frame(predictor_index = 1:6, predictor_name = paste(rep(snps, 3), coding, sep = "__"),
                      snp = rep(snps, 3), coding = coding)
baseline <- list(susie_add = make_fit(2, TRUE), susie_add_perm = make_fit(2, FALSE),
                 susie_mix = make_fit(6, TRUE), susie_mix_perm = make_fit(6, FALSE),
                 weighted_fit_mix = make_fit(6, TRUE, 3), # Old recessive fit must be replaced.
                 add_predictor_map = add_map, mix_predictor_map = mix_map,
                 min_pv = 1e-10, mean_read = 200, median_read = 150, n_ind = 100,
                 n_SNP = 2, n_mix_predictor = 6, n_add_rm = 0, n_rec_rm = 0, n_dom_rm = 0)
prior <- c(additive = .8, recessive = .05, dominant = .15)
current <- list(weighted_fit_mix = make_fit(6, TRUE, 6), mix_predictor_map = mix_map,
                mix_coding_prior = prior, weighted_mix_prior_weights = unname(prior[coding] / 2),
                mean_read = 200, median_read = 150, n_ind = 100)
manifest <- data.frame(chunk = 1L, gene = c("DOM", "ZERO", "FAIL"))
for (gene in manifest$gene) saveRDS(list(Tissue = baseline, Dropped = baseline),
                                  file.path(root, "results", paste0(gene, ".rds")))
history <- NULL
for (iteration in c(9L, 12L, 13L)) {
  d <- file.path(root, "results_em", sprintf("iteration_%03d", iteration))
  dir.create(file.path(d, "results"), recursive = TRUE)
  dir.create(file.path(d, "completed"))
  priors <- data.frame(iteration = iteration, tissue = "Tissue", pi_add = .8, pi_rec = .05, pi_dom = .15)
  history <- rbind(history, priors)
  write.csv(priors, file.path(d, "priors.csv"), row.names = FALSE)
  write.csv(manifest, file.path(d, "manifest.csv"), row.names = FALSE)
  if (iteration == 13) next # Prepared, not complete: latest must still be 12.
  for (gene in manifest$gene) {
    x <- current
    if (gene == "ZERO") x$weighted_fit_mix <- make_fit(6, FALSE)
    out <- if (gene == "FAIL") list(gene = gene, error = "Synthetic gene failure") else list(Tissue = x)
    attr(out, "em_iteration") <- iteration
    if (gene == "DOM") attr(out, "tissue_errors") <- list(Dropped = "Synthetic tissue failure")
    saveRDS(out, file.path(d, "results", paste0(gene, ".rds")))
  }
  saveRDS(list(iteration = iteration, chunk = 1L, n_genes = 3L), file.path(d, "completed/chunk_001.done"))
}
write.csv(history, file.path(root, "results_em/prior_history.csv"), row.names = FALSE)
expect_error <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}
context <- suppressWarnings(em_analysis_context(root))
stopifnot(context$iteration == 12L)
expect_error(em_analysis_context(root, 13), "unfinished")
expect_error(em_analysis_context(root, 11), "absent")
out <- suppressWarnings(em_generate_summary(root))
stopifnot(nrow(out$res_summary) == 2L, all(out$res_summary$em_iteration == 12L),
          out$res_summary$n_dom_weighted[out$res_summary$gene == "DOM"] == 1L,
          out$res_summary$ncs_weighted_fit_mix[out$res_summary$gene == "ZERO"] == 0L,
          !any(out$res_summary$n_rec_weighted), nrow(out$res_errors) == 3L,
          any(grepl("Synthetic tissue failure", out$res_errors$error)),
          !file.exists(file.path(root, "res_summary.RData")))
em_load_summary(context)
joined <- em_join_gene(file.path(context$iteration_dir, "results/DOM.rds"), context)
stopifnot(joined$Tissue$min_pv == baseline$min_pv,
          identical(joined$Tissue$weighted_fit_mix, current$weighted_fit_mix))
bad <- current
bad$mix_predictor_map$snp[1] <- snps[2]
bad_out <- list(Tissue = bad)
attr(bad_out, "em_iteration") <- 12L
bad_path <- file.path(context$iteration_dir, "results/DOM.rds")
good_out <- readRDS(bad_path)
saveRDS(bad_out, bad_path)
expect_error(em_join_gene(bad_path, context), "predictor maps differ")
bad_out <- good_out
attr(bad_out, "em_iteration") <- 9L
saveRDS(bad_out, bad_path)
expect_error(em_join_gene(bad_path, context), "wrong em_iteration")
saveRDS(good_out, bad_path)
plots <- em_plot_one_cs(root, 12, expression_plots = FALSE)
stopifnot(nrow(plots$comparisons) == 1L, plots$comparisons$distance_bp == 250000,
          plots$comparisons$lead_snp == snps[2], plots$comparisons$selected_cs_overlap_snps == 0,
          plots$summary$max_distance_kb == 250,
          file.info(plots$overview$output_file)$size > 10000)
empty <- em_plot_one_cs(root, 12, expected_coding = "recessive", expression_plots = FALSE)
stopifnot(nrow(empty$comparisons) == 0L, empty$summary$n_compared == 0L)
# Cover the additive no-CS fallback.
no_cs <- baseline$susie_add
no_cs$sets <- NULL
lead <- em_one_cs_lead(no_cs, add_map)
stopifnot(lead$snp == snps[1], grepl("no credible set", lead$selection))
if (requireNamespace("data.table", quietly = TRUE)) {
  tables <- run_weighted_descriptive_results(file.path(context$summary_dir, "res_summary.RData"),
                                             file.path(context$summary_dir, "res_cs_summary.RData"),
                                             file.path(context$iteration_dir, "descriptive_results"))
  stopifnot(nrow(tables$weighted_primary_analysis_set) == 2L,
            file.exists(file.path(context$iteration_dir, "descriptive_results/weighted_overall_summary.csv")))
} else message("SKIP descriptive integration: data.table not installed")
# Exercise the dedicated recessive Source entry against a real saved EM fit.
# Stub the GTEx-dependent renderer; test_one_cs_plots.R checks actual panels.
writeLines(c(
  'run_one_cs_fit_mix_plots <- function(cases, expected_coding, plot_dir, summary_filename,',
  '                                  result_reader, mixed_fit_name, ...) {',
  '  stopifnot(nrow(cases) == 1L, mixed_fit_name == "weighted_fit_mix")',
  '  saved <- result_reader(paste0(cases$gene[1], ".rds"))',
  '  stopifnot(!is.null(saved[[cases$tissue[1]]][[mixed_fit_name]]))',
  '  data.frame(status = "plotted", gene = cases$gene, tissue = cases$tissue,',
  '             coding = expected_coding, mixed_fit_name = mixed_fit_name)',
  '}'
), file.path(root, "script/analysis/fit_mix_expression_plot_utils.R"))
recessive_out <- good_out
recessive_out$Tissue$weighted_fit_mix <- make_fit(6, TRUE, 4L)
saveRDS(recessive_out, bad_path)
suppressWarnings(em_generate_summary(root, 12))
saved_project <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
recessive_env <- new.env(parent = globalenv())
suppressWarnings(source("script/analysis/finding_interesting_recessive_1cs_em.R", local = recessive_env))
if (is.na(saved_project)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else Sys.setenv(SUSIE_MIX_PROJECT_DIR = saved_project)
stopifnot(isTRUE(recessive_env$em_one_cs_settings$expression_plots),
          recessive_env$em_one_cs_results$expression_summary$status == "plotted",
          recessive_env$em_one_cs_results$expression_summary$coding == "recessive",
          recessive_env$em_one_cs_results$summary$n_compared == 1L,
          recessive_env$em_one_cs_results$comparisons$lead_coding == "recessive",
          recessive_env$em_one_cs_results$comparisons$distance_bp == 250000,
          nrow(recessive_env$em_one_cs_results$overview$annotated_cases) == 1L)
cat("PASS: latest complete iteration, EM/baseline joining, audit records, zero-CS fits, provenance,\n",
    "250-kb lead shift, empty selection, and downstream outputs. Fixtures: ", root, "\n", sep = "")
