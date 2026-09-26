# Rscript --vanilla script/analysis/tests/test_em_add1_mix2_plots.R
# Real EM discovery, joining, summary generation, Source entry point, and PNG
# rendering. Only GTEx/PLINK data reconstruction is replaced with synthetic data.
library(data.table)
source("script/analysis/em_analysis_utils.R")
source("script/analysis/generate_summary_results.R")
source("script/analysis/fit_mix_expression_plot_utils.R")
source("script/analysis/add1_mix2_plot_utils.R")

root <- tempfile("em-add1-mix2-")
dir.create(file.path(root, "script/scan_tissue_attempt"), recursive = TRUE)
dir.create(file.path(root, "script/analysis"), recursive = TRUE)
invisible(file.copy("script/scan_tissue_attempt/em_utils.R", file.path(root, "script/scan_tissue_attempt")))
invisible(file.copy(list.files("script/analysis", "\\.R$", full.names = TRUE), file.path(root, "script/analysis")))
dir.create(file.path(root, "results"))

make_fit <- function(p, leads, components) {
  alpha <- matrix(1 / p, nrow = 3, ncol = p)
  for (i in seq_along(leads)) {
    alpha[components[i], ] <- .02 / (p - 1L)
    alpha[components[i], leads[i]] <- .98
  }
  structure(list(alpha = alpha, pip = 1 - apply(1 - alpha, 2, prod),
    elbo = c(-120, -95, -100), KL = c(1, 2, 3), converged = TRUE,
    sets = list(cs = setNames(as.list(as.integer(leads)), paste0("L", components)),
      cs_index = as.integer(components), requested_coverage = .95,
      purity = data.frame(min.abs.corr = rep(1, length(leads)),
                          mean.abs.corr = 1, median.abs.corr = 1))), class = "susie")
}
p <- 60L
snps <- paste0("chr1_", 10000000L + seq_len(p) * 1000L, "_A_G_b38_A")
add_map <- data.frame(predictor_index = seq_len(p), predictor_name = snps,
                      snp = snps, coding = "additive")
codings <- rep(c("additive", "recessive", "dominant"), each = p)
mix_map <- data.frame(predictor_index = seq_len(3L * p),
                      predictor_name = paste(rep(snps, 3), codings, sep = "__"),
                      snp = rep(snps, 3), coding = codings)
fit_add <- make_fit(p, 2L, 3L)
old_one <- make_fit(3L * p, 7L, 3L)
old_two <- make_fit(3L * p, c(7L, 10L), c(3L, 1L))
em_two <- make_fit(3L * p, c(p + 2L, 2L * p + 40L), c(3L, 1L))
em_two$elbo <- c(-110, -90)
prior <- c(additive = .8, recessive = .05, dominant = .15)
genes <- c("EM_TWO", "OLD_TWO", "LOW_READ")
baseline <- list(susie_add = fit_add, susie_add_perm = fit_add,
                 susie_mix = old_one, susie_mix_perm = old_one,
                 weighted_fit_mix = old_two, add_predictor_map = add_map,
                 mix_predictor_map = mix_map, min_pv = 1e-10,
                 mean_read = 200, median_read = 150, n_ind = 180,
                 n_SNP = p, n_mix_predictor = 3L * p,
                 n_add_rm = 0, n_rec_rm = 0, n_dom_rm = 0)
iteration_dir <- file.path(root, "results_em/iteration_012")
dir.create(file.path(iteration_dir, "results"), recursive = TRUE)
dir.create(file.path(iteration_dir, "completed"))
priors <- data.frame(iteration = 12L, tissue = "Test tissue", pi_add = .8, pi_rec = .05, pi_dom = .15)
write.csv(priors, file.path(iteration_dir, "priors.csv"), row.names = FALSE)
write.csv(priors, file.path(root, "results_em/prior_history.csv"), row.names = FALSE)
write.csv(data.frame(chunk = 1L, gene = genes), file.path(iteration_dir, "manifest.csv"), row.names = FALSE)
saveRDS(list(iteration = 12L, chunk = 1L, n_genes = 3L), file.path(iteration_dir, "completed/chunk_001.done"))
for (gene in genes) {
  old <- baseline
  if (gene == "OLD_TWO") old$susie_mix <- old_two
  if (gene == "LOW_READ") old$mean_read <- 99
  saveRDS(list(`Test tissue` = old), file.path(root, "results", paste0(gene, ".rds")))
  current <- list(weighted_fit_mix = if (gene == "OLD_TWO") old_one else em_two,
                   mix_predictor_map = mix_map, mix_coding_prior = prior,
                   weighted_mix_prior_weights = unname(prior[codings] / p),
                   mean_read = old$mean_read, median_read = 150, n_ind = 180)
  out <- list(`Test tissue` = current)
  attr(out, "em_iteration") <- 12L
  saveRDS(out, file.path(iteration_dir, "results", paste0(gene, ".rds")))
}
summaries <- em_generate_summary(root)
selected <- select_add1_mix2_cases(summaries$res_summary,
  mixed_count_column = "ncs_weighted_fit_mix", overlap_column = "overlap_snp_add_weighted")
stopifnot(identical(selected$gene, "EM_TWO"),
          identical(select_add1_mix2_cases(summaries$res_summary)$gene, "OLD_TWO"))

# Stub only the external inputs; assert reconstruction receives the EM fit.
set.seed(92)
geno <- matrix(rbinom(180L * p, 2, .35), ncol = p, dimnames = list(NULL, snps))
y <- as.numeric(scale(2 * (geno[, 2] == 2) + 1.6 * (geno[, 40] >= 1) + rnorm(180, sd = .5)))
saveRDS(list(geno = geno, y = y, predictor_map = mix_map), file.path(root, "synthetic_inputs.rds"))
cat('\nload_fit_mix_expression_inputs <- function(gene_annotation_function, ...) {
  project <- dirname(dirname(dirname(gene_annotation_function)))
  readRDS(file.path(project, "synthetic_inputs.rds"))
}
prepare_fit_mix_gene_data <- function(gene_name, shared_inputs, results_dir, ...) {
  list(gene_results = readRDS(file.path(results_dir, paste0(gene_name, ".rds"))),
       synthetic = shared_inputs)
}
prepare_fit_mix_tissue_data <- function(gene_data, tissue_result, ...) {
  stopifnot(identical(tissue_result$susie_mix, tissue_result$weighted_fit_mix))
  x <- gene_data$synthetic
  list(y = x$y, geno_for_counts = x$geno, geno_mix = recode_fit_mix_predictors(x$geno),
       predictor_map = x$predictor_map)
}
', file = file.path(root, "script/analysis/fit_mix_expression_plot_utils.R"), append = TRUE)
cat('\noriginal_plot_four_panel_case <- plot_four_panel_case
plot_four_panel_case <- function(...) {
  args <- list(...)
  saveRDS(args, file.path(dirname(args$output_file), "render_inputs.rds"))
  original_plot_four_panel_case(...)
}
', file = file.path(root, "script/analysis/add1_mix2_plot_utils.R"), append = TRUE)

old_project <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
entry_env <- new.env(parent = globalenv())
sys.source("script/analysis/plot_additive_1cs_vs_mix_2cs_em.R", envir = entry_env)
if (is.na(old_project)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else Sys.setenv(SUSIE_MIX_PROJECT_DIR = old_project)
out <- entry_env$em_add1_mix2_results
audit <- out$plot_summary
plot_dir <- file.path(iteration_dir, "plot/additive_1cs_vs_mix_2cs")
stopifnot(out$context$iteration == 12L, nrow(audit) == 1L,
          audit$status == "plotted", audit$em_iteration == 12L,
          audit$mixed_fit_name == "weighted_fit_mix",
          audit$additive_lead_snp == snps[2], audit$cs1_lead_snp == snps[2],
          audit$cs2_lead_snp == snps[40], audit$cs1_lead_coding == "recessive",
          audit$cs2_lead_coding == "dominant", audit$cs1_component == 3L,
          audit$cs2_component == 1L, audit$closest_lead_distance_bp == 0,
          audit$cs1_to_cs2_distance_bp == 38000, audit$likelihood_statistic == 20,
          audit$cs1_adjusted_for == mix_map$predictor_name[2L * p + 40L],
          audit$cs2_adjusted_for == mix_map$predictor_name[p + 2L],
          file.info(audit$output_file)$size > 10000, grDevices::dev.cur() == 1L,
          !dir.exists(file.path(root, "plot")))
saved_audit <- fread(file.path(plot_dir, "additive_1cs_vs_mix_2cs_plot_summary.csv"))
stopifnot(saved_audit$em_iteration == 12L, saved_audit$cs2_lead_snp == snps[40],
          fread(file.path(plot_dir, "selected_cases.csv"))$gene == "EM_TWO")
rendered <- readRDS(file.path(plot_dir, "render_inputs.rds"))
x1 <- as.numeric(geno[, 2] == 2)
x2 <- as.numeric(geno[, 40] >= 1)
joint <- lm.fit(cbind(1, x1, x2), y)$coefficients
stopifnot(identical(rendered$fit_mix, em_two), rendered$mixed_label == "EM 12 SuSiE-mix",
          isTRUE(all.equal(rendered$conditional_expression_cs1, unname(y - joint[3] * x2))),
          isTRUE(all.equal(rendered$conditional_expression_cs2, unname(y - joint[2] * x1))))

preview_dir <- file.path("tmp", "em_add1_mix2_plot_validation")
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
preview_file <- file.path(preview_dir, "synthetic_em_four_panel.png")
stopifnot(file.copy(audit$output_file, preview_file, overwrite = TRUE))

# No selected cases should write an empty audit without needing GTEx input.
empty <- entry_env$em_plot_add1_mix2(root, minimum_mean_reads = 1000)
stopifnot(nrow(empty$plot_summary) == 0L,
          nrow(fread(file.path(plot_dir, "selected_cases.csv"))) == 0L,
          nrow(fread(file.path(plot_dir, "additive_1cs_vs_mix_2cs_plot_summary.csv"))) == 0L)

# Missing EM fits and stale CS counts are errors, never old-fit fallbacks.
run_env <- new.env(parent = globalenv())
sys.source("script/analysis/add1_mix2_plot_utils.R", envir = run_env)
run_env$load_fit_mix_expression_inputs <- function(...) NULL
run_env$prepare_fit_mix_gene_data <- function(...) list(gene_results = list())
bad_cases <- data.table(gene = c("Missing", "Wrong count"), tissue = "Test tissue")
bad <- suppressWarnings(run_env$run_add1_mix2_plots(bad_cases,
  plot_dir = file.path(root, "bad-plots"), project_dir = root,
  mixed_fit_name = "weighted_fit_mix", em_iteration = 12L,
  result_reader = function(file) {
    old <- baseline
    old$weighted_fit_mix <- if (basename(file) == "Missing.rds") NULL else old_one
    list(`Test tissue` = old)
  }))
stopifnot(all(bad$status == "error"), all(bad$em_iteration == 12L),
          grepl("no weighted_fit_mix", bad$message[1]),
          grepl("saved fits have 1 and 1", bad$message[2]))
cat("EM additive-1-CS / mixed-2-CS integration tests passed. Preview:", preview_file, "\n")
