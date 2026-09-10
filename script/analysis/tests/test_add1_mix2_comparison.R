# Run from the repository root with base R:
# Rscript --vanilla script/analysis/tests/test_add1_mix2_comparison.R
# Load only function definitions, without GTEx I/O or package loading. The
# PIP renderer is stubbed for the preview; the changed figure and shared
# expression renderer, metric and distance calculation are production code.
plot_env <- new.env(parent = globalenv())
load_functions <- function(path, wanted, envir) {
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) {
      eval(expr, envir)
    }
  }
  stopifnot(all(vapply(wanted, exists, logical(1), envir = envir, inherits = FALSE)))
}
load_functions(
  "script/analysis/plot_additive_1cs_vs_mix_2cs.R",
  c("get_add1_mix2_lead_distances", "plot_four_panel_case"), plot_env
)
load_functions(
  "script/analysis/fit_mix_expression_plot_utils.R",
  c("get_one_cs_likelihood_comparison", "finite_plot_limits", "open_four_panel_png",
    "add_mix_coding_boundaries", "plot_fit_mix_expression_panel"), plot_env
)
summary_env <- new.env(parent = globalenv())
load_functions("script/analysis/generate_summary_results.R", "get_log_lik_metric", summary_env)

add_lead <- data.frame(lead_snp = "chr1_1000000_A_G_b38_A")
leads <- data.frame(
  cs_name = c("L7", "L2"),
  lead_snp = c("chr1_1010000_C_T_b38_C", "chr1_999800_T_G_b38_T"),
  lead_coding = c("recessive", "dominant")
)
distances <- plot_env$get_add1_mix2_lead_distances(add_lead, leads)
stopifnot(distances$additive_to_cs1_distance_bp == 10000,
          distances$additive_to_cs2_distance_bp == 200,
          distances$cs1_to_cs2_distance_bp == 10200,
          distances$closest_lead_distance_bp == 200,
          distances$closest_mix_cs_names == "L2",
          distances$closest_mix_lead_snps == leads$lead_snp[2],
          distances$closest_mix_lead_codings == "dominant",
          !distances$closest_lead_tie,
          distances$additive_lead_chromosome == "chr1")
swapped <- plot_env$get_add1_mix2_lead_distances(add_lead, leads[2:1, ])
stopifnot(swapped$closest_lead_distance_bp == 200,
          swapped$closest_mix_cs_names == "L2", swapped$additive_to_cs1_distance_bp == 200)

# Tie handling must preserve both CS names, independent of their ordering.
tied_leads <- leads
tied_leads$lead_snp <- c("chr1_1001000_C_T_b38_C", "chr1_999000_T_G_b38_T")
tied <- plot_env$get_add1_mix2_lead_distances(add_lead, tied_leads)
stopifnot(tied$closest_lead_distance_bp == 1000, tied$closest_lead_tie,
          tied$closest_mix_cs_names == "L7;L2")
same_leads <- leads
same_leads$lead_snp[] <- add_lead$lead_snp
same <- plot_env$get_add1_mix2_lead_distances(add_lead, same_leads)
stopifnot(same$closest_lead_distance_bp == 0, same$cs1_to_cs2_distance_bp == 0,
          same$closest_lead_tie, same$closest_mix_cs_names == "L7;L2")
for (bad_snp in c(NA_character_, "invalid", "chr2_999800_T_G_b38_T",
                  "chr1_999800_T_G_b37_T", "chr1_0_T_G_b38_T")) {
  bad_leads <- leads
  bad_leads$lead_snp[2] <- bad_snp
  bad <- plot_env$get_add1_mix2_lead_distances(add_lead, bad_leads)
  stopifnot(bad$additive_to_cs1_distance_bp == 10000,
            is.na(bad$additive_to_cs2_distance_bp),
            is.na(bad$closest_lead_distance_bp), is.na(bad$closest_mix_cs_names),
            grepl("unavailable", bad$lead_distance_status))
}

fit_add <- list(elbo = c(-120, -95, -100), KL = c(1, NA_real_, 2),
                sets = list(cs = list(L1 = 2L)), pip = c(.02, .96, .02))
fit_mix <- list(elbo = c(-110, -90), KL = c(1, 2),
                sets = list(cs = list(L7 = 5L, L2 = 7L)),
                pip = c(.01, .01, .01, .01, .85, .01, .75, .01, .01))
comparison <- plot_env$get_one_cs_likelihood_comparison(fit_add, fit_mix)
stopifnot(comparison$log_lik_add == -97, comparison$log_lik_mix == -87,
          comparison$likelihood_statistic == 20,
          comparison$likelihood_statistic == 2 *
            (summary_env$get_log_lik_metric(fit_mix) - summary_env$get_log_lik_metric(fit_add)))
negative <- plot_env$get_one_cs_likelihood_comparison(fit_mix, fit_add)
stopifnot(negative$likelihood_statistic == -20)
missing_fit <- fit_mix
missing_fit$KL <- NULL
missing <- plot_env$get_one_cs_likelihood_comparison(fit_add, missing_fit)
stopifnot(is.na(missing$likelihood_statistic))

output_dir <- file.path("tmp", "add1_mix2_comparison_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Test that numeric distances and both tied CS labels survive CSV export.
audit_file <- file.path(output_dir, "comparison_fields.csv")
write.csv(as.data.frame(c(comparison, tied)), audit_file, row.names = FALSE)
saved <- read.csv(audit_file)
stopifnot(saved$closest_lead_distance_bp == 1000,
          saved$closest_mix_cs_names == "L7;L2", saved$likelihood_statistic == 20)

plot_env$susie_plot <- function(fit, y, main) {
  plot(seq_along(fit$pip), fit$pip, pch = 16, ylim = c(0, 1),
       xlab = "variable", ylab = "PIP", main = main)
}
captured_labels <- character()
plot_env$mtext <- function(text, ...) {
  captured_labels <<- c(captured_labels, as.character(text))
  graphics::mtext(text, ...)
}
predictor_map <- data.frame(predictor_index = 1:9,
                            coding = rep(c("additive", "recessive", "dominant"), each = 3))
set.seed(73)
g1 <- sample(rep(0:2, c(80, 60, 40)))
g2 <- sample(rep(0:2, c(75, 70, 35)))
e1 <- 1.6 * (g1 == 2) + rnorm(length(g1), sd = .6)
e2 <- 1.3 * (g2 >= 1) + rnorm(length(g2), sd = .6)
render <- function(name, lead_set, metric, lead_distances) {
  captured_labels <<- character()
  output <- file.path(output_dir, paste0(name, ".png"))
  plot_env$plot_four_panel_case(
    gene_name = "Synthetic comparison", tissue_name = "Test tissue",
    fit_add = fit_add, fit_mix = fit_mix, predictor_map = predictor_map,
    leads = lead_set, add_lead = add_lead,
    likelihood_comparison = metric, lead_distances = lead_distances,
    raw_genotype_cs1 = g1, raw_genotype_cs2 = g2,
    conditional_expression_cs1 = e1, conditional_expression_cs2 = e2,
    output_file = output
  )
  stopifnot(file.info(output)$size > 10000, grDevices::dev.cur() == 1L)
}
render("nearest_cs2", leads, comparison, distances)
stopifnot(any(grepl("additive -> L2 = 200 bp", captured_labels, fixed = TRUE)),
          any(grepl("= 20.00", captured_labels, fixed = TRUE)),
          any(grepl("no calibrated p-value", captured_labels, fixed = TRUE)))
render("tied_negative", tied_leads, negative, tied)
stopifnot(any(grepl("L7 / L2 = 1,000 bp (GRCh38) [tie]", captured_labels, fixed = TRUE)),
          any(grepl("= -20.00", captured_labels, fixed = TRUE)))
render("unavailable", leads, missing, bad)
stopifnot(any(grepl("Closest lead distance unavailable", captured_labels, fixed = TRUE)),
          any(grepl("Likelihood-ratio statistic unavailable", captured_labels, fixed = TRUE)))
cat("Additive-1-CS / mixed-2-CS distance, likelihood and figure tests passed.\n")
