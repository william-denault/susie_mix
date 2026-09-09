# Run from the repository root with data.table, matrixStats and susieR installed:
# Rscript --vanilla script/analysis/tests/test_one_cs_plots.R
# Synthetic fits/data exercise lead mapping, summary-metric agreement and the
# complete plotting/audit loop. No GTEx inputs or PLINK invocation are needed.

plot_env <- new.env(parent = globalenv())
sys.source("script/analysis/fit_mix_expression_plot_utils.R", envir = plot_env)

summary_env <- new.env(parent = globalenv())
for (expr in parse("script/analysis/generate_summary_results.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("get_log_lik_metric"))) {
    eval(expr, summary_env)
  }
}

expect_error <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}

make_fit <- function(p, lead, cs_index = 3L) {
  alpha <- matrix(1 / p, nrow = 3L, ncol = p)
  alpha[cs_index, ] <- 0.02 / (p - 1L)
  alpha[cs_index, lead] <- 0.98
  structure(list(
    alpha = alpha,
    pip = 1 - apply(1 - alpha, 2L, prod),
    elbo = c(-120, -95, -100),
    KL = c(1, NA_real_, 2),
    sets = list(
      cs = setNames(list(as.integer(lead)), paste0("L", cs_index)),
      cs_index = cs_index,
      purity = data.frame(min.abs.corr = 1, mean.abs.corr = 1,
                          median.abs.corr = 1),
      requested_coverage = 0.95
    )
  ), class = "susie")
}

p <- 60L
snps <- paste0("chr1_", 10000000L + seq_len(p) * 1000L, "_A_G_b38_A")
add_map <- data.frame(predictor_index = seq_len(p), predictor_name = snps,
                      snp = snps, coding = "additive")
codings <- rep(c("additive", "recessive", "dominant"), each = p)
mix_map <- data.frame(
  predictor_index = seq_len(3L * p),
  predictor_name = paste(rep(snps, 3L), codings, sep = "__"),
  snp = rep(snps, 3L), coding = codings
)
add_fit <- make_fit(p, 2L)
rec_fit <- make_fit(3L * p, p + 2L)
dom_fit <- make_fit(3L * p, 2L * p + 4L)
rec_fit$elbo <- c(-110, -90)
dom_fit$elbo <- c(-110, -92)

add_lead <- plot_env$get_fit_add_plot_lead(add_fit, add_map)
rec_lead <- plot_env$get_fit_mix_cs_leads(rec_fit, mix_map)
dom_lead <- plot_env$get_fit_mix_cs_leads(dom_fit, mix_map)
stopifnot(
  add_lead$lead_snp == snps[2], add_lead$component_index == 3L,
  rec_lead$lead_snp == snps[2], rec_lead$lead_coding == "recessive",
  dom_lead$lead_snp == snps[4], dom_lead$lead_coding == "dominant",
  grepl("Same lead SNP", plot_env$describe_one_cs_lead_change(add_lead, rec_lead)),
  grepl("Different lead SNP", plot_env$describe_one_cs_lead_change(add_lead, dom_lead)),
  identical(plot_env$get_fit_add_plot_lead(add_fit, add_map[p:1, ]), add_lead)
)

# A CS lead is defined by its component alpha, even if another member's PIP
# is larger because it also appears in other components.
alpha_fit <- add_fit
alpha_fit$sets$cs[[1]] <- c(2L, 3L)
alpha_fit$pip[3] <- 1
stopifnot(plot_env$get_fit_add_plot_lead(alpha_fit, add_map)$lead_snp == snps[2])

no_cs_fit <- add_fit
no_cs_fit$sets <- list(cs = NULL)
no_cs_lead <- plot_env$get_fit_add_plot_lead(no_cs_fit, add_map)
stopifnot(no_cs_lead$lead_snp == snps[2], is.na(no_cs_lead$cs_name),
          no_cs_lead$lead_selection == "Highest-PIP SNP (no credible set)")

multi_fit <- add_fit
multi_fit$alpha[1, ] <- 0.02 / (p - 1L)
multi_fit$alpha[1, 5] <- 0.98
multi_fit$pip[5] <- 0.999
multi_fit$sets$cs <- list(L3 = 2L, L1 = 5L)
multi_fit$sets$cs_index <- c(3L, 1L)
multi_fit$sets$purity <- multi_fit$sets$purity[c(1, 1), ]
multi_lead <- plot_env$get_fit_add_plot_lead(multi_fit, add_map)
stopifnot(multi_lead$lead_snp == snps[5], multi_lead$component_index == 1L)

bad_map <- add_map
bad_map$predictor_index[2] <- 1L
expect_error(plot_env$get_fit_add_plot_lead(add_fit, bad_map), "does not match")
expect_error(plot_env$get_fit_add_plot_lead(add_fit, NULL), "absent")

# Check against the actual summary helper, including decreasing ELBO, missing
# KL, non-finite ELBO, and a negative statistic (never truncate it to zero).
missing_kl <- rec_fit
missing_kl$KL <- NULL
nonfinite <- rec_fit
nonfinite$elbo <- c(NA_real_, Inf)
for (fit in list(rec_fit, dom_fit, missing_kl, nonfinite, add_fit)) {
  comparison <- plot_env$get_one_cs_likelihood_comparison(add_fit, fit)
  expected_add <- summary_env$get_log_lik_metric(add_fit)
  expected_mix <- summary_env$get_log_lik_metric(fit)
  stopifnot(
    identical(comparison$log_lik_add, expected_add),
    identical(comparison$log_lik_mix, expected_mix),
    identical(comparison$likelihood_statistic, -2 * (expected_add - expected_mix))
  )
}
stopifnot(plot_env$get_one_cs_likelihood_comparison(add_fit, rec_fit)$likelihood_statistic == 20,
          plot_env$get_one_cs_likelihood_comparison(rec_fit, add_fit)$likelihood_statistic == -20)

# Only input loading/reconstruction is stubbed. The production runner selects
# both leads, calls the installed susieR plotter, and writes the audit CSV.
set.seed(52)
geno <- matrix(rbinom(180L * p, 2, 0.35), ncol = p,
               dimnames = list(NULL, snps))
noise <- rnorm(nrow(geno), sd = 0.55)
fixture <- function(add, mix, coding) list(
  tissue_result = list(susie_add = add, susie_mix = mix,
                       add_predictor_map = add_map, mix_predictor_map = mix_map),
  y = as.numeric(scale(if (coding == "recessive") {
    2.5 * (geno[, 2] == 2L) + noise
  } else {
    1.8 * (geno[, 4] >= 1L) + noise
  }))
)
fixtures <- list(
  `Synthetic recessive` = fixture(add_fit, rec_fit, "recessive"),
  `Synthetic dominant` = fixture(add_fit, dom_fit, "dominant"),
  `Synthetic no additive CS` = fixture(no_cs_fit, rec_fit, "recessive"),
  `Synthetic multiple additive CSs` = fixture(multi_fit, dom_fit, "dominant"),
  `Synthetic missing metric` = fixture(add_fit, missing_kl, "recessive"),
  `Missing additive SNP` = fixture(add_fit, rec_fit, "recessive"),
  `Wrong coding` = fixture(add_fit, dom_fit, "dominant")
)
fixtures[["Missing additive SNP"]]$tissue_result$add_predictor_map$snp[2] <- "absent"
plot_env$load_fit_mix_expression_inputs <- function(...) NULL
plot_env$prepare_fit_mix_gene_data <- function(gene_name, ...) {
  list(gene_results = list(`Small Intestine` = fixtures[[gene_name]]$tissue_result))
}
plot_env$prepare_fit_mix_tissue_data <- function(gene_name, ...) {
  list(y = fixtures[[gene_name]]$y, geno_for_counts = geno,
       geno_mix = plot_env$recode_fit_mix_predictors(geno), predictor_map = mix_map)
}
output_dir <- file.path("tmp", "one_cs_plot_validation")
run_cases <- function(genes, coding) plot_env$run_one_cs_fit_mix_plots(
  cases = data.table::data.table(gene = genes, tissue = "Small Intestine"),
  expected_coding = coding,
  plot_dir = file.path(output_dir, coding),
  summary_filename = "plot_summary.csv",
  project_dir = output_dir
)
rec_rows <- run_cases(c("Synthetic recessive", "Synthetic no additive CS",
                        "Synthetic missing metric", "Missing additive SNP", "Wrong coding"),
                      "recessive")
dom_rows <- run_cases(c("Synthetic dominant", "Synthetic multiple additive CSs"), "dominant")
stopifnot(
  identical(rec_rows$status, c("plotted", "plotted", "plotted", "error", "error")),
  all(dom_rows$status == "plotted"),
  all(rec_rows$same_lead_snp[1:3]), !any(dom_rows$same_lead_snp),
  rec_rows$additive_lead_snp[1] == snps[2], dom_rows$lead_snp[1] == snps[4],
  rec_rows$likelihood_statistic[1] == 20,
  is.na(rec_rows$likelihood_statistic[3]),
  rec_rows$n_genotype_2[1] == sum(geno[, 2] == 2L),
  dom_rows$additive_n_genotype_0[1] == sum(geno[, 2] == 0L),
  all(file.info(c(rec_rows$output_file[1:3], dom_rows$output_file))$size > 10000),
  file.exists(file.path(output_dir, "recessive", "plot_summary.csv")),
  grDevices::dev.cur() == 1L
)
empty_rows <- plot_env$run_one_cs_fit_mix_plots(
  cases = data.table::data.table(gene = character(), tissue = character()),
  expected_coding = "recessive", plot_dir = file.path(output_dir, "empty"),
  summary_filename = "plot_summary.csv", project_dir = output_dir
)
stopifnot(nrow(empty_rows) == 0L)
cat("One-CS plot tests passed. Synthetic previews:", output_dir, "\n")
