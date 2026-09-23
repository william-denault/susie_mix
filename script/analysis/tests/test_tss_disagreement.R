# Rscript --vanilla script/analysis/tests/test_tss_disagreement.R
source("script/analysis/tss_disagreement_utils.R")
fixture <- data.frame(
  gene = c(rep("Partial", 4), rep("CodingOnly", 2), "AddOnly", "MixOnly",
           "CrossTissue", "CrossTissue", "Missing", "Missing"),
  tissue = c(rep("Tissue", 8), "Tissue", "Other", "Tissue", "Tissue"),
  model = c("SuSiE", "SuSiE", "SuSiE-mix", "SuSiE-mix", "SuSiE", "SuSiE-mix",
            "SuSiE", "SuSiE-mix", "SuSiE", "SuSiE-mix", "SuSiE", "SuSiE-mix"),
  lead_snp = c("s1", "s2", "s1", "s3", "s4", "s4", "s5", "s6", "s7", "s7", "s8", NA),
  lead_coding = c("additive", "additive", "dominant", "recessive", "additive", "dominant",
                  "additive", "recessive", "additive", "dominant", "additive", "dominant"),
  distance_to_tss_kb = c(0, -20, 0, 40, 0, 0, -300, NA, 30, 30, 5, 5))
selected <- select_tss_disagreement(fixture)
stopifnot(identical(selected$audit$tss_include,
                    c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)),
          all(selected$selected$lead_snp %in% c("s2", "s3", "s5", "s6", "s7")),
          identical(selected$audit$tss_agreement_status[11:12], c("unknown_other_model_lead", "unknown_lead")))
counts <- selected$summary[selected$summary$tissue == "All tissues", ]
stopifnot(all(counts$n_cs_agreed_excluded == 2L), all(counts$n_cs_selected == 3L),
          counts$n_selected_outside_window[counts$model == "SuSiE"] == 1L,
          counts$n_selected_missing_distance[counts$model == "SuSiE-mix"] == 1L)
empty <- select_tss_disagreement(fixture[FALSE, ])
stopifnot(nrow(empty$selected) == 0L, all(empty$summary$n_cs_selected == 0L))
all_agreed <- select_tss_disagreement(fixture[5:6, ])
stopifnot(nrow(all_agreed$selected) == 0L)

if (requireNamespace("data.table", quietly = TRUE)) {
  source("script/analysis/generate_summary_results.R")
  source("script/analysis/descriptive_results_weighted.R")
  project_source <- normalizePath(".", winslash = "/")
  root <- tempfile("tss-disagreement-")
  dir.create(file.path(root, "results"), recursive = TRUE)
  dir.create(file.path(root, "script/analysis"), recursive = TRUE)
  file.copy("script/analysis/tss_disagreement_utils.R", file.path(root, "script/analysis"))
  snps <- paste0("chr1_", c(100000, 130000, 170000), "_A_G_b38_A")
  add_map <- data.frame(predictor_index = 1:3, predictor_name = snps, snp = snps, coding = "additive")
  coding <- rep(c("additive", "recessive", "dominant"), each = 3)
  mix_map <- data.frame(predictor_index = 1:9, predictor_name = paste(rep(snps, 3), coding, sep = "__"),
                        snp = rep(snps, 3), coding = coding)
  make_fit <- function(p, leads) {
    alpha <- matrix(1 / p, nrow = 2, ncol = p)
    for (i in seq_along(leads)) { alpha[i, ] <- .02 / (p - 1); alpha[i, leads[i]] <- .98 }
    list(alpha = alpha, pip = 1 - apply(1 - alpha, 2L, prod), elbo = c(-20, -10), KL = c(1, 1),
         converged = TRUE, sets = if (length(leads)) list(
           cs = setNames(lapply(leads, as.integer), paste0("L", seq_along(leads))),
           cs_index = seq_along(leads)) else NULL)
  }
  scenarios <- list(Partial = list(c(1, 2), c(1, 9)), CodingOnly = list(1, 7),
                    Changed = list(1, 8), AddOnly = list(1, integer()),
                    MixOnly = list(integer(), 9), NoCS = list(integer(), integer()))
  for (gene in names(scenarios)) {
    leads <- scenarios[[gene]]
    x <- list(susie_add = make_fit(3, leads[[1]]), susie_mix = make_fit(9, leads[[2]]),
              weighted_fit_mix = make_fit(9, leads[[2]]), susie_add_perm = make_fit(3, integer()),
              susie_mix_perm = make_fit(9, integer()), add_predictor_map = add_map,
              mix_predictor_map = mix_map, min_pv = 1e-10, mean_read = 200, median_read = 150,
              n_ind = 100, n_SNP = 3, n_mix_predictor = 9, n_add_rm = 0, n_rec_rm = 0, n_dom_rm = 0,
              tss = 120000, chromosome = "chr1", strand = "+")
    saveRDS(list(Tissue = x), file.path(root, "results", paste0(gene, ".rds")))
  }
  generated <- generate_summary_results(file.path(root, "results"), root)
  weighted <- run_weighted_descriptive_results(file.path(root, "res_summary.RData"),
                                               file.path(root, "res_cs_summary.RData"),
                                               file.path(root, "weighted"), project_dir = root)
  check <- function(histogram, audit) {
    stopifnot(sum(histogram$count) == 6L, all(histogram$n_cs_total == 3L),
              all(histogram$n_cs_in_window == 3L),
              all(abs(tapply(histogram$proportion, histogram$model, sum) - 1) < 1e-10),
              sum(audit$tss_shared_lead) == 4L, sum(audit$tss_include) == 6L,
              sum(audit$tss_include[audit$gene == "Partial"]) == 2L,
              !any(audit$tss_include[audit$gene == "CodingOnly"]))
  }
  check(weighted$weighted_tss_distance_distribution, weighted$weighted_tss_cs_agreement_audit)
  # Execute the original script, including its overall and tissue TSS panels.
  saved_project <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
  old_env <- new.env(parent = globalenv())
  pdf(file.path(root, "interactive-previews.pdf"))
  source(file.path(project_source, "script/analysis/descriptive results.R"), local = old_env)
  dev.off()
  if (is.na(saved_project)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else Sys.setenv(SUSIE_MIX_PROJECT_DIR = saved_project)
  check(old_env$tss_distance_summary, old_env$tss_selection$audit)
  stopifnot(sum(old_env$tissue_tss_distance_summary$count) == 6L,
            sum(weighted$weighted_tissue_tss_distance_distribution$count) == 6L,
            nrow(old_env$cs_idx) == 10L,
            file.exists(file.path(root, "descriptive_results/tss_cs_agreement_audit.csv")))
  # All-agreed input must yield an explicitly empty plot rather than old curves.
  res_summary <- generated$res_summary[generated$res_summary$gene == "CodingOnly", , drop = FALSE]
  res_cs_summary <- generated$res_cs_summary[generated$res_cs_summary$gene == "CodingOnly", , drop = FALSE]
  save(res_summary, file = file.path(root, "res_summary.RData"))
  save(res_cs_summary, file = file.path(root, "res_cs_summary.RData"))
  all_shared <- run_weighted_descriptive_results(file.path(root, "res_summary.RData"),
                                                 file.path(root, "res_cs_summary.RData"),
                                                 file.path(root, "all_shared"), project_dir = root)
  stopifnot(all(all_shared$weighted_tss_distance_distribution$count == 0L),
            all(is.na(all_shared$weighted_tss_distance_distribution$proportion)))
  message("Descriptive integration output: ", root)
} else message("SKIP descriptive integration: data.table is unavailable")
cat("PASS: per-CS agreement, coding-only agreement, partial regions, one-model CSs,\n",
    "tissue isolation, unknown leads, distance denominators, and empty disagreement sets.\n")
