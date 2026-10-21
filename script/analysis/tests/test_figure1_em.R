# Base-R end-to-end tests. All data, including posterior alpha, are SYNTHETIC.
# Rscript --vanilla script/analysis/tests/test_figure1_em.R
source("script/analysis/em_analysis_utils.R")
source("script/analysis/em_cs_agreement_utils.R")
source("script/analysis/figure1_em_utils.R")
fails <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
outdir <- file.path("tmp", "figure1_em_v2_validation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
root <- tempfile("figure1-v2-")
dir.create(file.path(root, "script/analysis"), recursive = TRUE)
dir.create(file.path(root, "script/scan_tissue_attempt"), recursive = TRUE)
invisible(file.copy(file.path("script/analysis", c("em_analysis_utils.R", "em_cs_agreement_utils.R",
  "figure1_em_utils.R", "generate_summary_results.R")), file.path(root, "script/analysis")))
invisible(file.copy("script/scan_tissue_attempt/em_utils.R", file.path(root, "script/scan_tissue_attempt")))
idir <- file.path(root, "results_em/iteration_018")
for (path in c("results", "completed", "summary")) dir.create(file.path(idir, path), recursive = TRUE)
source("script/analysis/generate_summary_results.R")
tissues <- c("Adrenal Gland", "Adipose Tissue", "Blood", "Blood Vessel", "Brain", "Breast", "Colon",
  "Esophagus", "Heart", "Kidney", "Liver", "Lung", "Muscle", "Nerve", "Ovary", "Pancreas", "Pituitary",
  "Prostate", "Salivary Gland", "Skin", "Small Intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina")
priors <- data.frame(iteration = 18L, tissue = tissues, pi_add = .8, pi_rec = .05, pi_dom = .15)
pos <- c(1000000, 1002000, 1030000, 1150000, 1450000, 1000000)
snps <- paste0("chr19_", sprintf("%.0f", pos), "_A_G_b38_A"); snps[6] <- "chr19_1000000_C_T_b38_C"
classes <- c("additive", "recessive", "dominant")
map <- data.frame(predictor_index = 1:18, predictor_name = paste(rep(snps, 3), rep(classes, each = 6), sep = "__"),
  snp = rep(snps, 3), coding = rep(classes, each = 6))
regions <- sets <- list()
expected <- list()
save_case <- function(gene, tissue, nadd, nem, a_snps = seq_len(nadd), e_snps = seq_len(nem),
                      profile = 1L, pvalue = 1e-12, reads = 200) {
  regions[[length(regions) + 1L]] <<- data.frame(gene = gene, tissue = tissue, min_pv = pvalue, mean_count = reads,
    ncs_susie = nadd, ncs_weighted_fit_mix = nem, has_weighted_fit_mix = TRUE, em_iteration = 18L)
  # Components 4 then 1 report CSs. Component 2 has extreme recessive support
  # but no CS and MUST NOT enter panel D. Component numbering is not CS numbering.
  alpha <- matrix(0, 4, 18)
  alpha[2, 9] <- 1; alpha[3, 5] <- 1
  cssets <- list(); components <- head(c(4L, 1L), nem)
  for (j in seq_len(nem)) {
    probs <- if (profile == 0L) c(.6, .36, .04) else
      if (profile == 1L) c(.94, .02, .04) else
      if (profile == 2L) c(.03, .96, .01) else c(.03, .01, .96)
    indices <- e_snps[j] + c(0, 6, 12)
    alpha[components[j], indices] <- probs
    members <- if (profile %in% 0:1) indices[1:2] else indices[which.max(probs)]
    cssets[[j]] <- as.integer(members)
    expected[[length(expected) + 1L]] <<- data.frame(gene = gene, tissue = tissue, cs_number = j,
      additive = probs[1], recessive = probs[2], dominant = probs[3])
  }
  # Fill non-reported rows so the saved fit remains a probability model.
  for (j in which(rowSums(alpha) == 0)) alpha[j, 1] <- 1
  fit <- list(alpha = alpha, pip = 1 - apply(1 - alpha, 2, prod), converged = TRUE,
    sets = list(cs = setNames(cssets, if (length(components)) paste0("L", components) else character()), cs_index = components))
  for (model in c("susie_add", "weighted_fit_mix")) {
    nn <- if (model == "susie_add") nadd else nem
    for (j in seq_len(nn)) {
      idx <- if (model == "susie_add") a_snps[j] else get_cs_lead_predictor(fit, cssets, j)
      id <- if (model == "susie_add") snps[idx] else map$snp[idx]
      members <- if (model == "susie_add") id else unique(map$snp[cssets[[j]]])
      sets[[length(sets) + 1L]] <<- data.frame(gene = gene, tissue = tissue, model_key = model,
        cs_number = j, component = if (model == "susie_add") j else components[j],
        lead_snp = id, lead_coding = if (model == "susie_add") "additive" else map$coding[idx],
        lead_pip = if (model == "susie_add") .9 else fit$pip[idx],
        cs_snps = paste(members, collapse = ";"), cs_size_snps = length(members),
        distance_to_tss_kb = (em_agreement_coordinates(id)$pos - 1100000) / 1000, em_iteration = 18L)
    }
  }
  list(weighted_fit_mix = fit, mix_predictor_map = map,
       mix_coding_prior = c(additive = .8, recessive = .05, dominant = .15))
}
genes <- sprintf("G%02d", 1:12)
for (g in seq_along(genes)) {
  saved <- list()
  for (t in seq_along(tissues)) {
    na <- if (g <= (t %% 5)) 2L else 1L
    ne <- if (g == 6L && t %% 3L == 0) 2L else 1L
    if (t == 27L || g == 11L) ne <- 0L
    if (t == 1L && g == 12L) { na <- 0L; ne <- 1L }
    a <- head(c(1L, 5L), na)
    e <- head(c(if (g %% 3 == 0) 4L else 1L, 3L), ne)
    if (t == 1L && g == 2L) e <- 6L  # Other allele, same position.
    pr <- if (g == 1L && t == 1L) 0L else if (g == 7L && t %% 2 == 0) 2L else if (g == 8L) 3L else 1L
    saved[[tissues[t]]] <- save_case(genes[g], tissues[t], na, ne, a, e, profile = pr)
  }
  attr(saved, "em_iteration") <- 18L
  saveRDS(saved, file.path(idir, "results", paste0(genes[g], ".rds")))
}
for (gene in c("BoundaryP", "LowReads", "BothZero")) {
  saved <- list(Blood = save_case(gene, "Blood", if (gene == "BothZero") 0L else 1L,
    if (gene == "BothZero") 0L else 1L, pvalue = if (gene == "BoundaryP") 1e-8 else 1e-12,
    reads = if (gene == "LowReads") 99 else 200))
  # Excluded genes deliberately have invalid fits, proving D uses cohort selection.
  saved$Blood$weighted_fit_mix <- NULL
  attr(saved, "em_iteration") <- 18L
  saveRDS(saved, file.path(idir, "results", paste0(gene, ".rds")))
}
res_summary <- do.call(rbind, regions); res_cs_summary <- do.call(rbind, sets)
# Missing TSS in one extra additive CS must not change any cohort or D denominator.
res_cs_summary$distance_to_tss_kb[res_cs_summary$gene == "G01" & res_cs_summary$tissue == "Blood" &
  res_cs_summary$model_key == "susie_add" & res_cs_summary$cs_number == 2] <- NA_real_
res_errors <- cs_errors <- data.frame()
save(res_summary, res_errors, file = file.path(idir, "summary/res_summary.RData"))
save(res_cs_summary, cs_errors, file = file.path(idir, "summary/res_cs_summary.RData"))
saveRDS(list(iteration = 18L, priors = priors, generated_at = "SYNTHETIC"), file.path(idir, "summary/em_summary_metadata.rds"))
manifest <- data.frame(gene = c(genes, "BoundaryP", "LowReads", "BothZero"), chunk = 1L)
write.csv(manifest, file.path(idir, "manifest.csv"), row.names = FALSE)
write.csv(priors, file.path(idir, "priors.csv"), row.names = FALSE)
saveRDS(list(iteration = 18L, chunk = 1L, n_genes = nrow(manifest)), file.path(idir, "completed/chunk_001.done"))
future <- priors; future$iteration <- 19L
dir.create(file.path(root, "results_em/iteration_019"))
write.csv(manifest, file.path(root, "results_em/iteration_019/manifest.csv"), row.names = FALSE)
write.csv(rbind(priors, future), file.path(root, "results_em/prior_history.csv"), row.names = FALSE)

# Execute the actual Source entry, editing only its user-editable settings block.
entry <- parse("script/analysis/plot_figure1_em.R")
env <- new.env(parent = globalenv()); eval(entry[[1L]], env)
env$figure1_settings$project_dir <- root
env$figure1_settings$minimum_tissue_regions <- 10L
env$figure1_settings$output_dir <- outdir
env$figure1_settings$show_plot <- FALSE
env$figure1_settings$preview_label <- "SYNTHETIC VALIDATION - NOT GTEx RESULTS"
for (expr in as.list(entry)[-1L]) suppressWarnings(eval(expr, env))
z <- env$figure1_result
settings <- env$figure1_settings
context <- suppressWarnings(em_analysis_context(root))
stopifnot(z$iteration == 18L, z$metrics$n_regions == 324, z$metrics$n_tissues == 27,
  !dir.exists(file.path(root, "results")), is.null(z$tissue_coding_priors),
  nrow(z$selected_tissues) == 4L, nrow(z$coding_tissue_support) == 27L,
  all(z$cohort$gene %in% genes), nrow(z$coding_region_audit) == nrow(z$cohort),
  sum(z$cohort$n_em) == nrow(z$coding_cs_support),
  sum(z$coding_tissue_support$n_em_cs) == nrow(z$coding_cs_support),
  z$metrics$count_agreement_rate == mean(z$cohort$n_add == z$cohort$n_em),
  z$metrics$lead_retention_rate == mean(z$matched_cs_pairs$same_lead_snp),
  sum(z$tissue_count_table$n_regions) == nrow(z$cohort),
  all(z$matched_cs_pairs$region_id %in% z$cohort$region_id),
  all(z$unpaired_cs$region_id %in% z$cohort$region_id),
  all(z$coding_cs_support$region_id %in% z$cohort$region_id),
  file.info(file.path(outdir, "figure1_em.png"))$size > 10000,
  file.info(file.path(outdir, "figure1_em.pdf"))$size > 1000, dev.cur() == 1L)
# The zero-EM tissue remains present, with undefined support rather than 0%.
zero <- z$coding_tissue_support[z$coding_tissue_support$tissue == "Vagina", ]
stopifnot(zero$n_cohort_regions == 12L, zero$n_em_cs == 0, is.na(zero$mean_recessive),
  all(z$coding_region_audit$status[z$coding_region_audit$tissue == "Vagina"] == "no_em_cs"))
# D includes EM-only and multi-CS regions, and uses complete alpha probabilities.
target <- z$coding_cs_support[z$coding_cs_support$gene == "G01" & z$coding_cs_support$tissue == tissues[1], ]
stopifnot(target$component == 4L, abs(target$posterior_additive - .6) < 1e-12,
  abs(target$posterior_recessive - .36) < 1e-12, abs(target$posterior_dominant - .04) < 1e-12,
  abs(target$posterior_mass_inside_cs - .96) < 1e-12,
  all(abs(rowSums(z$coding_cs_support[paste0("posterior_", classes)]) - 1) < 1e-12),
  any(z$coding_cs_support$gene == "G12" & z$coding_cs_support$tissue == tissues[1]),
  all(z$coding_cs_support$component %in% c(4L, 1L)))
expected <- do.call(rbind, expected)
expected <- expected[expected$gene %in% genes, ]
key <- function(d) paste(d$gene, d$tissue, d$cs_number)
ix <- match(key(z$coding_cs_support), key(expected))
stopifnot(!anyNA(ix), isTRUE(all.equal(unname(as.matrix(z$coding_cs_support[paste0("posterior_", classes)])),
  unname(as.matrix(expected[ix, classes])), check.attributes = FALSE)))
for (t in tissues) {
  d <- z$coding_cs_support[z$coding_cs_support$tissue == t, ]
  row <- z$coding_tissue_support[z$coding_tissue_support$tissue == t, ]
  if (nrow(d)) stopifnot(abs(row$mean_recessive - mean(d$posterior_recessive)) < 1e-12)
}
stopifnot(abs(z$coding_overall_support$mean_recessive - mean(z$coding_cs_support$posterior_recessive)) < 1e-12)
strict_tables <- settings; strict_tables$minimum_tissue_regions <- 100L
strict <- suppressWarnings(summarize_figure1_em(res_summary, res_cs_summary, priors, strict_tables))
stopifnot(nrow(strict$selected_tissues) == 0L, identical(strict$cohort, z$cohort), identical(strict$metrics, z$metrics))
# Preserve the preferred B/C groups and distinguish nearest distances from pairing.
ref <- summarize_em_cs_agreement(res_summary, res_cs_summary, minimum_tissue_regions = 10, em_priors = priors)
stopifnot(identical(z$distances, ref$distances), identical(z$tss_distances, ref$tss_distances),
  identical(z$selected_tissues, ref$selected_tissues), identical(z$tissue_cs_count_tables, ref$tissue_cs_count_tables),
  identical(z$matched_cs_pairs, ref$matched_cs_pairs),
  all(z$tss_distance_summary$n_panel_b == z$lead_distance_summary$n_distance[2:3]),
  sum(z$tss_distance_summary$n_missing_tss) == 1L)
p <- z$matched_cs_pairs[z$matched_cs_pairs$gene == "G02" & z$matched_cs_pairs$tissue == tissues[1], ]
stopifnot(!p$same_lead_snp, p$distance_bp == 0)
steps <- fig1_ecdf_steps(c(0, 0, 10, NA), 6)
stopifnot(identical(steps$x[1:2], c(0, 0)), steps$y[2] == 2/3, tail(steps$y, 1) == 1,
  nrow(fig1_ecdf_steps(numeric(), 6)) == 0)

# Fail loudly on stale summaries, wrong iteration/prior, or malformed alpha.
one <- z$cohort[z$cohort$gene == "G01" & z$cohort$tissue == tissues[1], ]
fpath <- file.path(idir, "results/G01.rds"); original <- readRDS(fpath)
bad <- original; attr(bad, "em_iteration") <- 17L; saveRDS(bad, fpath)
fails(fig1_collect_coding_support(context, one, res_cs_summary), "Invalid selected-iteration")
bad <- original; bad[[1]]$mix_coding_prior <- c(additive = .9, recessive = .05, dominant = .05); saveRDS(bad, fpath)
fails(fig1_collect_coding_support(context, one, res_cs_summary), "prior does not match")
bad <- original; bad[[1]]$weighted_fit_mix$alpha[4, ] <- 0; saveRDS(bad, fpath)
fails(fig1_collect_coding_support(context, one, res_cs_summary), "summing to one")
bad <- original; bad[[1]]$weighted_fit_mix$sets$cs_index <- 2L; saveRDS(bad, fpath)
fails(fig1_collect_coding_support(context, one, res_cs_summary), "summary disagree")
saveRDS(original, fpath)
badcs <- res_cs_summary; badcs$component[badcs$gene == "G01" & badcs$tissue == tissues[1]] <- 2L
fails(fig1_collect_coding_support(context, one, badcs), "summary disagree")
badcs <- res_cs_summary; badcs$em_iteration[1] <- 17L
fails(summarize_figure1_em(res_summary, badcs, priors, settings), "one EM iteration")
fails(em_analysis_context(root, 19L), "unfinished")
saveRDS(list(iteration = 17L, priors = priors), file.path(idir, "summary/em_summary_metadata.rds"))
fails(em_load_summary(context), "provenance mismatch")
empty_settings <- settings; empty_settings$association_threshold <- 1e-20
fails(suppressWarnings(summarize_figure1_em(res_summary, res_cs_summary, priors, empty_settings)), "No regions pass")

# Exported v2 data replot even with invalid raw metadata; v1 data are rejected.
replot <- settings; replot$plot_data_file <- normalizePath(file.path(outdir, "figure1_plot_data.rds"), winslash = "/")
replot$output_dir <- file.path(outdir, "replot")
again <- run_figure1_em(replot)
stopifnot(identical(again$metrics, z$metrics), identical(again$coding_cs_support, z$coding_cs_support),
  identical(again$preview_label, z$preview_label))
replot$iteration <- 19L
fails(run_figure1_em(replot), "differs from the exported")
legacy <- tempfile(fileext = ".rds"); saveRDS(list(schema = "figure1_em_v1"), legacy)
replot$plot_data_file <- legacy; replot$iteration <- "latest"
fails(run_figure1_em(replot), "Older/incompatible")

# Empty ECDFs and no EM components remain valid; no phantom probabilities.
single_res <- res_summary[res_summary$gene == "G01" & res_summary$tissue == "Vagina", ]
single_cs <- res_cs_summary[res_cs_summary$gene == "G01" & res_cs_summary$tissue == "Vagina", ]
single <- suppressWarnings(summarize_figure1_em(single_res, single_cs, priors, settings))
single <- fig1_add_coding_summary(single, fig1_collect_coding_support(context, single$cohort, single_cs))
single$preview_label <- settings$preview_label
plot_figure1_em(single, file.path(outdir, "no_pairs_no_em.png"))
stopifnot(single$metrics$n_pairs == 0, is.na(single$metrics$lead_retention_rate),
  single$coding_overall_support$n_em_cs == 0, is.na(single$coding_overall_support$mean_additive))
pdf(file.path(outdir, "small_rstudio_preview.pdf"), width = 7, height = 5)
plot_figure1_em(z); invisible(dev.off())
stopifnot(dev.cur() == 1L)
cat("PASS: v2 Source entry, shared cohort, headlines, unchanged ECDFs, full-alpha coding probabilities, provenance failures, empty groups, replot and PNG/PDF rendering.\n")
cat("Synthetic preview:", normalizePath(file.path(outdir, "figure1_em.png"), winslash = "/"), "\n")

