# Run from the project root; no GTEx or model fitting required.
source("script/analysis/slide_analysis_utils.R")
source("script/analysis/slide_descriptive_utils.R")
source("script/analysis/slide_one_cs_utils.R")
project <- normalizePath(".", winslash = "/")
root <- tempfile("slide-analysis-")
dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
snps <- sprintf("chr1_%d_A_G_b38_A", c(100000, 350000, 400000))
map <- data.frame(predictor_index = 1:3, predictor_name = snps, snp = snps, coding = "additive")
make_fit <- function(leads = integer(), components = seq_along(leads), deltas = rep(0, length(leads))) {
  alpha <- matrix(.01, 3, 3, dimnames = list(NULL, snps))
  delta <- matrix(0, 3, 3)
  for (i in seq_along(leads)) {
    alpha[components[i], leads[i]] <- .98
    delta[components[i], leads[i]] <- deltas[i]
  }
  list(alpha = alpha, pip = 1 - apply(1 - alpha, 2, prod), delta = delta,
    delta_forced = c(FALSE, FALSE, TRUE), elbo = c(-20, -10), KL = c(1, 1, 1), converged = TRUE,
    sets = list(cs = setNames(lapply(leads, as.integer), if (length(leads)) paste0("L", components) else character()), cs_index = components,
                coverage = rep(.98, length(leads))))
}
make_tissue <- function(add_leads = 1L, slide_leads = 2L, components = 3L, deltas = 1) {
  list(susie_add = make_fit(add_leads), susie_add_perm = make_fit(),
    fit_slide = make_fit(slide_leads, components, deltas), fit_slide_perm = make_fit(),
    add_predictor_map = map, min_pv = 1e-10, mean_read = 200, median_read = 150, n_ind = 100,
    # These offsets encode TSS=150000; the slide lead has distance +200 kb.
    susie_add_lead_snp_tss_distance = setNames(c(-50000), snps[1]))
}
fixtures <- list(
  Dominant = make_tissue(), Recessive = make_tissue(deltas = -1),
  PartialDominant = make_tissue(deltas = .4), PartialRecessive = make_tissue(deltas = -.6),
  SameLead = make_tissue(slide_leads = 1L),
  PartialAgreement = make_tissue(c(1L, 2L), c(1L, 3L), c(3L, 1L), c(1, -1)),
  NoAddCS = make_tissue(integer()), MultipleAddCS = make_tissue(c(1L, 3L)),
  NoSlideCS = make_tissue(slide_leads = integer(), components = integer(), deltas = numeric()),
  NoCS = make_tissue(integer(), integer(), integer(), numeric()),
  MissingSlide = make_tissue(), MissingPerm = make_tissue(), BadOrder = make_tissue())
fixtures$MissingSlide$fit_slide <- NULL
fixtures$MissingPerm$fit_slide_perm <- NULL
colnames(fixtures$BadOrder$fit_slide$alpha) <- rev(snps)
for (gene in names(fixtures)) saveRDS(list(Tissue = fixtures[[gene]]), file.path(root, "results", paste0(gene, ".rds")))
saveRDS(list(gene = "Failed", error = "Synthetic failure"), file.path(root, "results/Failed.rds"))
generated <- generate_summary_results_slide(project, file.path(root, "results"),
                                            file.path(root, "summary"), gtf_file = "")
res <- generated$res_summary
cs <- generated$res_cs_summary
stopifnot(nrow(res) == 11L, nrow(generated$res_errors) == 4L,
          !any(res$gene %in% c("MissingSlide", "BadOrder", "Failed")),
          is.na(res$perm_cs_slide[res$gene == "MissingPerm"]),
          res$ncs_slide[res$gene == "NoCS"] == 0L,
          all(res$tss_position == 150000),
          res$n_partial_dominant_slide[res$gene == "PartialDominant"] == 1L,
          res$n_partial_recessive_slide[res$gene == "PartialRecessive"] == 1L)
reordered <- cs[cs$gene == "PartialAgreement" & cs$model_key == "fit_slide", ]
stopifnot(all(reordered$component == c(3L, 1L)),
          identical(reordered$lead_delta, c(1, -1)),
          identical(reordered$lead_coding, c("dominant", "recessive")),
          identical(reordered$cs_number, 1:2), !anyDuplicated(reordered$cs_id))
dom_cs <- cs[cs$gene == "Dominant" & cs$model_key == "fit_slide", ]
stopifnot(dom_cs$lead_snp == snps[2], dom_cs$distance_to_tss_kb == 200,
          dom_cs$tss_source == "recovered_from_workhorse_signed_offsets")
annot <- data.frame(gene_name = "Dominant", strand = "-", start = 50000, end = 160000, chromosome = "chr1")
annotated <- slide_summarize_tissue(fixtures$Dominant, "Dominant", "Tissue", "Dominant.rds", annot)
stopifnot(annotated$cs$distance_to_tss_kb[annotated$cs$model_key == "fit_slide"] == -190)
bad_tss <- fixtures$Dominant
bad_tss$susie_mix_lead_snp_tss_distance <- setNames(100, snps[2])
stopifnot(is.na(slide_tss_metadata(bad_tss, "Dominant")$tss_position))
missing_tss <- fixtures$Dominant
missing_tss$susie_add_lead_snp_tss_distance <- NULL
stopifnot(is.na(slide_tss_metadata(missing_tss, "Dominant")$tss_position))
stopifnot(identical(slide_coding(c(-1, -.1, 0, .1, 1)),
                    c("recessive", "partial_recessive", "additive", "partial_dominant", "dominant")))

tables <- run_slide_descriptive_results(project, file.path(root, "summary"), file.path(root, "descriptive"))
audit <- tables$tss_cs_agreement_audit
stopifnot(!any(audit$tss_include[audit$gene == "SameLead"]),
          sum(audit$tss_include[audit$gene == "PartialAgreement"]) == 2L,
          !any(audit$tss_include[audit$gene == "PartialAgreement" & audit$lead_snp == snps[1]]),
          all(audit$tss_include[audit$gene %in% c("NoAddCS", "NoSlideCS")]),
          all(abs(tapply(tables$tss_distance_distribution$proportion,
                         tables$tss_distance_distribution$model, sum) - 1) < 1e-10),
          all(tables$permutation_summary$n_pairs_with_both_permutations == 10L))
dom <- plot_one_cs_slide(project, "dominant", file.path(root, "results"), file.path(root, "summary"),
                         file.path(root, "dominant"))
rec <- plot_one_cs_slide(project, "recessive", file.path(root, "results"), file.path(root, "summary"),
                         file.path(root, "recessive"))
stopifnot(dom$summary$n_compared == 5L, rec$summary$n_compared == 1L,
          dom$comparisons$distance_kb[dom$comparisons$gene == "Dominant"] == 250,
          dom$comparisons$distance_kb[dom$comparisons$gene == "SameLead"] == 0,
          dom$comparisons$additive_lead_selection[dom$comparisons$gene == "NoAddCS"] ==
            "Highest-PIP SNP (no credible set)",
          !any(dom$comparisons$gene == "PartialDominant"))
direction <- slide_one_cs_comparisons(res, file.path(root, "results"), coding_selection = "direction")
stopifnot("PartialDominant" %in% direction$gene,
          direction$lead_delta[direction$gene == "PartialDominant"] == .4)
strict <- slide_one_cs_comparisons(res, file.path(root, "results"), both_one_cs = TRUE)
stopifnot(!any(strict$gene %in% c("MultipleAddCS", "NoAddCS")))
empty <- plot_one_cs_slide(project, "recessive", file.path(root, "results"), file.path(root, "summary"),
                           file.path(root, "empty"), minimum_mean_reads = 1000)
stopifnot(empty$summary$n_candidates == 0L, empty$summary$n_compared == 0L)
run_slide_descriptive_results(project, file.path(root, "summary"), file.path(root, "empty_descriptive"),
                               minimum_mean_reads = 1000)
for (gene in c("SameLead", "NoCS")) {
  d <- file.path(root, gene)
  dir.create(d, showWarnings = FALSE)
  res_summary <- res[res$gene == gene, ]
  res_cs_summary <- cs[cs$gene == gene, ]
  save(res_summary, file = file.path(d, "res_summary.RData"))
  save(res_cs_summary, file = file.path(d, "res_cs_summary.RData"))
  result <- run_slide_descriptive_results(project, d, file.path(d, "descriptive"))
  stopifnot(all(result$tss_distance_distribution$count == 0L),
            all(is.na(result$tss_distance_distribution$proportion)))
}
cat("PASS: slide summary, component-specific coding, permutations, TSS disagreement, distances, and empty plots.\n")
