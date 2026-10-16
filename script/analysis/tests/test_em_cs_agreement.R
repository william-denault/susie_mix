# Rscript --vanilla script/analysis/tests/test_em_cs_agreement.R
# Only base R; no real-data inputs or model fitting.
source("script/analysis/em_cs_agreement_utils.R")
output_dir <- tempfile("em-cs-agreement-")
dir.create(output_dir, recursive = TRUE)

# Check rectangular optimization against all injective assignments, in both
# directions, including duplicates/ties. The short side cannot reuse a lead.
permutations <- function(x, n) {
  if (!n) return(matrix(integer(), nrow = 1L, ncol = 0L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i], n - 1L))))
}
set.seed(17)
for (n in 1:4) for (m in 1:4) for (replicate in 1:8) {
  a <- sample(1:20, n, replace = TRUE); b <- sample(1:20, m, replace = TRUE)
  chosen <- em_rectangular_assignment(a, b)
  short <- if (n <= m) a else b; long <- if (n <= m) b else a
  possibilities <- permutations(seq_along(long), length(short))
  optimum <- min(apply(possibilities, 1, function(p) sum(abs(short - long[p]))))
  stopifnot(nrow(chosen) == min(n, m), !anyDuplicated(chosen[, 1]), !anyDuplicated(chosen[, 2]),
            sum(abs(a[chosen[, 1]] - b[chosen[, 2]])) == optimum)
}
stopifnot(nrow(em_rectangular_assignment(numeric(), 100)) == 0L)

variant <- function(pos) paste0("chr1_", sprintf("%.0f", pos), "_A_G_b38_A")
region_rows <- cs_rows <- list()
add_region <- function(gene, tissue, add, em, codings = rep("dominant", length(em)),
                       pip = rep(.8, length(em)), min_pv = 1e-10, mean_count = 200) {
  region_rows[[length(region_rows) + 1L]] <<- data.frame(gene = gene, tissue = tissue,
    ncs_susie = length(add), ncs_weighted_fit_mix = length(em), has_weighted_fit_mix = TRUE,
    min_pv = min_pv, mean_count = mean_count, em_iteration = 18L)
  for (model in c("susie_add", "weighted_fit_mix")) {
    positions <- if (model == "susie_add") add else em
    if (!length(positions)) next
    cs_rows[[length(cs_rows) + 1L]] <<- data.frame(gene = gene, tissue = tissue,
      model_key = model, cs_number = seq_along(positions), lead_snp = variant(positions),
      lead_coding = if (model == "susie_add") "additive" else codings,
      lead_pip = if (model == "susie_add") .9 else pip,
      distance_to_tss_kb = (positions - 200000) / 1000, em_iteration = 18L)
  }
}
# Five tissues, equal denominators, with unambiguous disagreement ranks.
for (t in 1:5) for (g in 1:25) {
  add_region(paste0("Rank", t, "_", g), paste0("Tissue", t), 100000,
             if (g <= 5L * t) 100000 else c(100000, 300000))
}
# Small sample with complete disagreement must not outrank eligible tissues.
add_region("Tiny", "Tiny tissue", 100000, c(100000, 500000))
# Explicit 1-vs-2 example: retain exact lead and report the other as extra.
add_region("Example", "Edge cases", 100000, c(100000, 350000), c("additive", "recessive"), c(.95, .7))
add_region("ExactFirst", "Edge cases", c(100, 200), c(200, 300))
add_region("NearestTrap", "Edge cases", c(12, 112), c(5, 17, 500))
add_region("ExtraAdd", "Edge cases", c(100000, 350000), 100000)
add_region("EMOnly", "Edge cases", numeric(), 200000, "dominant", NA_real_)
add_region("AddOnly", "Edge cases", 200000, numeric())
add_region("BothZero", "Edge cases", numeric(), numeric())
add_region("DuplicateLead", "Edge cases", c(100000, 100000), 100000)
add_region("BadCoordinates", "Edge cases", 100000, 100000)
add_region("Incomplete", "Edge cases", 100000, c(100000, 200000))
add_region("FilteredP", "Edge cases", 100000, 100000, min_pv = 1e-4)
add_region("FilteredReads", "Edge cases", 100000, 100000, mean_count = 20)
add_region("SamePosition", "Edge cases", 100000, 100000)
add_region("ExtraAtTSS", "Edge cases", 100000, c(100000, 200000))
res <- do.call(rbind, region_rows); cs <- do.call(rbind, cs_rows)
cs$lead_snp[cs$gene == "BadCoordinates"] <- "invalid"
cs <- cs[!(cs$gene == "Incomplete" & cs$model_key == "weighted_fit_mix" & cs$cs_number == 2L), ]
cs$lead_snp[cs$gene == "SamePosition" & cs$model_key == "weighted_fit_mix"] <- "chr1_100000_C_T_b38_C"
cs$distance_to_tss_kb[cs$gene == "Example" & cs$model_key == "weighted_fit_mix" & cs$cs_number == 2L] <- -150
cs$distance_to_tss_kb[cs$gene == "ExtraAdd" & cs$model_key == "susie_add" & cs$cs_number == 2L] <- NA_real_
priors <- data.frame(iteration = 18L, source_iteration = 17L, tissue = sort(unique(res$tissue)),
  pi_add = .8, pi_rec = .05, pi_dom = .15, prior_status = "estimated")
result <- summarize_em_cs_agreement(res, cs, minimum_tissue_regions = 15L, em_priors = priors)
stopifnot(identical(result$selected_tissues$tissue, c("Tissue1", "Tissue2", "Tissue3", "Tissue4")))
audit <- result$region_audit
status <- setNames(audit$matching_status, audit$gene)
stopifnot(status["Incomplete"] == "incomplete_cs_summary",
          status["BadCoordinates"] == "invalid_or_inconsistent_coordinates",
          audit$lead_status[audit$gene == "BadCoordinates"] == "ok",
          all(!audit$count_eligible[audit$gene %in% c("BothZero", "FilteredP", "FilteredReads")]))
p <- result$matched_cs_pairs; u <- result$unpaired_cs
example <- p[p$gene == "Example", ]; extra <- u[u$gene == "Example", ]
stopifnot(nrow(example) == 1L, example$same_lead_snp, example$distance_bp == 0,
          nrow(extra) == 1L, extra$nearest_other_distance_bp == 250000,
          extra$lead_coding == "recessive", extra$model == "EM SuSiE-mix")
exact <- p[p$gene == "ExactFirst", ]
stopifnot(sum(exact$same_lead_snp) == 1L, identical(sort(exact$distance_bp), c(0, 200)),
          p$distance_bp[p$gene == "SamePosition"] == 0, !p$same_lead_snp[p$gene == "SamePosition"],
          sum(p$gene == "DuplicateLead") == 1L, sum(u$gene == "DuplicateLead") == 1L,
          is.na(u$nearest_other_distance_bp[u$gene == "EMOnly"]),
          !u$other_model_has_cs[u$gene == "EMOnly"],
          !any(p$gene %in% c("Incomplete", "BadCoordinates")),
          result$extra_em_coding_summary$n_without_additive_cs[result$extra_em_coding_summary$coding == "dominant"] == 1L)
stopifnot(sum(result$lead_retention_summary$n_regions) == sum(audit$count_eligible & audit$lead_status == "ok"))
# C uses the exact extra-CS populations in B; its distance is to the gene TSS,
# not the other lead. Missing TSSs are audited, negative distances become
# absolute distances, zero remains valid, and model-only regions stay out.
stopifnot(extra$distance_to_tss_kb == -150, extra$absolute_distance_to_tss_bp == 150000,
          all(u$tss_panel_status[u$gene %in% c("EMOnly", "AddOnly")] == "no_counterpart_in_panel_b"),
          u$tss_panel_status[u$gene == "ExtraAdd"] == "missing_tss_distance",
          u$absolute_distance_to_tss_bp[u$gene == "ExtraAtTSS"] == 0,
          all(result$tss_distance_summary$n_panel_b == result$lead_distance_summary$n_distance[2:3]),
          sum(result$tss_distance_summary$n_missing_tss) == 1L,
          0 %in% result$tss_distances[["Extra EM"]],
          nrow(result$tissue_coding_priors) == nrow(priors),
          all(result$tissue_coding_priors$pi_rec == .05),
          all(result$tissue_coding_priors$pi_dom == .15),
          all(abs(result$tissue_coding_priors$pi_nonadditive - .2) < 1e-10))
for (tissue in result$selected_tissues$tissue) {
  tab <- result$tissue_cs_count_tables[result$tissue_cs_count_tables$tissue == tissue, ]
  t <- result$tissue_ranking[result$tissue_ranking$tissue == tissue, ]
  stopifnot(sum(tab$n_regions) == t$n_regions,
            sum(tab$n_regions[tab$additive_cs_count == tab$em_cs_count]) == t$same_cs_count)
}
# CS row order cannot change exact matches, extras, or distances.
reordered <- summarize_em_cs_agreement(res, cs[rev(seq_len(nrow(cs))), ], minimum_tissue_regions = 15L, em_priors = priors)
pair_signature <- function(x) sort(paste(x$additive_cs_id, x$em_cs_id, x$distance_bp))
stopifnot(identical(pair_signature(result$matched_cs_pairs), pair_signature(reordered$matched_cs_pairs)),
          setequal(result$unpaired_cs$cs_id, reordered$unpaired_cs$cs_id))
write_em_cs_agreement(result, output_dir, show_plot = FALSE, preview_label = "SYNTHETIC VALIDATION")
stopifnot(file.info(file.path(output_dir, "em_cs_agreement_4panel.png"))$size > 10000,
          file.info(file.path(output_dir, "em_cs_agreement_4panel.pdf"))$size > 1000,
          grDevices::dev.cur() == 1L)
# Empty selection, a single region, all-zero distances, missing PIPs.
empty <- suppressWarnings(summarize_em_cs_agreement(res, cs, association_threshold = 1e-20, minimum_tissue_regions = 1L, em_priors = priors))
plot_em_cs_agreement(empty, file.path(output_dir, "empty.png"), "SYNTHETIC VALIDATION")
one <- suppressWarnings(summarize_em_cs_agreement(res[res$gene == "Example", ], cs[cs$gene == "Example", ], minimum_tissue_regions = 1L, em_priors = priors))
plot_em_cs_agreement(one, file.path(output_dir, "single.png"), "SYNTHETIC VALIDATION")
bad_iteration <- cs; bad_iteration$em_iteration[1] <- 17L
stopifnot(inherits(tryCatch(summarize_em_cs_agreement(res, bad_iteration, em_priors = priors), error = identity), "error"))
bad_priors <- priors; bad_priors$iteration <- 17L
stopifnot(inherits(tryCatch(summarize_em_cs_agreement(res, cs, em_priors = bad_priors), error = identity), "error"))
bad_priors <- priors; bad_priors$pi_add <- .9
stopifnot(inherits(tryCatch(summarize_em_cs_agreement(res, cs, em_priors = bad_priors), error = identity), "error"))
bad_priors <- priors[-1L, ]
stopifnot(inherits(tryCatch(summarize_em_cs_agreement(res, cs, em_priors = bad_priors), error = identity), "error"))
missing_tss <- cs; missing_tss$distance_to_tss_kb <- NULL
stopifnot(inherits(tryCatch(summarize_em_cs_agreement(res, missing_tss, em_priors = priors), error = identity), "error"))

# Exercise the actual RStudio Source entry with the existing context/summary
# loaders. Build a minimal complete EM project with matching provenance.
project <- file.path(output_dir, "project")
dir.create(file.path(project, "script/analysis"), recursive = TRUE)
dir.create(file.path(project, "script/scan_tissue_attempt"), recursive = TRUE)
file.copy(file.path("script/analysis", c("em_cs_agreement_utils.R", "em_analysis_utils.R")), file.path(project, "script/analysis"))
file.copy("script/scan_tissue_attempt/em_utils.R", file.path(project, "script/scan_tissue_attempt"))
iteration_dir <- file.path(project, "results_em/iteration_018")
dir.create(file.path(iteration_dir, "summary"), recursive = TRUE)
dir.create(file.path(iteration_dir, "results"))
dir.create(file.path(iteration_dir, "completed"))
write.csv(priors, file.path(project, "results_em/prior_history.csv"), row.names = FALSE)
write.csv(priors, file.path(iteration_dir, "priors.csv"), row.names = FALSE)
write.csv(data.frame(chunk = 1L, gene = "Fixture"), file.path(iteration_dir, "manifest.csv"), row.names = FALSE)
saveRDS(list(iteration = 18L, chunk = 1L, n_genes = 1L), file.path(iteration_dir, "completed/chunk_001.done"))
saveRDS(list(), file.path(iteration_dir, "results/Fixture.rds"))
res_summary <- res; res_cs_summary <- cs
res_errors <- cs_errors <- data.frame()
save(res_summary, res_errors, file = file.path(iteration_dir, "summary/res_summary.RData"))
save(res_cs_summary, cs_errors, file = file.path(iteration_dir, "summary/res_cs_summary.RData"))
saveRDS(list(iteration = 18L, priors = priors), file.path(iteration_dir, "summary/em_summary_metadata.rds"))
old_project <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
tryCatch({
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = project)
  env <- new.env(parent = globalenv())
  suppressWarnings(source("script/analysis/plot_cs_agreement_em.R", local = env))
  stopifnot(env$em_agreement_results$iteration == 18L,
            file.exists(file.path(iteration_dir, "descriptive_results/cs_agreement/em_cs_agreement_4panel.png")))
}, finally = {
  if (is.na(old_project)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else Sys.setenv(SUSIE_MIX_PROJECT_DIR = old_project)
})
cat("PASS: assignments, four disagreement ranks, exact extra-lead TSS groups, prior validation, audits, empty selections, Source entry, and PNG/PDF exports.\n")
# Optional persistent preview for visual QA; the label prevents confusion with real results.
preview_dir <- Sys.getenv("EM_CS_AGREEMENT_TEST_OUTPUT", unset = "")
if (nzchar(preview_dir)) write_em_cs_agreement(result, preview_dir, FALSE, "SYNTHETIC VALIDATION")

# Publication-layout stress case: broad CS counts, long tissue names, many
# changed leads, and varied coding-specific PIPs. Also tests all plot groups.
region_rows <- cs_rows <- list()
set.seed(57)
tissue_names <- c("Adipose Tissue", "Blood Vessel", "Small Intestine", "Salivary Gland")
for (t in seq_along(tissue_names)) for (g in 1:110) {
  nadd <- sample(0:10, 1)
  nem <- if (runif(1) < c(.3, .4, .8, .9)[t]) nadd else sample(0:10, 1)
  if (nadd + nem == 0) nem <- 1L
  a <- sort(sample(100000:900000, nadd))
  m <- sort(sample(100000:900000, nem))
  if (min(nadd, nem) && runif(1) < .8) {
    shared <- sample(seq_len(min(nadd, nem)), 1)
    m[seq_len(shared)] <- a[seq_len(shared)]
  }
  add_region(paste0("Dense", t, "_", g), tissue_names[t], a, m,
              sample(c("additive", "recessive", "dominant"), nem, replace = TRUE), runif(nem))
}
dense_priors <- data.frame(iteration = 18L, source_iteration = 17L,
  tissue = c(tissue_names, paste("Other tissue", 1:23)),
  pi_rec = seq(.001, .02, length.out = 27), pi_dom = seq(.15, .02, length.out = 27))
dense_priors$pi_add <- 1 - dense_priors$pi_rec - dense_priors$pi_dom
dense <- summarize_em_cs_agreement(do.call(rbind, region_rows), do.call(rbind, cs_rows), em_priors = dense_priors)
stopifnot(nrow(dense$selected_tissues) == 4L, nrow(dense$matched_cs_pairs) > 100L,
          all(dense$lead_distance_summary$n_distance > 0L))
dense_path <- file.path(if (nzchar(preview_dir)) preview_dir else output_dir, "synthetic_dense_4panel.png")
plot_em_cs_agreement(dense, dense_path, "SYNTHETIC VALIDATION")
stopifnot(grDevices::dev.cur() == 1L)
cat("PASS: dense multi-CS layout, all distance groups, and 27 tissue priors.\n")
