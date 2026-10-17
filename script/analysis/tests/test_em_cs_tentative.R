# Run from the project root with Rscript --vanilla. No model fitting or GTEx inputs.
source("script/analysis/em_cs_agreement_utils.R")
source("script/analysis/em_cs_comparison_tentative_utils.R")
variant <- function(p) {
  if (!length(p)) return(character())
  paste0("chr1_", sprintf("%.0f", p), "_A_G_b38_A")
}
regions <- sets <- list()
add_case <- function(gene, tissue, add, em, coding = rep("additive", length(em)),
                     add_members = lapply(add, variant), em_members = lapply(em, variant),
                     pips = rep(.9, length(em))) {
  regions[[length(regions) + 1L]] <<- data.frame(gene = gene, tissue = tissue,
    min_pv = 1e-10, mean_count = 200, ncs_susie = length(add), ncs_weighted_fit_mix = length(em),
    has_weighted_fit_mix = TRUE, em_iteration = 18L)
  for (model in c("susie_add", "weighted_fit_mix")) {
    pos <- if (model == "susie_add") add else em
    members <- if (model == "susie_add") add_members else em_members
    if (!length(pos)) next
    sets[[length(sets) + 1L]] <<- data.frame(gene = gene, tissue = tissue, model_key = model,
      cs_number = seq_along(pos), lead_snp = variant(pos),
      lead_coding = if (model == "susie_add") "additive" else coding,
      lead_pip = if (model == "susie_add") .8 else pips,
      cs_snps = vapply(members, paste, collapse = ";", FUN.VALUE = character(1)),
      cs_size_snps = lengths(members), distance_to_tss_kb = (pos - 500000) / 1000, em_iteration = 18L)
  }
}
for (t in 1:5) for (g in 1:20) add_case(paste0("G", t, "_", g), paste0("Tissue", t),
  100000, if (g <= 4 * t) 100000 else c(100000, 350000))
add_case("SharedPlusExtra", "Edges", 100000, c(100000, 350000), c("dominant", "recessive"),
         em_members = list(variant(100000), variant(c(350000, 350001))), pips = c(.99, .4))
add_case("SameSetNewLead", "Edges", 100000, 350000, "dominant",
         list(variant(c(100000, 350000))), list(variant(c(350000, 100000))))
add_case("Partial", "Edges", 100000, 350000, "recessive",
         list(variant(c(100000, 350000, 400000))), list(variant(c(350000, 400000))))
add_case("Disjoint", "Edges", 100000, 350000, "dominant")
add_case("NoAdd", "Edges", numeric(), 350000, "recessive")
add_case("NoEM", "Edges", 100000, numeric())
add_case("NoCS", "Edges", numeric(), numeric())
add_case("DuplicateLead", "Edges", c(100000, 100000), 100000)
add_case("MissingMembership", "Edges", 100000, 100000)
add_case("WrongSize", "Edges", 100000, 100000)
add_case("WrongLead", "Edges", 100000, 100000)
add_case("BadCoordinate", "Edges", 100000, 100000)
add_case("SamePositionOtherAllele", "Edges", 100000, 100000)
add_case("DuplicateEncoding", "Edges", 100000, 100000,
         add_members = list(variant(c(100000, 100000, 350000))))
res <- do.call(rbind, regions); cs <- do.call(rbind, sets)
cs$cs_snps[cs$gene == "MissingMembership" & cs$model_key == "weighted_fit_mix"] <- NA_character_
cs$cs_size_snps[cs$gene == "WrongSize" & cs$model_key == "weighted_fit_mix"] <- 12L
cs$cs_snps[cs$gene == "WrongLead" & cs$model_key == "weighted_fit_mix"] <- variant(900000)
cs$lead_snp[cs$gene == "BadCoordinate"] <- "unparseable"
other <- cs$gene == "SamePositionOtherAllele" & cs$model_key == "weighted_fit_mix"
cs$lead_snp[other] <- cs$cs_snps[other] <- "chr1_100000_C_T_b38_C"
cs$cs_size_snps[cs$gene == "DuplicateEncoding" & cs$model_key == "susie_add"] <- 2L
priors <- data.frame(iteration = 18L, tissue = unique(res$tissue), pi_add = .9, pi_rec = .02, pi_dom = .08)
out <- summarize_em_cs_tentative(res, cs, priors, minimum_tissue_regions = 15)
stopifnot(identical(out$selected_tissues$tissue, c("Tissue1", "Tissue2", "Tissue5", "Tissue4")))
p <- out$matched_cs_pairs; e <- out$unpaired_cs
stopifnot(all(c("SharedPlusExtra", "SameSetNewLead", "Partial", "Disjoint") %in% p$gene),
          p$same_lead_snp[p$gene == "SharedPlusExtra"],
          e$nearest_other_distance_bp[e$gene == "SharedPlusExtra"] == 250000,
          p$jaccard[p$gene == "SameSetNewLead"] == 1,
          p$jaccard[p$gene == "Partial"] == 2/3,
          p$size_ratio_em_to_additive[p$gene == "Partial"] == 2/3,
          p$jaccard[p$gene == "Disjoint"] == 0,
          p$jaccard[p$gene == "DuplicateEncoding"] == .5,
          p$distance_bp[p$gene == "SamePositionOtherAllele"] == 0,
          !p$same_lead_snp[p$gene == "SamePositionOtherAllele"],
          all(is.na(p$jaccard[p$gene %in% c("MissingMembership", "WrongSize", "WrongLead")])),
          !"BadCoordinate" %in% p$gene,
          is.na(e$nearest_other_distance_bp[e$gene == "NoAdd"]),
          sum(p$gene == "DuplicateLead") == 1L, sum(e$gene == "DuplicateLead") == 1L,
          !anyDuplicated(p$additive_cs_id), !anyDuplicated(p$em_cs_id),
          all(out$tissue_ranking$agreement_low95 <= out$tissue_ranking$agreement_rate),
          all(out$tissue_ranking$agreement_high95 >= out$tissue_ranking$agreement_rate))
for (t in out$selected_tissues$tissue) stopifnot(
  sum(out$tissue_cs_count_tables$n_regions[out$tissue_cs_count_tables$tissue == t]) ==
    out$tissue_ranking$n_regions[out$tissue_ranking$tissue == t])
reverse <- summarize_em_cs_tentative(res, cs[nrow(cs):1, ], priors, minimum_tissue_regions = 15)
signature <- function(x) sort(paste(x$additive_cs_id, x$em_cs_id, x$jaccard, x$size_ratio_em_to_additive))
stopifnot(identical(signature(p), signature(reverse$matched_cs_pairs)))
missing_members <- cs; missing_members$cs_snps <- NULL
wrong_iteration <- cs; wrong_iteration$em_iteration[1] <- 17L
stopifnot(inherits(tryCatch(summarize_em_cs_tentative(res, missing_members, priors), error = identity), "error"),
          inherits(tryCatch(summarize_em_cs_tentative(res, wrong_iteration, priors), error = identity), "error"))

output <- Sys.getenv("EM_TENTATIVE_TEST_OUTPUT", unset = tempfile("em-tentative-"))
write_em_cs_tentative(out, output, preview_label = "SYNTHETIC VALIDATION")
empty <- suppressWarnings(summarize_em_cs_tentative(res, cs, priors, association_threshold = 1e-20))
plot_em_cs_tentative(empty, file.path(output, "empty.png"), "SYNTHETIC VALIDATION")
single <- suppressWarnings(summarize_em_cs_tentative(res[res$gene == "SharedPlusExtra", ],
  cs[cs$gene == "SharedPlusExtra", ], priors, minimum_tissue_regions = 1))
plot_em_cs_tentative(single, file.path(output, "single.png"), "SYNTHETIC VALIDATION")
stopifnot(file.info(file.path(output, "em_cs_comparison_tentative_4panel.png"))$size > 10000,
          file.info(file.path(output, "em_cs_comparison_tentative_4panel.pdf"))$size > 1000,
          grDevices::dev.cur() == 1L,
          readRDS(file.path(output, "tentative_plot_data.rds"))$preview_label == "SYNTHETIC VALIDATION")

# Exercise the actual entry script and its provenance-checked summary loader.
project <- file.path(output, "fixture_project")
dir.create(file.path(project, "script/analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project, "script/scan_tissue_attempt"), recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.copy(file.path("script/analysis", c("em_analysis_utils.R", "em_cs_agreement_utils.R",
  "em_cs_comparison_tentative_utils.R")), file.path(project, "script/analysis"), overwrite = TRUE)))
stopifnot(file.copy("script/scan_tissue_attempt/em_utils.R", file.path(project, "script/scan_tissue_attempt"),
                   overwrite = TRUE))
it <- file.path(project, "results_em/iteration_018")
for (part in c("summary", "results", "completed")) dir.create(file.path(it, part), recursive = TRUE, showWarnings = FALSE)
write.csv(priors, file.path(project, "results_em/prior_history.csv"), row.names = FALSE)
write.csv(priors, file.path(it, "priors.csv"), row.names = FALSE)
write.csv(data.frame(chunk = 1L, gene = "Fixture"), file.path(it, "manifest.csv"), row.names = FALSE)
saveRDS(list(iteration = 18L, chunk = 1L, n_genes = 1L), file.path(it, "completed/chunk_001.done"))
saveRDS(list(), file.path(it, "results/Fixture.rds"))
res_summary <- res; res_cs_summary <- cs; res_errors <- cs_errors <- data.frame()
save(res_summary, res_errors, file = file.path(it, "summary/res_summary.RData"))
save(res_cs_summary, cs_errors, file = file.path(it, "summary/res_cs_summary.RData"))
saveRDS(list(iteration = 18L, priors = priors), file.path(it, "summary/em_summary_metadata.rds"))
old_project <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
tryCatch({
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = project)
  env <- new.env(parent = globalenv())
  suppressWarnings(source("script/analysis/plot_cs_comparison_em_tentative.R", local = env))
  stopifnot(env$em_tentative_result$iteration == 18L)
}, finally = {
  if (is.na(old_project)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else Sys.setenv(SUSIE_MIX_PROJECT_DIR = old_project)
})

# Dense layout: several CS counts, overlapping sets, all genetic codings, extras.
regions <- sets <- list(); set.seed(18)
tissues <- c("Adipose Tissue", "Blood Vessel", "Small Intestine", "Salivary Gland", "Thyroid")
for (t in seq_along(tissues)) for (g in 1:130) {
  na <- sample(0:8, 1); ne <- if (runif(1) < c(.4, .5, .85, .9, .7)[t]) na else sample(0:8, 1)
  if (na + ne == 0) ne <- 1L
  a <- sort(sample(1000:999999, na)); b <- sort(sample(1000:999999, ne))
  shared <- if (min(na, ne)) sample(0:min(na, ne), 1) else 0L
  if (shared) b[seq_len(shared)] <- a[seq_len(shared)]
  am <- lapply(a, function(v) variant(unique(c(v, sample(1000:999999, sample(0:40, 1))))))
  bm <- lapply(seq_along(b), function(j) {
    overlap <- if (j <= length(am)) sample(am[[j]], sample(0:length(am[[j]]), 1)) else character()
    unique(c(variant(b[j]), overlap, variant(sample(1000:999999, sample(0:20, 1)))))
  })
  add_case(paste0("Dense", t, "_", g), tissues[t], a, b,
    sample(c("additive", "dominant", "recessive"), ne, replace = TRUE, prob = c(.8, .15, .05)),
    am, bm, runif(ne))
}
res <- do.call(rbind, regions); cs <- do.call(rbind, sets)
priors <- data.frame(iteration = 18L, tissue = tissues, pi_add = .9, pi_rec = .02, pi_dom = .08)
dense <- summarize_em_cs_tentative(res, cs, priors)
stopifnot(all(dense$cs_membership_audit$membership_status == "ok"))
write_em_cs_tentative(dense, file.path(output, "dense"), preview_label = "SYNTHETIC VALIDATION")
cat("PASS: tissue extremes, exact-first pairing, CS overlap/size, invalid membership, empty/single cohorts, entry script, and PNG/PDF exports.\n")
