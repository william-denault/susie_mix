# Run from the repository root with base R:
# Rscript --vanilla script/analysis/tests/test_matched_cs_lead_distances.R
source("script/analysis/plot_matched_cs_lead_distances.R")

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i]))))
}

# Independent exhaustive oracle: compare both the total absolute distance
# and the secondary squared-distance criterion against every assignment.
set.seed(251)
for (k in 1:6) {
  perms <- permutations(seq_len(k))
  for (trial in 1:35) {
    a <- sample(1:25, k, replace = TRUE)
    m <- sample(1:25, k, replace = TRUE)
    costs <- t(apply(perms, 1, function(p) {
      d <- abs(a - m[p])
      c(total = sum(d), squared = sum(d^2))
    }))
    result <- minimum_cs_distance_assignment(a, m)
    best <- min(costs[, "total"])
    stopifnot(result$minimum_total_distance_bp == best,
              sum(result$distance_bp^2) == min(costs[costs[, "total"] == best, "squared"]),
              !anyDuplicated(result$add_index), !anyDuplicated(result$mix_index))
    # A reported tie must really have an alternative minimum-total assignment.
    if (result$minimum_total_has_tie) stopifnot(sum(costs[, "total"] == best) > 1L)
  }
}

# A greedy nearest-available approach gives 5 + 107 = 112 bp here; optimum 102.
greedy_trap <- minimum_cs_distance_assignment(c(12, 112), c(5, 17))
stopifnot(identical(greedy_trap$distance_bp, c(7, 95)))
# Independent nearest neighbours would reuse lead 200; sorted assignment uses
# it once and resolves the (0,200) versus (100,100) minimum-total tie explicitly.
tied <- minimum_cs_distance_assignment(c(100, 200), c(200, 300))
stopifnot(identical(tied$distance_bp, c(100, 100)), tied$minimum_total_has_tie)

make_cs <- function(gene, tissue, add, mix, coding = "dominant") {
  do.call(rbind, lapply(c("susie_add", "susie_mix"), function(model) {
    pos <- if (model == "susie_add") add else mix
    n <- length(pos)
    data.frame(
      gene = gene, tissue = tissue, result_file = paste0(gene, ".rds"), model_key = model,
      cs_number = seq_len(n), cs_name = paste0("L", seq_len(n) + 3L), component = seq_len(n) + 3L,
      lead_snp = paste0("chr1_", sprintf("%.0f", pos), "_A_G_b38_A"), lead_position = pos,
      lead_coding = if (model == "susie_add") rep("additive", n) else rep(coding, length.out = n),
      lead_pip = rep(c(.9, .3, .7), length.out = n), lead_alpha = rep(.8, n),
      stringsAsFactors = FALSE
    )
  }))
}
res_summary <- data.frame(
  gene = c("THREE", "ONE", "SAME", "LOCUS", "MISSINGPIP", "INCOMPLETE", "DUPCS", "BADBP", "CHROM", "BUILD", "STALE", "ZERO", "UNEQUAL", "BADCOUNT", "NOFILTER", "TIE"),
  tissue = c("Blood", "Liver", "Blood", "Blood", "Liver", rep("Blood", 11)),
  ncs_susie = c(3, 1, 1, 2, 1, 2, 2, 1, 1, 1, 1, 0, 1, NA, 1, 2),
  ncs_susie_mix = c(3, 1, 1, 2, 1, 2, 2, 1, 1, 1, 1, 0, 2, 1, 1, 2),
  min_pv = 1, mean_count = 1, converged_add = FALSE, converged_mix = FALSE,
  stringsAsFactors = FALSE
)
res_summary$result_file <- paste0(res_summary$gene, ".rds")
res_cs_summary <- rbind(
  make_cs("THREE", "Blood", c(300000, 100000, 200000), c(201000, 310000, 100000), c("recessive", "dominant", "additive")),
  make_cs("ONE", "Liver", 100, 101),
  make_cs("SAME", "Blood", 1000, 1000),
  make_cs("LOCUS", "Blood", c(1000, 1000), c(1000, 1000)),
  make_cs("MISSINGPIP", "Liver", 100, 1000000),
  make_cs("INCOMPLETE", "Blood", c(100, 200), 100),
  make_cs("DUPCS", "Blood", c(100, 200), c(100, 200)),
  make_cs("BADBP", "Blood", 100, 0),
  make_cs("CHROM", "Blood", 100, 200),
  make_cs("BUILD", "Blood", 100, 200),
  make_cs("STALE", "Blood", 100, 200),
  make_cs("NOFILTER", "Blood", 500, 600),
  make_cs("TIE", "Blood", c(100, 200), c(200, 300))
)
is_mix <- res_cs_summary$model_key == "susie_mix"
res_cs_summary$lead_pip[res_cs_summary$gene == "MISSINGPIP" & is_mix] <- NA_real_
res_cs_summary$cs_number[res_cs_summary$gene == "DUPCS" & is_mix] <- 1L
res_cs_summary$lead_snp[res_cs_summary$gene == "CHROM" & is_mix] <- "chr2_200_A_G_b38_A"
res_cs_summary$lead_snp[res_cs_summary$gene == "BUILD" & is_mix] <- "chr1_200_A_G_b37_A"
res_cs_summary$result_file[res_cs_summary$gene == "STALE" & is_mix] <- "different_fit.rds"
# Positional metadata from a different workhorse SNP must not override lead_snp.
res_cs_summary$lead_position[res_cs_summary$gene == "ONE" & is_mix] <- 9999
# Same base position, different variant identity: zero_distance != same_lead_snp.
loc <- which(res_cs_summary$gene == "LOCUS" & is_mix)[1]
res_cs_summary$lead_snp[loc] <- "chr1_1000_C_T_b38_C"
weighted <- make_cs("THREE", "Blood", 1, 1)[1, ]
weighted$model_key <- "weighted_fit_mix"
orphan <- make_cs("ORPHAN", "Blood", 100, 200)
res_cs_summary <- rbind(res_cs_summary, weighted, orphan)

result <- match_equal_cs_regions(res_summary, res_cs_summary)
pairs <- result$pairs
audit <- result$region_audit
status <- setNames(audit$status, audit$gene)
stopifnot(
  nrow(pairs) == 11L, nrow(result$orphan_cs) == 2L,
  identical(pairs$distance_bp[pairs$gene == "THREE"], c(0, 1000, 10000)),
  nrow(pairs[pairs$gene == "THREE", ]) == 3L,
  all(tapply(pairs$region_weight, pairs$region_id, sum) == 1),
  status["INCOMPLETE"] == "incomplete_or_inconsistent_cs_summary",
  status["DUPCS"] == "invalid_or_duplicate_cs_number",
  status["BADBP"] == "invalid_lead_coordinates",
  status["BUILD"] == "invalid_lead_coordinates",
  status["CHROM"] == "inconsistent_chromosomes",
  status["STALE"] == "inconsistent_result_files",
  status["ZERO"] == "equal_zero_cs", status["UNEQUAL"] == "unequal_cs_counts",
  status["BADCOUNT"] == "invalid_cs_count", status["NOFILTER"] == "matched",
  pairs$distance_bp[pairs$gene == "ONE"] == 1,
  pairs$distance_bp[pairs$gene == "MISSINGPIP"] == 999900,
  !pairs$scatter_eligible[pairs$gene == "MISSINGPIP"],
  sum(pairs$zero_distance) == 4L, sum(pairs$same_lead_snp) == 3L,
  all(pairs$minimum_total_has_tie[pairs$gene == "TIE"]),
  !any(pairs$converged_mix), all(pairs$mean_count == 1)
)
shuffled <- match_equal_cs_regions(res_summary[sample(nrow(res_summary)), ], res_cs_summary[sample(nrow(res_cs_summary)), ])
stopifnot(identical(pairs, shuffled$pairs))
dup_region <- try(match_equal_cs_regions(rbind(res_summary, res_summary[1, ]), res_cs_summary), silent = TRUE)
stopifnot(inherits(dup_region, "try-error"))
summary <- matched_cs_distance_summary(pairs)
stopifnot(summary$n_pairs == 11, summary$n_regions == 7, summary$n_zero_distance == 4,
          summary$n_missing_or_invalid_mixed_pip == 1, summary$median_distance_bp == 100)
empty <- match_equal_cs_regions(res_summary[res_summary$gene == "ZERO", ], data.frame())
stopifnot(nrow(empty$pairs) == 0, empty$region_audit$status == "equal_zero_cs")

out <- "tmp/matched_cs_distance_validation"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save(res_summary, file = file.path(out, "res_summary.RData"))
save(res_cs_summary, file = file.path(out, "res_cs_summary.RData"))
write.csv(res_summary, file.path(out, "res_summary.csv"), row.names = FALSE)
write.csv(res_cs_summary, file.path(out, "res_cs_summary.csv"), row.names = FALSE)
run <- suppressWarnings(run_matched_cs_lead_distances(c(
  paste0("--summary-file=", out, "/res_summary.RData"),
  paste0("--cs-summary-file=", out, "/res_cs_summary.RData"),
  paste0("--output-dir=", out, "/all")
)))
run_csv <- suppressWarnings(run_matched_cs_lead_distances(c(
  paste0("--summary-file=", out, "/res_summary.csv"),
  paste0("--cs-summary-file=", out, "/res_cs_summary.csv"),
  paste0("--output-dir=", out, "/csv")
)))
stopifnot(isTRUE(all.equal(run$pairs, run_csv$pairs)))
blood <- run_matched_cs_lead_distances(c(
  paste0("--pairs-file=", out, "/all/matched_cs_pairs.csv"), "--tissue=Blood"
))
blood_summary <- read.csv(file.path(blood$output_dir, "distance_summary.csv"))
stopifnot(blood_summary$n_pairs == 9, blood_summary$n_regions == 5)
plot_count <- plot_matched_cs_distances(pairs, file.path(out, "fixture_scatter.png"))
stopifnot(plot_count$n_zero == 4, plot_count$n_positive == 6, plot_count$n_missing_pip == 1)
plot_matched_cs_distances(empty$pairs, file.path(out, "empty_scatter.png"))
zero_count <- plot_matched_cs_distances(pairs[pairs$distance_bp == 0, ], file.path(out, "zero_scatter.png"))
stopifnot(zero_count$n_zero == 4, zero_count$n_positive == 0)
plot_matched_cs_distances(pairs[pairs$gene == "ONE", ], file.path(out, "one_scatter.png"))

# Dense synthetic preview exercises overplotting/legends without presenting it
# as real biological output. The full summaries are not stored in this checkout.
dense <- pairs[rep(which(pairs$gene == "THREE"), 500), ]
dense$gene <- paste0("Synthetic", seq_len(nrow(dense)))
dense$region_id <- matched_cs_key(dense$gene, dense$tissue)
dense$distance_bp <- c(rep(0, 400), round(10^runif(nrow(dense) - 400, 0, 6)))
dense$distance_kb <- dense$distance_bp / 1000
dense$mixed_lead_pip <- runif(nrow(dense))^3
dense$n_cs[seq_len(400)] <- 1L
plot_matched_cs_distances(dense, file.path(out, "synthetic_dense_scatter.png"))
cat("All matched-CS tests passed, including exhaustive assignment checks and end-to-end RData/CSV/tissue runs.\n")
