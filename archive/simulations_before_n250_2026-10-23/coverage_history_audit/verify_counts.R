out <- "tmp/coverage_history_audit"
checks <- list()
for (archive in c("old", "current")) {
  base <- if (archive == "old") "simulation results" else "simulation results/slide_prior_v1"
  summary <- read.csv(file.path(base, "figures/metric_summary.csv"))
  summary <- summary[summary$scenario == "Additive only" & summary$method %in% c("SuSiE", "SuSiE-mix"), ]
  cells <- unique(summary[c("pve", "K")])
  for (i in seq_len(nrow(cells))) {
    k <- cells$K[i]; pve <- cells$pve[i]
    pattern <- paste0("^additive_add", k, "_rec0_dom0", if (archive != "old") "_prec0_pdom0" else "",
                      "_n500_L10_pve", format(pve, scientific = FALSE, trim = TRUE), "_seed.*_reps400_chunk1\\.RData$")
    files <- list.files(file.path(base, "chunks"), pattern, full.names = TRUE)
    stopifnot(length(files) == 1L)
    e <- new.env(parent = emptyenv()); load(files, e)
    successful <- Filter(function(x) is.null(x$error), e$results)
    stopifnot(!anyDuplicated(vapply(successful, function(x) x$seed, numeric(1))))
    for (method in c("SuSiE", "SuSiE-mix")) {
      counts <- t(vapply(successful, function(x) {
        sets <- if (method == "SuSiE") x$susie_cs$cs else x$cs_mix_as_additive_indices
        c(n_cs = length(sets), covered_cs = sum(vapply(sets, function(s) any(s %in% x$true_pos), logical(1))))
      }, numeric(2)))
      totals <- colSums(counts)
      row <- summary[summary$pve == pve & summary$K == k & summary$method == method, ]
      stopifnot(nrow(row) == 1L, row$n_cs == totals["n_cs"], row$covered_cs == totals["covered_cs"],
                row$n_replicates == length(successful), abs(row$coverage - totals["covered_cs"] / totals["n_cs"]) < 1e-12)
      checks[[length(checks) + 1L]] <- data.frame(archive, pve, K = k, method,
        n_replicates = length(successful), n_cs = totals["n_cs"], covered_cs = totals["covered_cs"],
        coverage = totals["covered_cs"] / totals["n_cs"], verified = TRUE, row.names = NULL)
    }
  }
  cat(archive, ":", nrow(cells), "additive cells checked against saved credible sets\n")
}
verified <- do.call(rbind, checks)
write.csv(verified, file.path(out, "verified_counts.csv"), row.names = FALSE)
wide <- reshape(verified[c("archive", "pve", "K", "method", "coverage")],
                idvar = c("archive", "pve", "K"), timevar = "method", direction = "wide")
wide$mix_minus_susie_pp <- 100 * (wide[["coverage.SuSiE-mix"]] - wide$coverage.SuSiE)
old <- wide[wide$archive == "old", -1L]
current <- wide[wide$archive == "current", -1L]
comparison <- merge(old, current, by = c("pve", "K"), all = TRUE, suffixes = c("_old", "_current"))
write.csv(comparison, file.path(out, "coverage_comparison_wide.csv"), row.names = FALSE)
common <- current[current$pve >= .1, ]
cat("Current comparable cells: mix higher =", sum(common$mix_minus_susie_pp > 1e-10),
    "; tie =", sum(abs(common$mix_minus_susie_pp) < 1e-10), "; lower =", sum(common$mix_minus_susie_pp < -1e-10), "\n")
cat("Publication snapshot exactly matches current metric summary:",
    identical(read.csv("publication_plots/input_data/metric_summary.csv"),
              read.csv("simulation results/slide_prior_v1/figures/metric_summary.csv")), "\n")
cat("All", nrow(verified), "method/cell summaries verified; no fits rerun.\n")
