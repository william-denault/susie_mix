# Run from the project root with Rscript --vanilla.
options(susie_mix.lead_pip_plots.run = FALSE)
source("script/sim/plot_lead_pip_coding.R")
equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))
toy <- function(pip = c(.8, .1, 0, .1, .2, .1, .2, 0)) {
  list(settings = list(L_add = 1, L_rec = 0, L_dom = 0, all_additive = FALSE,
    n = 500, L = 10, pve = .2, min_maf = .05, hwe_thresh = 1e-8, min_n_rec = 5),
    seed = 1, true_pos = 1L, true_pos_mix = 1L, causal_coding = "additive",
    susie_pip = c(.9, .1, .2), susie_mix_pip = pip, susie_mix_pip_snp = c(.9, .3, .4),
    mix_to_add = c(1L, 2L, 3L, 1L, 3L, 1L, 2L, 3L),
    metrics = data.frame(method = c("SuSiE", "SuSiE-mix"), converged = TRUE))
}
outcome <- function(x) lp_outcomes[lp_assess_replicate(x)$metrics["outcome_id"]]
x <- toy()
equal(lp_assess_replicate(x)$metrics[c("lead_index", "lead_snp_index", "lead_coding_id")], c(1, 1, 1))
stopifnot(outcome(x) == "exact")
wrong_code <- x; wrong_code$susie_mix_pip[6] <- .9
stopifnot(outcome(wrong_code) == "wrong_coding")
wrong_snp <- x; wrong_snp$susie_mix_pip[7] <- .95
stopifnot(outcome(wrong_snp) == "wrong_snp")
# A lead at an unrelated SNP is selected even if it lies outside any credible set.
wrong_snp$susie_mix_cs <- list(cs = list(L1 = 1L))
equal(lp_assess_replicate(wrong_snp)$metrics["lead_index"], 7)
tie <- x; tie$susie_mix_pip[6] <- .8
stopifnot(outcome(tie) == "ambiguous", is.na(lp_assess_replicate(tie)$metrics["lead_index"]))
equal(lp_assess_replicate(tie)$metrics[c("n_top", "n_top_exact")], c(2, 1))
near <- tie; near$susie_mix_pip[6] <- .8 - 5e-13
stopifnot(outcome(near) == "ambiguous")
near$susie_mix_pip[6] <- .8 - 5e-11
stopifnot(outcome(near) == "exact")
zero <- x; zero$susie_mix_pip[] <- 0
stopifnot(outcome(zero) == "no_support")
equal(lp_assess_replicate(zero)$metrics["n_top"], 0)
low <- x; low$susie_mix_pip <- low$susie_mix_pip / 1000
stopifnot(outcome(low) == "exact") # Ranking does not impose a detection cutoff.
multiple <- x; multiple$true_pos <- c(1L, 3L); multiple$true_pos_mix <- c(1L, 5L)
multiple$causal_coding <- c("additive", "recessive"); multiple$susie_mix_pip[5] <- .99
stopifnot(outcome(multiple) == "exact")
equal(lp_assess_replicate(multiple)$metrics["lead_index"], 5)
rounded <- x; rounded$susie_mix_pip[1] <- 1 + .Machine$double.eps
equal(lp_assess_replicate(rounded)$metrics["lead_pip"], 1)
stopifnot(lp_assess_replicate(rounded)$n_clamped == 1)
bad <- x; bad$susie_mix_pip[1] <- 1.01
stopifnot(inherits(tryCatch(lp_assess_replicate(bad), error = identity), "error"))

# The integrated reader must preserve the sampling unit, duplicate policy,
# true all-additive controls, and separate denominators for tied/unique leads.
dir.create("tmp", showWarnings = FALSE)
scratch <- normalizePath("tmp", winslash = "/")
test_dir <- tempfile("lead_pip_test_", tmpdir = scratch)
dir.create(file.path(test_dir, "chunks"), recursive = TRUE)
save_checkpoint <- function(results, prefix, requested) {
  save(results, file = file.path(test_dir, "chunks",
    paste0(prefix, "_n500_L10_pve0.2_seed100_reps", requested, "_chunk1.RData")))
}
cases <- list(x, wrong_code, wrong_snp, tie, zero, rounded)
for (i in seq_along(cases)) cases[[i]]$seed <- i
save_checkpoint(c(cases, list(list(seed = 7, error = "Failed simulation")), list(bad)), "scenario_add1_rec0_dom0", 8)
save_checkpoint(cases[1:2], "scenario_add1_rec0_dom0", 2)
control <- wrong_code
control$settings$L_add <- 0; control$settings$L_rec <- 1; control$settings$all_additive <- TRUE
save_checkpoint(list(control), "control_add0_rec1_dom0", 1)
config <- lead_pip_settings
config$chunk_dir <- file.path(test_dir, "chunks")
config$output_dir <- file.path(test_dir, "figures")
config$pve_values <- .2
analysis <- run_lead_pip_plots(config)
stopifnot(nrow(analysis$replicates) == 7, sum(analysis$audit$duplicates) == 2,
  sum(analysis$audit$errors) == 1, sum(analysis$audit$invalid) == 1,
  all(analysis$replicates$scenario == "Additive only"),
  length(list.files(config$output_dir, "[.]pdf$")) == 4,
  length(list.files(config$output_dir, "[.]png$")) == 4)
summary <- analysis$summary[analysis$summary$scope == "All simulations", ]
equal(summary[c("n_runs", "n_unique", "n_exact", "n_wrong_coding", "n_wrong_snp", "n_ambiguous", "n_no_support")],
  data.frame(n_runs = 7, n_unique = 5, n_exact = 2, n_wrong_coding = 2, n_wrong_snp = 1, n_ambiguous = 1, n_no_support = 1))
equal(summary$exact_accuracy_unique, 2 / 5)
equal(summary$fraction_exact, 2 / 7)
equal(summary$share_unique_dominant, 3 / 5)
stopifnot(all(abs(rowSums(analysis$bins[paste0("fraction_", lp_outcomes)]) - 1) < 1e-12))
bins <- analysis$bins[analysis$bins$scope == "All simulations", ]
equal(bins$n_no_support[bins$bin == 1], 1)
equal(bins$n_exact[bins$bin == 10], 1) # PIP=1 stays in the last bin.
equal(sum(analysis$accuracy$n_leads[analysis$accuracy$scope == "All simulations"]), 5)
equal(sum(analysis$accuracy$n_exact[analysis$accuracy$scope == "All simulations"]), 2)

# No unique leads means undefined accuracy/shares, not perfect performance.
ambiguous <- analysis$replicates[analysis$replicates$outcome %in% c("ambiguous", "no_support"), ]
empty <- lp_summarize(ambiguous, config)
stopifnot(nrow(empty$accuracy) == 0, all(is.na(empty$summary$exact_accuracy_unique)),
  all(is.na(empty$additive$share_unique_dominant)))
lp_render_plots(c(list(replicates = ambiguous), empty), config)

resolved <- normalizePath(test_dir, winslash = "/", mustWork = TRUE)
stopifnot(startsWith(resolved, paste0(scratch, "/lead_pip_test_")))
unlink(resolved, recursive = TRUE)
cat("All lead selection, truth matching, tie, denominator, checkpoint, and plotting checks passed.\n")
