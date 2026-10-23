# Saved-fit diagnostics: independently verify ranks, pooled counts, influence,
# CS denominators and selections; no external data or fitting packages needed.
source("script/sim/diagnose_additive_init.R")
near <- function(x, y) stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-12)))
expect_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(is.character(err), grepl(pattern, err, fixed = TRUE))
}
pip <- c(.2, .9, .2, 0, 1, .5)
truth <- c(1L, 2L)
dif <- outer(pip[truth], pip[-truth], `-`)
near(additive_diag_auc(pip, truth), mean((dif > 0) + .5 * (dif == 0)))
near(additive_diag_auc(rep(.5, 6L), truth), .5)
# Known curves (ascending thresholds, final Inf): diagonal, perfect, reversed.
near(additive_diag_pauc(c(2, 1, 0), c(2, 1, 0), 2, 2, .25), .125)
near(additive_diag_pauc(c(2, 2, 0), c(2, 0, 0), 2, 2, .25), 1)
near(additive_diag_pauc(c(2, 0, 0), c(2, 2, 0), 2, 2, .25), 0)

test_dir <- tempfile("init_diagnostics_")
preview_dir <- file.path(Sys.getenv("TEMP", tempdir()), "codex-additive-init-diagnostics-validation")
methods <- c("SuSiE", "SuSiE-slide", "SuSiE-init-slide")
make_record <- function(i) {
  p <- 8L + i %% 3L
  a <- c(.3, .2, .8, .7, rep(.1, p - 4L))
  b <- c(.31, .21, .8, .7, rep(.1, p - 4L))
  if (i == 1L) b <- c(.95, .9, rep(.01, p - 2L))
  if (i == 2L) b <- c(.01, .02, .8, .7, rep(.1, p - 4L))
  sets_a <- if (i == 1L) list(c(3L, 4L)) else list()
  sets_b <- if (i == 1L) list(1L, 2L) else list()
  sets <- function(x) list(cs = x, purity = data.frame(min.abs.corr = rep(.8, length(x))))
  fits <- list(list(pip = a, sets = sets(sets_a), converged = TRUE),
    list(pip = b, sets = sets(sets_b), converged = i != 12L),
    list(pip = a, sets = sets(sets_a), converged = TRUE))
  names(fits) <- methods
  list(seed = 1000000L + i, p = p, true_pos = truth, fits = fits,
    gene = paste0("GENE", i), causal_snps = c("chr1_100_A_G", "chr1_200_A_G"),
    causal_r = .1, beta_standardized = c(.1, -.1))
}
records <- lapply(1:12, make_record)
for (pve in c(.025, .05)) {
  path <- file.path(test_dir, sprintf("pve_%g/chunks", pve))
  dir.create(path, recursive = TRUE)
  settings <- list(schema = "additive_slide_init_v3", chunk = 1L, pve = pve,
    K = 2L, n = 500L, fit_L = 10L, reps_per_chunk = 14L,
    package_versions = c(susieR = "old.saved.version"),
    genotype_inputs = list(mode = "gtex_fallback", plink_memory = if (pve == .025) 2000L else 8000L))
  saveRDS(list(settings = settings, results = c(records, list(list(seed = 1000013, error = "test failure"), NULL))),
    file.path(path, "additive_init_chunk001.rds"))
}
files <- list.files(test_dir, recursive = TRUE, full.names = TRUE)
before <- tools::md5sum(files)
options(susie.init.roc.autorun = TRUE)
out <- diagnose_additive_initialization(test_dir, preview_dir, top_examples = 2L)
stopifnot(identical(before, tools::md5sum(files)), isTRUE(getOption("susie.init.roc.autorun")),
  nrow(out$per_dataset) == 48L, nrow(out$causal_snps) == 144L,
  all(out$overview$n_replicates == 12), all(out$cs_summary$n_replicates == 12),
  all(file.exists(file.path(preview_dir, c("diagnostics.pdf", "examples.pdf", "overview.csv")))))
audit <- read.csv(file.path(preview_dir, "roc_pip_replicates.csv"))
stopifnot(sum(audit$status == "included_nonconverged") == 2,
  sum(audit$status == "failed") == 2, sum(audit$status == "pending") == 2)
summary <- out$cs_summary[out$cs_summary$pve == .025, ]
near(summary$power, c(0, 2 / 24, 0))
near(summary$coverage, c(0, 1, 0))
near(summary$cs_fdr, c(1, 0, 1))
stopifnot(all(out$per_dataset$within_auc_gain[out$per_dataset$method == methods[3L]] == 0))
clear <- out$examples[out$examples$selection == "cs_win_no_extra_false_cs", ]
stopifnot(nrow(clear) == 2L, all(clear$seed == 1000001), all(clear$extra_recovered == 2))
loss <- out$examples[out$examples$selection == "within_dataset_ranking_loss", ]
stopifnot(nrow(loss) == 2L, all(loss$seed == 1000002))

# Brute-force SNP counts validate every removal subset, using the removed
# seeds from its audit. Different p across regions must retain pooled weights.
for (subset in unique(out$roc_sensitivity$subset)) for (method in methods) {
  s <- out$roc_sensitivity[out$roc_sensitivity$pve == .025 &
    out$roc_sensitivity$subset == subset & out$roc_sensitivity$method == method, ]
  remove <- if (nzchar(s$removed_seeds)) as.numeric(strsplit(s$removed_seeds, ";", fixed = TRUE)[[1L]]) else numeric()
  kept <- records[!vapply(records, function(r) r$seed %in% remove, TRUE)]
  pos <- unlist(lapply(kept, function(r) r$fits[[method]]$pip[truth]))
  neg <- unlist(lapply(kept, function(r) r$fits[[method]]$pip[-truth]))
  c <- out$roc_sensitivity_curves
  c <- c[c$pve == .025 & c$subset == subset & c$method == method, ]
  tp <- vapply(c$threshold, function(t) sum(pos >= t), 0L)
  fp <- vapply(c$threshold, function(t) sum(neg >= t), 0L)
  near(c$tp, tp); near(c$fp, fp)
  near(c$tpr, tp / length(pos)); near(c$fpr, fp / length(neg))
  near(s$mean_tpr_over_fpr_interval, additive_diag_pauc(tp, fp, length(pos), length(neg), .25))
}
influence <- out$roc_influence[out$roc_influence$pve == .025 & out$roc_influence$method == methods[2L], ]
stopifnot(influence$seed[which.max(influence$influence)] == 1000001)
removed <- out$roc_sensitivity[out$roc_sensitivity$pve == .025 &
  out$roc_sensitivity$subset == "remove_top_influence_1" & out$roc_sensitivity$method == methods[2L], ]
near(removed$gain_over_susie, influence$gap_without_dataset[influence$seed == 1000001])

excluded <- diagnose_additive_initialization(test_dir, file.path(preview_dir, "converged_only"),
  exclude_nonconverged = TRUE, top_examples = 0L)
stopifnot(all(excluded$overview$n_replicates == 11), nrow(excluded$examples) == 0L,
  identical(before, tools::md5sum(files)))
# One successful dataset, no CSs, no selected examples: undefined influence
# and coverage must stay NA rather than become spurious zero estimates.
single_dir <- tempfile("init_diagnostics_single_")
dir.create(single_dir)
single <- readRDS(files[1L]); single$results <- list(make_record(3L))
single$settings$reps_per_chunk <- 1L
saveRDS(single, file.path(single_dir, "additive_init_chunk001.rds"))
one <- diagnose_additive_initialization(single_dir, file.path(preview_dir, "single"), pves = .025)
stopifnot(all(is.na(one$cs_summary$coverage)), all(is.na(one$cs_summary$cs_fdr)),
  all(is.na(one$roc_influence$influence)), all(one$cs_summary$power == 0), nrow(one$examples) == 0L)
saved <- readRDS(files[1L]); saved$results[[1L]]$fits[[1L]]$sets <- NULL
saveRDS(saved, files[1L])
expect_error(diagnose_additive_initialization(test_dir, file.path(preview_dir, "missing_sets")),
  "Missing saved credible-set information")
cat("PASS: exact tied AUC, partial AUC, pooled removal curves, influence, CS power/coverage, example selection, convergence pairing, read-only checkpoints and plots.\n")
cat("PREVIEW: ", normalizePath(preview_dir, winslash = "/"), "\n", sep = "")
