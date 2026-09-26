# ROC definitions, pooling, partial checkpoints and paired method selection.
options(susie.init.roc.autorun = FALSE)
source("script/sim/plot_additive_init_roc.R")
expect_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(is.character(err), grepl(pattern, err, fixed = TRUE))
}
pip <- c(0, .2, .2, .5, 1)
truth <- c(2L, 5L)
thresholds <- c(0, .2, .5, 1)
counts <- additive_roc_counts(pip, truth, thresholds)
brute <- t(vapply(c(thresholds, Inf), function(t)
  c(tp = sum(pip[truth] >= t), fp = sum(pip[-truth] >= t)), c(tp = 0, fp = 0)))
stopifnot(identical(unname(counts), unname(brute)))
expect_error(additive_roc_counts(c(.1, NA), 1L), "Invalid saved SNP PIPs")
expect_error(additive_roc_counts(pip, c(1L, 1L)), "Invalid causal SNP indices")

test_project <- tempfile("init_roc_project_")
test_dir <- file.path(test_project, "simulation results/additive_slide_init_v3")
preview_dir <- file.path(Sys.getenv("TEMP", tempdir()), "codex-additive-roc-validation")
make_record <- function(seed, pips, truth, converged = c(TRUE, TRUE, TRUE)) {
  fits <- lapply(seq_along(pips), function(m) list(pip = pips[[m]], converged = converged[m]))
  names(fits) <- additive_roc_methods
  list(seed = seed, p = length(pips[[1L]]), true_pos = truth, fits = fits)
}
first <- make_record(1000001, list(
  c(.7, .8, .2, .6, .1, .01, .3, .55, .1, 0),
  c(.3, .85, .1, .12, .02, 0, .05, .8, .04, 0),
  c(.65, .9, .15, .11, .02, .1, .05, .6, .01, 0)), c(2L, 8L))
second <- make_record(1000002, list(
  c(.5, .8, .1, .3, .1, .01, 0),
  c(.8, .6, .01, .7, .1, 0, 0),
  c(.8, .55, .1, .5, .01, 0, 0)), c(1L, 4L), c(TRUE, FALSE, TRUE))
for (pve in c(.025, .05)) {
  path <- file.path(test_dir, sprintf("pve_%g/chunks", pve))
  dir.create(path, recursive = TRUE)
  settings <- list(schema = "additive_slide_init_v3", chunk = 1L, pve = pve,
    K = 2L, n = 500L, fit_L = 10L, reps_per_chunk = 4L,
    genotype_inputs = list(mode = "gtex_fallback", plink_memory = if (pve == .025) 2000L else 8000L))
  saveRDS(list(settings = settings, results = list(first, second,
    list(seed = 1000003, error = "PLINK test failure"), NULL)),
    file.path(path, "additive_init_chunk001.rds"))
}
out <- plot_additive_init_roc(test_dir, preview_dir)
stopifnot(all(file.exists(out$files)), nrow(out$replicates) == 8L,
  sum(out$replicates$status == "pending") == 2L,
  sum(out$replicates$status == "failed") == 2L,
  sum(out$replicates$status == "included_nonconverged") == 2L,
  all(out$curves$n_causal == 4L), all(out$curves$n_noncausal == 13L),
  all(out$curves$n_replicates == 2L))
# Pool actual causal and null counts, not the mean of per-replicate FPRs.
for (method in additive_roc_methods) {
  d <- out$curves[out$curves$pve == .025 & out$curves$method == method, ]
  causal <- c(first$fits[[method]]$pip[first$true_pos], second$fits[[method]]$pip[second$true_pos])
  null <- c(first$fits[[method]]$pip[-first$true_pos], second$fits[[method]]$pip[-second$true_pos])
  stopifnot(all(d$tp == vapply(d$threshold, function(t) sum(causal >= t), 0L)),
            all(d$fp == vapply(d$threshold, function(t) sum(null >= t), 0L)),
            all(d$tpr == d$tp / 4), all(d$fpr == d$fp / 13))
}
excluded <- plot_additive_init_roc(test_dir, file.path(preview_dir, "converged_only"),
                                 exclude_nonconverged = TRUE, write_png = FALSE)
stopifnot(all(excluded$curves$n_replicates == 1L), all(excluded$curves$n_noncausal == 8L))
# Exercise the exact source()-and-plot entry point used in RCC RStudio.
previous_root <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
Sys.setenv(SUSIE_MIX_PROJECT_DIR = test_project)
options(susie.init.roc.autorun = TRUE)
loaded <- new.env(parent = globalenv())
source("script/sim/plot_additive_init_roc.R", local = loaded)
stopifnot(all(file.exists(loaded$additive_init_roc$files)))
options(susie.init.roc.autorun = FALSE)
if (is.na(previous_root)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else
  Sys.setenv(SUSIE_MIX_PROJECT_DIR = previous_root)
single <- file.path(test_dir, "pve_0.025/chunks/additive_init_chunk001.rds")
duplicate <- file.path(dirname(single), "additive_init_chunk002.rds")
file.copy(single, duplicate)
expect_error(plot_additive_init_roc(test_dir, preview_dir, write_png = FALSE), "Duplicate PVE/seed")
unlink(duplicate)
saved <- readRDS(single)
saved$results <- list(list(seed = 1000001, error = "failed"), NULL, NULL, NULL)
saveRDS(saved, single)
partial <- plot_additive_init_roc(test_dir, file.path(preview_dir, "partial"), write_png = FALSE)
stopifnot(all(partial$curves$pve == .05))
expect_error(plot_additive_init_roc(test_dir, file.path(preview_dir, "no_success"),
             pves = .025, write_png = FALSE), "No successful replicates")
cat("PASS: tied PIPs, pooled SNP denominators, three paired methods, failed/pending records, convergence filtering, duplicate detection and plots.\n")
cat("PREVIEW: ", normalizePath(file.path(preview_dir, "roc_pip_zoom.png"), winslash = "/"), "\n", sep = "")
