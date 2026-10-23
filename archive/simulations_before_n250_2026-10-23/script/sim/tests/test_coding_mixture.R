# From the project root: Rscript --vanilla script/sim/tests/test_coding_mixture.R
options(susie_mix.coding_plots.run = FALSE)
source("script/sim/plot_coding_mixture.R")

equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))
expect_error <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
toy <- function(pip = c(.8, .1, 0, .1, .9, .1, .2, 0), coding_labels = FALSE) {
  x <- list(settings = list(L_add = 1, L_rec = 1, L_dom = 0, all_additive = FALSE,
                            n = 500, L = 10, pve = .2, min_maf = .05, hwe_thresh = 1e-8, min_n_rec = 5),
            seed = 1, true_pos = c(1L, 3L), true_pos_mix = c(1L, 5L),
            causal_coding = c("additive", "recessive"),
            susie_pip = c(.9, .1, .9), susie_mix_pip = pip,
            susie_mix_pip_snp = c(.95, .3, .95),
            mix_to_add = c(1L, 2L, 3L, 1L, 3L, 1L, 2L, 3L),
            metrics = data.frame(method = c("SuSiE", "SuSiE-mix"), converged = TRUE))
  if (coding_labels) x$mix_coding <- rep(cm_classes, c(3, 2, 3))
  x
}

x <- toy()
equal(cm_coding_map(x)$coding, rep(cm_classes, c(3, 2, 3)))
stopifnot(cm_coding_map(x)$origin == "reconstructed_blocks")
thresholds <- c(0, .1, .5, .8, .9, 1)
a <- cm_assess_replicate(x, .9, thresholds)
equal(a$metrics[c("mass_add", "mass_rec", "mass_dom")], c(.9, 1, .3))
equal(a$metrics[c("true_add", "true_rec", "true_dom")], c(1, 1, 0))
equal(a$metrics[c("conf_add_additive", "conf_rec_recessive")], c(1, 1))

# Exact threshold comparisons (including equality and zero) and decomposition
# into wrong coding at a causal SNP versus a wrong biological SNP.
labels <- cm_coding_map(x)$coding
truth <- seq_along(x$susie_mix_pip) %in% x$true_pos_mix
at_causal <- x$mix_to_add %in% x$true_pos
for (cl in seq_along(cm_classes)) for (t in seq_along(thresholds)) {
  selected <- labels == cm_classes[cl] & x$susie_mix_pip >= thresholds[t]
  direct <- c(sum(selected), sum(selected & truth), sum(selected & at_causal & !truth),
              sum(selected & !at_causal))
  equal(a$discovery[(cl - 1) * length(thresholds) + t, ], direct)
}
equal(sum(a$calibration[, "n"]), length(x$susie_mix_pip))
equal(sum(a$calibration[, "n_true"]), 2)
equal(sum(a$calibration[, "sum_pip"]), sum(x$susie_mix_pip))
equal(sum(a$calibration[, "sum_brier"]), sum((x$susie_mix_pip - truth)^2))
endpoints <- x; endpoints$susie_mix_pip[c(1, 2)] <- c(1, 0)
end <- cm_assess_replicate(endpoints)
equal(end$calibration[10, "n_true"], 1) # PIP=1 stays in the final additive bin.

# Ties never default to additive, and undetected causal SNPs stay in the
# denominator instead of inflating accuracy by being silently removed.
y <- x
y$susie_mix_pip[c(1, 6)] <- .5
y$susie_mix_pip_snp[3] <- .1
b <- cm_assess_replicate(y)
equal(b$metrics[c("conf_add_ambiguous", "conf_rec_undetected")], c(1, 1))
zero <- x; zero$susie_mix_pip[] <- 0; zero$susie_mix_pip_snp[] <- 0
z <- cm_assess_replicate(zero)
equal(z$metrics[c("conf_add_undetected", "conf_rec_undetected")], c(1, 1))

# Reject an R/D split that cannot be inferred from legacy compact data.
ambiguous <- x
ambiguous$true_pos <- 1L; ambiguous$true_pos_mix <- 1L; ambiguous$causal_coding <- "additive"
ambiguous$mix_to_add <- c(1L, 2L, 3L, 1L, 2L, 3L)
ambiguous$susie_mix_pip <- rep(.1, 6)
expect_error(cm_coding_map(ambiguous), "Ambiguous")
ambiguous$mix_coding <- rep(c("additive", "dominant"), each = 3)
equal(cm_coding_map(ambiguous)$coding, ambiguous$mix_coding)
bad <- x; bad$causal_coding[2] <- "dominant"
expect_error(cm_coding_map(bad), "disagree")
bad <- x; bad$susie_mix_pip[1] <- NA_real_
expect_error(cm_assess_replicate(bad), "Invalid")
roundoff <- x; roundoff$susie_mix_pip_snp[1] <- 1 + 2 * .Machine$double.eps
stopifnot(cm_assess_replicate(roundoff)$n_clamped == 1)
bad <- x; bad$susie_mix_pip_snp[1] <- 1.01
expect_error(cm_assess_replicate(bad), "Invalid")

# A small end-to-end checkpoint collection tests duplicate handling,
# all-additive controls, source errors, and every figure.
dir.create("tmp", showWarnings = FALSE)
scratch <- normalizePath("tmp", winslash = "/")
test_dir <- tempfile("coding_test_", tmpdir = scratch)
dir.create(file.path(test_dir, "chunks"), recursive = TRUE)
save_checkpoint <- function(results, prefix, requested) {
  file <- file.path(test_dir, "chunks", paste0(prefix, "_n500_L10_pve0.2_seed100_reps", requested, "_chunk1.RData"))
  save(results, file = file)
}
mixed <- lapply(1:4, function(seed) { z <- x; z$seed <- seed; z })
save_checkpoint(mixed, "additive_recessive_add1_rec1_dom0", 4)
save_checkpoint(mixed[1:2], "additive_recessive_add1_rec1_dom0", 2)
pure <- lapply(1:4, function(seed) {
  z <- toy(c(.8, .1, .9, .1, .2, .4, .2, .6), TRUE)
  z$seed <- seed; z$true_pos_mix <- c(1L, 3L); z$causal_coding <- rep("additive", 2)
  z$settings$L_add <- 2; z$settings$L_rec <- 0
  z
})
save_checkpoint(c(pure, list(list(seed = 5, error = "Known failed draw"))), "additive_add2_rec0_dom0", 5)
control <- mixed
for (i in seq_along(control)) {
  control[[i]]$settings$all_additive <- TRUE
  control[[i]]$causal_coding <- rep("additive", 2)
  control[[i]]$true_pos_mix <- c(1L, 3L)
}
save_checkpoint(control, "control_add1_rec1_dom0", 6)
config <- coding_plot_settings
config$chunk_dir <- file.path(test_dir, "chunks")
config$output_dir <- file.path(test_dir, "figures")
config$pve_values <- .2
config$bootstrap_reps <- 50L
config$write_png <- FALSE
analysis <- run_coding_plots(config)
stopifnot(nrow(analysis$replicates) == 12, sum(analysis$audit$duplicates) == 2,
          sum(analysis$audit$errors) == 1, sum(analysis$audit$invalid) == 0,
          length(list.files(config$output_dir, "[.]pdf$")) == 8)
equal(analysis$mixture$estimated_share[analysis$mixture$design == "add1_rec1_dom0" &
       !analysis$mixture$all_additive], rep(c(.9, 1, .3) / 2.2, 1))
stopifnot(all(analysis$mixture$true_share[analysis$mixture$all_additive &
                                      analysis$mixture$coding == "additive"] == 1))

# Pool raw PIPs across replicates, not replicate-normalized proportions.
r <- analysis$replicates[1:2, ]
r$scenario <- "Additive only"; r$configuration <- r$design <- "add1_rec0_dom0"; r$all_additive <- FALSE
r$mass_add <- c(.9, .01); r$mass_rec <- c(.1, .09); r$mass_dom <- 0
r$true_add <- 1; r$true_rec <- r$true_dom <- 0
summary <- cm_mixture_summary(r, 50)
equal(summary$estimated_share, c(.91, .19, 0) / 1.1)
twice <- cm_mixture_summary(rbind(r, r), 50)
equal(twice[c("estimated_share", "lower", "upper")], summary[c("estimated_share", "lower", "upper")])
no_calls <- analysis$discovery[analysis$discovery$threshold == 1, ]
stopifnot(all(is.na(no_calls$precision)), all(is.na(no_calls$fdp)))
no_truth <- analysis$discovery[analysis$discovery$scenario == "Additive only" &
                              analysis$discovery$coding != "additive", ]
stopifnot(all(is.na(no_truth$recall)))

resolved <- normalizePath(test_dir, winslash = "/", mustWork = TRUE)
stopifnot(startsWith(resolved, paste0(scratch, "/coding_test_")))
unlink(resolved, recursive = TRUE)
cat("All coding-map, classification, calibration, mixture, duplicate, and plotting checks passed.\n")
