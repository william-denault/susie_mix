# Objective accounting on known fits, separate from optimizer behavior.
for (f in c("em_utils.R", "slide_prior_em_utils.R")) source(file.path("script/scan_tissue_attempt", f))
fails <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
repo <- normalizePath(".", winslash = "/")
dir.create(file.path(repo, "tmp"), showWarnings = FALSE)
root <- tempfile("slide_objective_", tmpdir = file.path(repo, "tmp")); dir.create(root)
model <- function(prior, elbo, V = 1) {
  post <- prior * c(2, rep(1, 16)); post <- post / sum(post)
  joint <- array(0, c(1L, 2L, 17L)); joint[1, 1, ] <- .7 * post; joint[1, 2, 9] <- .3
  structure(list(alpha = matrix(c(.7, .3), 1L), alpha_delta = joint,
    delta_prior_counts = matrix(.7 * post, 1L), V = V, delta_forced = c(FALSE, TRUE),
    input_p = 2L, delta_grid = spe_grid(), delta_prior = prior, converged = TRUE,
    elbo = c(elbo - 10, elbo)), class = c("susie_slide", "susie"))
}
setup <- function(iteration, priors, values) {
  d <- file.path(root, "results_slide_prior_em", sprintf("iteration_%03d", iteration))
  dir.create(file.path(d, "results"), recursive = TRUE)
  priors$iteration <- iteration
  write.csv(priors, file.path(d, "priors.csv"), row.names = FALSE)
  saveRDS(list(Brain = list(fit_slide_prior = model(spe_prior(priors, "Brain"), values[1])),
               Liver = list(fit_slide_prior = model(spe_prior(priors, "Liver"), values[2]))),
          file.path(d, "results/A.rds"))
  # A zero-active-component fit still contributes its full ELBO.
  saveRDS(list(Brain = list(fit_slide_prior = model(spe_prior(priors, "Brain"), values[3], V = 0))),
          file.path(d, "results/B.rds"))
  saveRDS(list(error = "Fixture frozen exclusion"), file.path(d, "results/C.rds"))
  files <- list.files(file.path(d, "results"), full.names = TRUE)
  list(dir = d, files = files, pooled = spe_pool(files, priors))
}
w <- matrix(1/17, 2L, 17L, dimnames = list(NULL, spe_columns()))
zero <- setup(0L, data.frame(tissue = c("Brain", "Liver"), w), c(-100, -50, -25))
s0 <- spe_record_objective(root, zero$dir, zero$pooled)
stopifnot(s0$elbo[1] == -175, s0$n_fits[1] == 3L, s0$n_genes[1] == 2L,
          s0$elbo[s0$tissue == "Brain"] == -125, is.na(s0$delta_elbo[1]),
          is.na(s0$prior_max_change[1]), sum(zero$pooled$priors$n_active) == 2L)
history_file <- file.path(root, "results_slide_prior_em/objective_history.csv")
spe_record_objective(root, zero$dir, zero$pooled)
stopifnot(nrow(read.csv(history_file)) == 3L)
one <- setup(1L, zero$pooled$priors, c(-95, -48, -25))
s1 <- spe_record_objective(root, one$dir, one$pooled)
stopifnot(s1$elbo[1] == -168, s1$delta_elbo[1] == 7,
          abs(s1$relative_delta_elbo[1] - 7/175) < 1e-12,
          abs(s1$delta_elbo_per_fit[1] - 7/3) < 1e-12,
          s1$prior_max_change[1] > 0, !s1$elbo_decreased[1],
          s1$next_m_step_q_gain[1] == sum(one$pooled$priors$m_step_q_gain))
hash <- tools::md5sum(history_file)
# Equal fit counts are insufficient: the exact gene/tissue cohort must match.
bad <- one$pooled; bad$audit$gene[1] <- "DIFFERENT"
fails(spe_record_objective(root, one$dir, bad), "cohort changed")
bad <- one$pooled; bad$audit <- bad$audit[-1, ]
fails(spe_record_objective(root, one$dir, bad), "cohort changed")
stopifnot(identical(hash, tools::md5sum(history_file)))
# A changed earlier total cannot silently invalidate already-recorded deltas.
bad <- zero$pooled; bad$audit$elbo[1] <- -1000
fails(spe_record_objective(root, zero$dir, bad), "earlier objective")
spe_record_objective(root, zero$dir, zero$pooled)
two <- setup(2L, one$pooled$priors, c(-107, -48, -25))
warned <- FALSE
s2 <- withCallingHandlers(spe_record_objective(root, two$dir, two$pooled), warning = function(w) {
  warned <<- grepl("ELBO decreased", conditionMessage(w)); invokeRestart("muffleWarning")
})
stopifnot(warned, s2$delta_elbo[1] == -12, s2$elbo_decreased[1],
          nrow(read.csv(history_file)) == 9L,
          file.exists(file.path(two$dir, "objective_summary.csv")),
          !dir.exists(file.path(root, "results_slide_prior_em/iteration_003")))
bad <- readRDS(two$files[1]); bad$Brain$fit_slide_prior$elbo <- NA_real_
saveRDS(bad, two$files[1])
fails(spe_pool(two$files, one$pooled$priors), "Missing finite ELBO")
# The real CLI refuses to prepare beyond a target even on a repeated "new".
scripts <- file.path(root, "script/scan_tissue_attempt"); dir.create(scripts, recursive = TRUE)
invisible(file.copy(file.path(repo, "script/scan_tissue_attempt",
  c("em_utils.R", "slide_prior_em_utils.R", "slide_prior_recovery.R")), scripts))
write.csv(data.frame(iteration = 20L), file.path(root, "results_slide_prior_em/prior_history.csv"), row.names = FALSE)
old_target <- Sys.getenv("SUSIE_SLIDE_TARGET_ITERATION", unset = NA_character_)
Sys.setenv(SUSIE_SLIDE_TARGET_ITERATION = "20")
cli <- file.path(repo, "script/scan_tissue_attempt/recover_slide_prior_em.R")
run_r <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
answer <- suppressWarnings(system2(run_r, c("--vanilla", shQuote(cli), "new", shQuote(root), "20"),
                                  stdout = TRUE, stderr = TRUE))
if (is.na(old_target)) Sys.unsetenv("SUSIE_SLIDE_TARGET_ITERATION") else
  Sys.setenv(SUSIE_SLIDE_TARGET_ITERATION = old_target)
stopifnot(attr(answer, "status") == 1L, any(grepl("Target iteration already prepared", answer)),
          !dir.exists(file.path(root, "results_slide_prior_em/iteration_021")))
stopifnot(identical(dirname(normalizePath(root, winslash = "/")), paste0(repo, "/tmp")))
unlink(root, recursive = TRUE)
cat("PASS: full ELBO sums, zero-active fits, cohort identity, correct iteration/deltas, prior movement, idempotence, decreases and final-iteration recording.\n")
