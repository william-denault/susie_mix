# Run from the project root. Base R; never reads or writes real scan results.
for (file in c("em_utils.R", "slide_prior_em_utils.R", "prepare_slide_prior_em_iteration.R",
               "run_slide_prior_em_chunk.R")) source(file.path("script/scan_tissue_attempt", file))
fails <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
uniform <- rep(1 / 17, 17)
make_fit <- function(prior = uniform, uninformative = FALSE) {
  a <- cbind(free = c(.2, .8, .5), forced = c(.5, .1, .3), null = c(.3, .1, .2))
  joint <- array(0, c(3L, 3L, 17L))
  q <- matrix(0, 3, 17); q[1, 1] <- 1; q[2, 17] <- 1; q[3, 9] <- 1
  if (uninformative) q <- matrix(prior, 3, 17, byrow = TRUE)
  for (l in 1:3) {
    joint[l, 1, ] <- a[l, 1] * q[l, ]
    joint[l, 2:3, 9] <- a[l, 2:3]
  }
  counts <- a[, 1] * q; counts[3, ] <- 0
  structure(list(alpha = a, alpha_delta = joint, delta_prior_counts = counts,
                 delta_forced = c(FALSE, TRUE, TRUE), input_p = 2L,
                 delta_grid = spe_grid(), delta_prior = prior, V = c(1, 1, 0),
                 converged = TRUE, elbo = -100, sets = NULL), class = c("susie_slide", "susie"))
}
fit <- make_fit()
counts <- spe_fit_counts(fit, uniform)
stopifnot(all.equal(colSums(counts)[c(1, 17)], c(.2, .8)), sum(counts) == 1,
          all(rowSums(counts) == c(.2, .8, 0)))
expected <- numeric(17); expected[c(1, 17)] <- c(.2, .8)
stopifnot(isTRUE(all.equal(spe_update_weights(colSums(counts), uniform), expected)))
uninformative <- spe_fit_counts(make_fit(uninformative = TRUE), uniform)
stopifnot(isTRUE(all.equal(spe_update_weights(colSums(uninformative), uniform), uniform)),
          identical(spe_update_weights(numeric(17), uniform), uniform))
bad <- fit; bad$delta_grid <- seq(-1, 1, length.out = 9)
fails(spe_fit_counts(bad, uniform), "frozen grid")
bad <- fit; bad$converged <- FALSE
fails(spe_fit_counts(bad, uniform), "Unconverged")
bad <- fit; bad$delta_prior_counts <- NULL
fails(spe_fit_counts(bad, uniform), "delta_prior_counts")
bad <- fit; bad$delta_prior_counts[1, ] <- bad$delta_prior_counts[1, ] / .2
fails(spe_fit_counts(bad, uniform), "free-SNP")
old_budget <- Sys.getenv("SUSIE_SLIDE_MAX_ITER", unset = NA_character_)
Sys.setenv(SUSIE_SLIDE_MAX_ITER = "5000")
stopifnot(spe_fit_budget(1000L) == 5000L)
Sys.setenv(SUSIE_SLIDE_MAX_ITER = "500")
fails(spe_fit_budget(1000L), "at least")
if (is.na(old_budget)) Sys.unsetenv("SUSIE_SLIDE_MAX_ITER") else
  Sys.setenv(SUSIE_SLIDE_MAX_ITER = old_budget)

root <- tempfile("slider_em_test_")
dir.create(file.path(root, "results_em/iteration_012"), recursive = TRUE)
mix_history <- data.frame(iteration = 12L, tissue = c("Brain", "Liver"))
write.csv(mix_history, file.path(root, "results_em/prior_history.csv"), row.names = FALSE)
write.csv(data.frame(chunk = 1L, gene = c("A", "B", "C")),
          file.path(root, "results_em/iteration_012/manifest.csv"), row.names = FALSE)
mix_hash <- tools::md5sum(list.files(file.path(root, "results_em"), recursive = TRUE, full.names = TRUE))
prep <- function(mode = "new") spe_prepare_iteration(root, mode, n_chunks = 1L, package_version = "fixture")
init <- prep()
stopifnot(basename(init$iteration_dir) == "iteration_000")
prior0 <- read.csv(file.path(init$iteration_dir, "priors.csv"))
stopifnot(max(abs(as.matrix(prior0[spe_columns()]) - 1 / 17)) < 1e-14)
frozen <- tools::md5sum(file.path(init$iteration_dir, c("priors.csv", "manifest.csv", "settings.rds")))
fails(prep(), "unfinished")
resume <- prep("resume")
stopifnot(identical(init, resume), identical(frozen, tools::md5sum(names(frozen))))
calls <- character()
nonconverged <- TRUE
fake <- function(target_gene, tissue_priors, previous_result, ...) {
  calls <<- c(calls, target_gene)
  if (target_gene == "C") stop("Fixture baseline gene failure")
  f <- make_fit(spe_prior(tissue_priors, "Brain"))
  if (target_gene == "B" && nonconverged) f$converged <- FALSE
  list(Brain = list(fit_slide_prior = f, n_ind = 51, min_pv = .9, median_read = 1))
}
fails(spe_run_chunk(root, init$iteration_dir, 1L, fake), "Unconverged")
stopifnot(length(em_pending_chunks(init$iteration_dir, data.frame(chunk = 1L))) == 1)
nonconverged <- FALSE
spe_run_chunk(root, init$iteration_dir, 1L, fake)
stopifnot(sum(calls == "A") == 1L, sum(calls == "B") == 2L)
init_hash <- tools::md5sum(list.files(init$iteration_dir, recursive = TRUE, full.names = TRUE))
one <- prep()
prior1 <- read.csv(file.path(one$iteration_dir, "priors.csv"))
stopifnot(isTRUE(all.equal(spe_prior(prior1, "Brain"), expected)),
          isTRUE(all.equal(spe_prior(prior1, "Liver"), uniform)),
          prior1$n_fits[1] == 2L, prior1$n_active[1] == 4L,
          prior1$free_mass[1] == 2, prior1$m_step_q_gain[1] > 0)
# Losing a previously fitted tissue cannot silently change the EM cohort.
fails(spe_run_chunk(root, one$iteration_dir, 1L, function(...) list()), "cohort changed")
spe_run_chunk(root, one$iteration_dir, 1L, fake)
two <- prep()
stopifnot(basename(two$iteration_dir) == "iteration_002",
          identical(init_hash, tools::md5sum(names(init_hash))),
          identical(mix_hash, tools::md5sum(names(mix_hash))))
# Incomplete output cannot be used to learn another prior, even with markers.
em_atomic_write(list(), file.path(two$iteration_dir, "completed/chunk_001.done"))
fails(prep(), "Missing:")
balanced <- em_rechunk_manifest(data.frame(chunk = 1L, gene = paste0("G", 1:18468)), 298L)
stopifnot(length(unique(balanced$chunk)) == 298L, min(table(balanced$chunk)) == 61L,
          max(table(balanced$chunk)) == 62L)
cat("PASS: counts, forced/null/zero-V exclusions, uncertain fits, frozen cohort, initialization, resume and immutable snapshots.\n")
cat("Temporary test fixture:", root, "\n")
