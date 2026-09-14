# Run from the repository root:
# Rscript --vanilla script/scan_tissue_attempt/tests/test_em_iteration.R
# Synthetic two-iteration integration test; no GTEx, Slurm, or R packages needed.
source("script/scan_tissue_attempt/em_utils.R")
source("script/scan_tissue_attempt/prepare_em_iteration.R")
source("script/scan_tissue_attempt/run_em_chunk.R")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected error matching: ", pattern)
  }
}
equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))

make_tissue <- function(pip, coding = c("additive", "recessive", "dominant"),
                        weighted = pip, converged = TRUE) {
  predictor <- paste0("chr1_", seq_along(pip), "_A_C_b38_A__", coding)
  make_fit <- function(p) {
    a <- matrix(if (sum(p) > 0) p / sum(p) else rep(1/length(p), length(p)), 1,
                dimnames = list(NULL, predictor))
    structure(list(pip = setNames(p, predictor), converged = converged, alpha = a,
                   V = as.numeric(sum(p) > 0), mu = a * 0, mu2 = a * 0 + .1,
                   sigma2 = 1, pi = rep(1/length(p), length(p))), class = "susie")
  }
  list(susie_mix = make_fit(pip), weighted_fit_mix = make_fit(weighted),
       mix_coding = coding,
       mix_predictor_map = data.frame(predictor_index = seq_along(pip),
                                      predictor_name = predictor, coding = coding))
}

run_tests <- function() {
  dir.create("tmp", showWarnings = FALSE)
  test_base <- normalizePath("tmp", winslash = "/", mustWork = TRUE)
  project <- tempfile("em_test_", tmpdir = test_base)
  dir.create(file.path(project, "data/temp_index"), recursive = TRUE)
  dir.create(file.path(project, "results"))
  on.exit({
    resolved <- normalizePath(project, winslash = "/", mustWork = TRUE)
    stopifnot(startsWith(resolved, paste0(test_base, "/em_test_")))
    unlink(resolved, recursive = TRUE)
  }, add = TRUE)
  writeLines(c("G1", "G2"), file.path(project, "data/temp_index/chunk_001_genes.txt"))
  writeLines("G3", file.path(project, "data/temp_index/chunk_185_genes.txt"))
  # Supply synthetic annotation; production reads the same GTF as the worker.
  prepare <- function(mode = "new") em_prepare_iteration(project, mode,
    gene_annotations = data.frame(gene_name = c("G1", "G2", "G3"), chromosome = c("chr1", "2", "22")))
  # The directory is enumerated, not assumed to contain 100 or 185 chunks.
  g1 <- list(Brain = make_tissue(c(.6, .3, .1), weighted = c(.01, .01, .98)),
             Liver = make_tissue(c(.1, .2, .7)))
  g2 <- list(Brain = make_tissue(c(.1, .1, .1, .1),
                                c("additive", "additive", "recessive", "dominant")),
             Liver = make_tissue(c(.2, .4, .4)))
  attr(g2, "tissue_errors") <- list(Heart = "Too few usable predictors")
  saveRDS(g1, file.path(project, "results/G1.rds"))
  saveRDS(g2, file.path(project, "results/G2.rds"))
  # Missing gene outputs must block initial preparation, leaving no history.
  expect_error(prepare(), "Missing: 1")
  history_file <- file.path(project, "results_em/prior_history.csv")
  stopifnot(!file.exists(history_file))
  saveRDS(list(gene = "G3", error = "No SNPs after QC"), file.path(project, "results/G3.rds"))
  first <- prepare()
  equal(first$chunks, 1:2)
  p1 <- read.csv(history_file)
  stopifnot(all(p1$iteration == 1), all(p1$source_fit == "susie_mix"))
  # Use component assignment probabilities, with unequal coding block sizes.
  equal(em_prior_for_tissue(p1, "Brain"), c(1.1, .55, .35) / 2)
  equal(em_prior_for_tissue(p1, "Liver"), c(.3, .6, 1.1) / 2)
  equal(p1$n_source_gene_errors, c(1, 1))
  audit <- read.csv(file.path(first$iteration_dir, "source_audit.csv"))
  stopifnot(all(c("gene_error", "tissue_error") %in% audit$issue))
  stopifnot(all(p1$update_method == "susie_active_component_alpha_v2"),
            all(p1$n_active_components == p1$n_components),
            all(p1$n_active_fits == p1$n_fits))
  # Simulate the schema of a pre-correction iteration. New preparation must
  # append tagged rows without rewriting its frozen snapshot or old values.
  legacy_columns <- c("iteration", "source_iteration", "source_fit", "tissue", "pi_add", "pi_rec", "pi_dom",
                      "pip_add", "pip_rec", "pip_dom", "pip_total", "n_fits", "n_nonconverged", "prior_status",
                      "n_source_files", "n_source_gene_errors", "n_source_tissue_errors", "created_at")
  p1 <- p1[legacy_columns]
  write.csv(p1, history_file, row.names = FALSE)
  write.csv(p1, file.path(first$iteration_dir, "priors.csv"), row.names = FALSE)
  history_before <- readLines(history_file)
  expect_error(prepare(), "unfinished")
  stopifnot(identical(readLines(history_file), history_before))

  # Class prior mass survives unequal predictor counts, and absent classes
  # redistribute their mass only among the available classes.
  weights <- em_predictor_weights(c("additive", "recessive", "additive", "dominant"),
                                   c(dominant = .1, recessive = .3, additive = .6))
  equal(weights, c(.3, .3, .3, .1))
  equal(em_predictor_weights(c("dominant", "additive", "additive"),
                            c(additive = .6, recessive = .3, dominant = .1)), c(1, 3, 3) / 7)
  expect_error(em_predictor_weights("dominant", c(additive = 1, recessive = 0, dominant = 0)),
               "zero total prior")
  expect_error(em_prior_for_tissue(p1, "Missing"), "No estimated")
  expect_error(em_validate_priors(rbind(p1, p1[1, ])), "unique tissue")

  # Exercise the actual chunk runner with a deterministic stand-in for GTEx.
  calls <- character()
  fake_run <- function(target_gene, tissue_priors, project_dir, temp_dir, previous_result) {
    calls <<- c(calls, target_gene)
    equal(tissue_priors[c("tissue", "pi_add", "pi_rec", "pi_dom")],
          p1[c("tissue", "pi_add", "pi_rec", "pi_dom")])
    if (target_gene == "G3") stop("Still no SNPs after QC")
    stopifnot(!is.null(previous_result$Brain$susie_mix$alpha))
    list(Brain = make_tissue(c(.99, .005, .005), weighted = c(.2, .3, .5)),
         Liver = make_tissue(c(.99, .005, .005), weighted = c(.4, .4, .2)))
  }
  em_run_chunk(project, first$iteration_dir, 1L, fake_run)
  equal(calls, c("G1", "G2"))
  timing <- readRDS(file.path(first$iteration_dir, "completed/chunk_001.done"))
  stopifnot(grepl("UTC$", timing$started_at), grepl("UTC$", timing$finished_at),
            is.finite(timing$elapsed_seconds), timing$elapsed_seconds >= 0)
  # Simulate a killed worker after saving its genes but before its marker.
  unlink(file.path(first$iteration_dir, "completed/chunk_001.done"))
  resume <- prepare("resume")
  equal(resume$chunks, 1:2)
  stopifnot(identical(readLines(history_file), history_before))
  em_run_chunk(project, first$iteration_dir, 1L, fake_run)
  equal(calls, c("G1", "G2")) # Existing successful results were reused.
  equal(prepare("resume")$chunks, 2L)
  expect_error(prepare(), "unfinished")
  em_run_chunk(project, first$iteration_dir, 2L, fake_run)
  equal(calls, c("G1", "G2", "G3"))
  expect_error(prepare("resume"), "complete")

  second <- prepare()
  history <- read.csv(history_file)
  p2 <- subset(history, iteration == 2)
  stopifnot(nrow(history) == 4L, all(p2$source_fit == "weighted_fit_mix"))
  equal(em_prior_for_tissue(p2, "Brain"), c(.2, .3, .5))
  equal(em_prior_for_tissue(p2, "Liver"), c(.4, .4, .2))
  equal(subset(history, iteration == 1)[names(p1)], p1)
  stopifnot(all(subset(history, iteration == 1)$update_method == "legacy_pip_share"),
            all(p2$update_method == "susie_active_component_alpha_v2"))
  stopifnot(dir.exists(file.path(second$iteration_dir, "results")))
  # Original scan and immutable first-iteration snapshot were preserved.
  stopifnot(identical(readRDS(file.path(project, "results/G1.rds")), g1))
  equal(read.csv(file.path(first$iteration_dir, "priors.csv")), p1)

  # Invalid/misaligned PIPs must never be silently omitted from pooled totals.
  fixture <- file.path(project, "bad.rds")
  bad <- g1
  bad$Brain$susie_mix$pip[1] <- NA_real_
  saveRDS(bad, fixture)
  expect_error(em_estimate_coding_priors(fixture, "susie_mix"), "invalid")
  bad <- g1
  bad$Brain$susie_mix$pip <- rev(bad$Brain$susie_mix$pip)
  saveRDS(bad, fixture)
  expect_error(em_estimate_coding_priors(fixture, "susie_mix"), "name mismatch")
  saveRDS(list(Brain = make_tissue(c(0, 0, 0))), fixture)
  zero <- em_estimate_coding_priors(fixture, "susie_mix", p1)$priors
  equal(em_prior_for_tissue(zero, "Brain"), em_prior_for_tissue(p1, "Brain"))
  equal(em_prior_for_tissue(zero, "Liver"), em_prior_for_tissue(p1, "Liver"))
  stopifnot(zero$n_zero_variance_components[zero$tissue == "Brain"] == 1L,
            all(zero$n_active_components == 0L), all(zero$n_active_fits == 0L),
            all(zero$alpha_total == 0), all(zero$mstep_q_gain == 0),
            all(zero$prior_status == "carried_forward_no_active_components"))
  initial_null <- em_estimate_coding_priors(fixture, "susie_mix")$priors
  equal(em_prior_for_tissue(initial_null, "Brain"), rep(1/3, 3))
  stopifnot(initial_null$prior_status == "initialized_uniform_no_active_components")

  # A signal fit plus an all-null fit uses only the signal's active rows.
  # CS labels deliberately disagree with activity: V, not CS/PIP, selects rows.
  signal <- make_tissue(c(.7, .2, .1))
  signal$susie_mix$alpha <- rbind(signal$susie_mix$alpha, c(.1, .3, .6), c(.99, .005, .005))
  signal$susie_mix$V <- c(1, 1e-12, 0) # Tiny positive V is still included.
  signal$susie_mix$sets <- list(cs = list(L3 = 1L), cs_index = 3L)
  null <- make_tissue(c(0, 0, 0))
  null$susie_mix$alpha <- null$susie_mix$alpha[rep(1L, 10), , drop = FALSE]
  null$susie_mix$V <- 0 # Scalar zero broadcasts to every row.
  null$susie_mix$sets <- list(cs = NULL)
  signal_file <- file.path(project, "signal.rds")
  null_file <- file.path(project, "null.rds")
  saveRDS(list(Brain = signal), signal_file)
  saveRDS(list(Brain = null), null_file)
  active_only <- em_estimate_coding_priors(c(signal_file, null_file), "susie_mix", p1)$priors
  brain <- active_only[active_only$tissue == "Brain", ]
  equal(em_prior_for_tissue(active_only, "Brain"), c(.4, .25, .35))
  stopifnot(brain$n_fits == 2L, brain$n_active_fits == 1L,
            brain$n_components == 13L, brain$n_active_components == 2L,
            brain$n_zero_variance_components == 11L, brain$alpha_total == 2)
  # All positive rows without a CS remain eligible; scalar V also broadcasts.
  signal$susie_mix$V <- 1
  signal$susie_mix$sets <- list(cs = NULL)
  saveRDS(list(Brain = signal), signal_file)
  all_positive <- em_estimate_coding_priors(signal_file, "susie_mix")$priors
  equal(em_prior_for_tissue(all_positive, "Brain"), colSums(signal$susie_mix$alpha) / 3)
  for (invalid_v in list(c(1, NA, 0), c(1, -1, 0), c(1, Inf, 0), c(1, 0))) {
    signal$susie_mix$V <- invalid_v
    saveRDS(list(Brain = signal), signal_file)
    expect_error(em_estimate_coding_priors(signal_file, "susie_mix"), "prior variances")
  }

  # Null rows must not connect otherwise missing coding classes in the M-step.
  ar <- make_tissue(c(.8, .2), c("additive", "recessive"))
  rd <- make_tissue(c(.2, .8), c("recessive", "dominant"))
  saveRDS(list(Brain = ar), signal_file)
  saveRDS(list(Brain = rd), fixture)
  incomplete <- em_estimate_coding_priors(c(signal_file, fixture, null_file), "susie_mix")
  stopifnot(max(abs(em_prior_for_tissue(incomplete$priors, "Brain") - c(4, 1, 4) / 9)) < 1e-6)
  stopifnot(sum(incomplete$component_counts$Brain[7, ]) == 0,
            incomplete$priors$n_active_components == 2L)
  bad <- g1
  bad$Brain$susie_mix$alpha <- NULL
  saveRDS(bad, fixture)
  expect_error(em_estimate_coding_priors(fixture, "susie_mix"), "PIPs alone")
  bad <- g1
  bad$Brain$susie_mix$converged <- FALSE
  saveRDS(bad, fixture)
  expect_error(em_estimate_coding_priors(fixture, "susie_mix"), "not confirmed converged")
  bad <- g1
  colnames(bad$Brain$susie_mix$alpha) <- rev(colnames(bad$Brain$susie_mix$alpha))
  saveRDS(bad, fixture)
  expect_error(em_estimate_coding_priors(fixture, "susie_mix"), "Alpha/map name mismatch")
  with_elbo <- g1
  with_elbo$Brain$susie_mix$elbo <- c(-12, -10)
  saveRDS(with_elbo, fixture)
  diagnostics <- em_estimate_coding_priors(fixture, "susie_mix")$priors
  stopifnot(diagnostics$source_elbo_sum[diagnostics$tissue == "Brain"] == -10,
            diagnostics$n_source_elbo[diagnostics$tissue == "Brain"] == 1L,
            is.na(diagnostics$source_elbo_sum[diagnostics$tissue == "Liver"]))
  saveRDS(list(error = "failed"), fixture)
  expect_error(em_estimate_coding_priors(fixture, "weighted_fit_mix", p1), "No usable")
  writeLines("not an RDS", fixture)
  expect_error(em_estimate_coding_priors(fixture, "susie_mix"), "unknown input format")

  # An unconverged E-step cannot create a done marker or feed an automatic
  # next iteration. Resume retries the affected gene instead of reusing it.
  fail_fit <- TRUE
  nonconverged_run <- function(target_gene, ...) {
    list(Brain = make_tissue(c(.2, .3, .5), converged = !fail_fit))
  }
  expect_error(em_run_chunk(project, second$iteration_dir, 1L, nonconverged_run), "unconverged")
  stopifnot(!file.exists(file.path(second$iteration_dir, "completed/chunk_001.done")))
  equal(prepare("resume")$chunks, 1:2)
  fail_fit <- FALSE
  em_run_chunk(project, second$iteration_dir, 1L, nonconverged_run)
  equal(prepare("resume")$chunks, 2L)
  invalid_model <- function(...) em_stop_fit("Invalid warm-start predictor order")
  expect_error(em_run_chunk(project, second$iteration_dir, 2L, invalid_model), "Invalid warm-start")
  stopifnot(!file.exists(file.path(second$iteration_dir, "completed/chunk_002.done")))
  em_run_chunk(project, second$iteration_dir, 2L, nonconverged_run)

  # Simulate a completed all-row-alpha v1 iteration. Additional diagnostic
  # columns may be present as NA in history, but not in its frozen snapshot.
  p2_v1 <- p2[setdiff(names(p2), c("n_active_components", "n_active_fits"))]
  p2_v1$update_method <- "susie_component_alpha_v1"
  write.csv(p2_v1, file.path(second$iteration_dir, "priors.csv"), row.names = FALSE)
  history_v1 <- em_bind_history(read.csv(history_file)[1:2, ], p2_v1)
  bad_history <- history_v1
  bad_history$n_active_components[bad_history$iteration == 2] <- 1
  write.csv(bad_history, history_file, row.names = FALSE)
  expect_error(prepare(), "differ from")
  write.csv(history_v1, history_file, row.names = FALSE)
  third <- prepare()
  p3 <- read.csv(file.path(third$iteration_dir, "priors.csv"))
  stopifnot(all(p3$source_iteration == 2L), all(p3$update_method == "susie_active_component_alpha_v2"))
  equal(read.csv(file.path(second$iteration_dir, "priors.csv")), p2_v1)
  migrated <- subset(read.csv(history_file), iteration == 2)
  equal(migrated[names(p2_v1)], p2_v1)
  stopifnot(all(is.na(migrated$n_active_components)), all(is.na(migrated$n_active_fits)))
  equal(read.csv(file.path(first$iteration_dir, "priors.csv")), p1)

  # Confirm the EM workhorse parses and contains only one (weighted) SuSiE call.
  code <- parse("script/scan_tissue_attempt/workhorse_em.R")
  susie_calls <- list()
  visit <- function(expr) {
    if (missing(expr)) return(invisible(NULL))
    if (is.call(expr) && identical(expr[[1]], as.name("susie"))) {
      susie_calls[[length(susie_calls) + 1L]] <<- expr
    }
    if (is.call(expr) || is.expression(expr)) for (item in as.list(expr)) visit(item)
  }
  visit(code)
  stopifnot(length(susie_calls) == 1L, "prior_weights" %in% names(susie_calls[[1]]))
  cat("All EM iteration tests passed.\n")
}
run_tests()
