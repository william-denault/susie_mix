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
  list(susie_mix = list(pip = setNames(pip, predictor), converged = converged),
       weighted_fit_mix = list(pip = setNames(weighted, predictor), converged = converged),
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
  # The directory is enumerated, not assumed to contain 100 or 185 chunks.
  g1 <- list(Brain = make_tissue(c(.6, .3, .1), weighted = c(.01, .01, .98)),
             Liver = make_tissue(c(.1, .2, .7)))
  g2 <- list(Brain = make_tissue(c(.1, .1, .1, .1),
                                c("additive", "additive", "recessive", "dominant")),
             Liver = make_tissue(c(.2, .4, .4), converged = FALSE))
  attr(g2, "tissue_errors") <- list(Heart = "Too few usable predictors")
  saveRDS(g1, file.path(project, "results/G1.rds"))
  saveRDS(g2, file.path(project, "results/G2.rds"))
  # Missing gene outputs must block initial preparation, leaving no history.
  expect_error(em_prepare_iteration(project), "Missing: 1")
  history_file <- file.path(project, "results_em/prior_history.csv")
  stopifnot(!file.exists(history_file))
  saveRDS(list(gene = "G3", error = "No SNPs after QC"), file.path(project, "results/G3.rds"))
  first <- em_prepare_iteration(project)
  equal(first$chunks, 1:2)
  p1 <- read.csv(history_file)
  stopifnot(all(p1$iteration == 1), all(p1$source_fit == "susie_mix"))
  # Pool raw PIP masses, not normalized gene proportions or the weighted fits.
  equal(em_prior_for_tissue(p1, "Brain"), c(.8, .4, .2) / 1.4)
  equal(em_prior_for_tissue(p1, "Liver"), c(.3, .6, 1.1) / 2)
  equal(p1$n_source_gene_errors, c(1, 1))
  audit <- read.csv(file.path(first$iteration_dir, "source_audit.csv"))
  stopifnot(all(c("gene_error", "tissue_error", "nonconverged") %in% audit$issue))
  history_before <- readLines(history_file)
  expect_error(em_prepare_iteration(project), "unfinished")
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
  fake_run <- function(target_gene, tissue_priors, project_dir, temp_dir) {
    calls <<- c(calls, target_gene)
    equal(tissue_priors[c("tissue", "pi_add", "pi_rec", "pi_dom")],
          p1[c("tissue", "pi_add", "pi_rec", "pi_dom")])
    if (target_gene == "G3") stop("Still no SNPs after QC")
    list(Brain = make_tissue(c(.99, .005, .005), weighted = c(.2, .3, .5)),
         Liver = make_tissue(c(.99, .005, .005), weighted = c(.4, .4, .2)))
  }
  em_run_chunk(project, first$iteration_dir, 1L, fake_run)
  equal(calls, c("G1", "G2"))
  # Simulate a killed worker after saving its genes but before its marker.
  unlink(file.path(first$iteration_dir, "completed/chunk_001.done"))
  resume <- em_prepare_iteration(project, "resume")
  equal(resume$chunks, 1:2)
  stopifnot(identical(readLines(history_file), history_before))
  em_run_chunk(project, first$iteration_dir, 1L, fake_run)
  equal(calls, c("G1", "G2")) # Existing successful results were reused.
  equal(em_prepare_iteration(project, "resume")$chunks, 2L)
  expect_error(em_prepare_iteration(project), "unfinished")
  em_run_chunk(project, first$iteration_dir, 2L, fake_run)
  equal(calls, c("G1", "G2", "G3"))
  expect_error(em_prepare_iteration(project, "resume"), "complete")

  second <- em_prepare_iteration(project)
  history <- read.csv(history_file)
  p2 <- subset(history, iteration == 2)
  stopifnot(nrow(history) == 4L, all(p2$source_fit == "weighted_fit_mix"))
  equal(em_prior_for_tissue(p2, "Brain"), c(.2, .3, .5))
  equal(em_prior_for_tissue(p2, "Liver"), c(.4, .4, .2))
  equal(subset(history, iteration == 1), p1)
  stopifnot(dir.exists(file.path(second$iteration_dir, "results")))
  # Original scan and immutable first-iteration snapshot were preserved.
  stopifnot(identical(readRDS(file.path(project, "results/G1.rds")), g1))
  equal(read.csv(file.path(first$iteration_dir, "priors.csv")), p1)

  # Invalid/misaligned PIPs must never be silently omitted from pooled totals.
  fixture <- file.path(project, "bad.rds")
  bad <- g1
  bad$Brain$susie_mix$pip[1] <- NA_real_
  saveRDS(bad, fixture)
  expect_error(em_sum_pips(fixture, "susie_mix"), "invalid")
  bad <- g1
  bad$Brain$susie_mix$pip <- rev(bad$Brain$susie_mix$pip)
  saveRDS(bad, fixture)
  expect_error(em_sum_pips(fixture, "susie_mix"), "name mismatch")
  saveRDS(list(Brain = make_tissue(c(0, 0, 0))), fixture)
  expect_error(em_sum_pips(fixture, "susie_mix"), "undefined")
  zero <- em_sum_pips(fixture, "susie_mix", p1)$priors
  equal(em_prior_for_tissue(zero, "Brain"), em_prior_for_tissue(p1, "Brain"))
  stopifnot(all(zero$prior_status == "carried_forward_zero_pip"))
  saveRDS(list(error = "failed"), fixture)
  expect_error(em_sum_pips(fixture, "weighted_fit_mix", p1), "No usable")
  writeLines("not an RDS", fixture)
  expect_error(em_sum_pips(fixture, "susie_mix"), "unknown input format")

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
