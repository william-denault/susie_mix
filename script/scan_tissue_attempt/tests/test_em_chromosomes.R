# Rscript --vanilla script/scan_tissue_attempt/tests/test_em_chromosomes.R
# Verify the scientific exclusion through preparation and chunk execution.
source("script/scan_tissue_attempt/em_utils.R")
source("script/scan_tissue_attempt/prepare_em_iteration.R")
source("script/scan_tissue_attempt/run_em_chunk.R")

run_tests <- function() {
  equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))
  expect_error <- function(expr, pattern) {
    e <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
  }
  dir.create("tmp", showWarnings = FALSE)
  base <- normalizePath("tmp", winslash = "/")
  project <- tempfile("em_chromosomes_", tmpdir = base)
  dir.create(file.path(project, "data/temp_index"), recursive = TRUE)
  dir.create(file.path(project, "results"))
  on.exit({
    resolved <- normalizePath(project, winslash = "/", mustWork = TRUE)
    stopifnot(startsWith(resolved, paste0(base, "/em_chromosomes_")))
    unlink(resolved, recursive = TRUE)
  }, add = TRUE)
  annotations <- data.frame(gene_name = c("AUTO", "GX", "GY", "GMT"),
                            chromosome = c("chr22", "chrX", "24", "chrM"))
  writeLines(c("AUTO", "GX"), file.path(project, "data/temp_index/chunk_001_genes.txt"))
  writeLines(c("GY", "GMT"), file.path(project, "data/temp_index/chunk_002_genes.txt"))
  equal(em_is_autosome(c("1", "chr22", "Chr2", "X", "chrY", "MT", "M", "23", "24", "26", "XY", "PAR1", NA)),
        c(TRUE, TRUE, TRUE, rep(FALSE, 10)))
  tissue <- function(p) {
    list(susie_mix = list(alpha = matrix(p, nrow = 1), pip = p, V = 1, converged = TRUE),
         weighted_fit_mix = list(alpha = matrix(p, nrow = 1), pip = p, V = 1, converged = TRUE),
         mix_coding = c("additive", "recessive", "dominant"))
  }
  auto <- list(Brain = tissue(c(.2, .3, .5)))
  saveRDS(auto, file.path(project, "results/AUTO.rds"))
  # An extreme successful X result must not affect priors. A corrupt Y result
  # must not be opened, and a missing MT result must not block preparation.
  saveRDS(list(Brain = tissue(c(1, 0, 0))), file.path(project, "results/GX.rds"))
  writeLines("unreadable", file.path(project, "results/GY.rds"))
  prepare <- function() em_prepare_iteration(project, gene_annotations = annotations)
  bad_annotations <- annotations[-1, ]
  expect_error(em_prepare_iteration(project, gene_annotations = bad_annotations), "Missing chromosome annotation")
  first <- prepare()
  p <- read.csv(file.path(first$iteration_dir, "priors.csv"))
  equal(em_prior_for_tissue(p, "Brain"), c(.2, .3, .5))
  stopifnot(p$n_source_files == 1, p$n_fits == 1, p$n_source_chromosome_exclusions == 3,
            p$n_source_gene_errors == 0, p$analysis_chromosomes == "autosomes_1_22")
  audit <- read.csv(file.path(first$iteration_dir, "source_audit.csv"))
  stopifnot(nrow(audit) == 3, all(audit$issue == "excluded_chromosome"))
  equal(first$chunks, 1:2) # Keep a chunk even if every gene is excluded.
  calls <- character()
  runner <- function(target_gene, ...) {
    calls <<- c(calls, target_gene)
    stopifnot(target_gene == "AUTO")
    auto
  }
  # Even an existing successful same-iteration X fit must be replaced by an
  # explicit exclusion before the ordinary saved-result reuse branch.
  x_old <- list(Brain = tissue(c(1, 0, 0)))
  attr(x_old, "em_iteration") <- p$iteration
  saveRDS(x_old, file.path(first$iteration_dir, "results/GX.rds"))
  s1 <- em_run_chunk(project, first$iteration_dir, 1L, runner)
  s2 <- em_run_chunk(project, first$iteration_dir, 2L, runner)
  equal(calls, "AUTO")
  equal(s1$status, c("finished", "excluded_chromosome"))
  stopifnot(all(s2$status == "excluded_chromosome"),
            file.exists(file.path(first$iteration_dir, "completed/chunk_002.done")))
  excluded <- readRDS(file.path(first$iteration_dir, "results/GX.rds"))
  stopifnot(identical(attr(excluded, "em_exclusion"), "non_autosomal_chromosome"))
  pooled <- em_estimate_coding_priors(file.path(first$iteration_dir, "results", c("AUTO.rds", "GX.rds")),
                                     "weighted_fit_mix")
  equal(em_prior_for_tissue(pooled$priors, "Brain"), c(.2, .3, .5))
  stopifnot(pooled$audit$issue == "excluded_chromosome")
  # Simulate migration from a completed pre-filter iteration: no chromosome
  # column in its frozen manifest, and a successful non-autosomal source fit.
  manifest_file <- file.path(first$iteration_dir, "manifest.csv")
  manifest <- read.csv(manifest_file)
  write.csv(manifest[c("chunk", "gene")], manifest_file, row.names = FALSE)
  saveRDS(x_old, file.path(first$iteration_dir, "results/GX.rds"))
  snapshot_before <- readLines(file.path(first$iteration_dir, "priors.csv"))
  manifest_before <- readLines(manifest_file)
  second <- prepare()
  p2 <- read.csv(file.path(second$iteration_dir, "priors.csv"))
  equal(em_prior_for_tissue(p2, "Brain"), c(.2, .3, .5))
  stopifnot(identical(snapshot_before, readLines(file.path(first$iteration_dir, "priors.csv"))),
            identical(manifest_before, readLines(manifest_file)))
  # Missing autosomal files still block advancement; exclusions never mask them.
  expect_error(em_check_result_files(file.path(project, "results"),
    data.frame(gene = c("AUTO", "MISSING")), allowed_extra_genes = c("GX", "GY")), "Missing: 1")
  # Direct workhorse calls also skip before touching expression or PLINK.
  # Only package loading and annotation are stubbed: the actual entry point runs.
  worker <- new.env(parent = globalenv())
  worker$library <- function(...) invisible(NULL)
  sys.source("script/scan_tissue_attempt/workhorse_em.R", worker)
  annotation_script <- file.path(project, "annotation.R")
  writeLines('get_gene_annotations <- function(...) data.frame(gene_name="GX", chromosome="chrX")', annotation_script)
  direct <- worker$run_susie_gene(target_gene = "GX", tissue_priors = p,
    project_dir = normalizePath(".", winslash = "/"), gene_annot_fun = annotation_script,
    subject_pheno_file = "MUST_NOT_BE_READ", expr_file = "MUST_NOT_BE_READ")
  stopifnot(identical(attr(direct, "em_exclusion"), "non_autosomal_chromosome"))
  cat("All EM chromosome exclusion tests passed.\n")
}
run_tests()
