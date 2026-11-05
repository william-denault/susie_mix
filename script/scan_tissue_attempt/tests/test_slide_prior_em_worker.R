# Exercise the real GTEx worker with tiny synthetic input files and mocked
# PLINK extraction. The actual data preparation and slider fitting run.
extra_lib <- Sys.getenv("SLIDER_PRIOR_TEST_LIB")
if (nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))
for (file in c("em_utils.R", "slide_prior_em_utils.R")) source(file.path("script/scan_tissue_attempt", file))
spe_check_packages()
worker_env <- new.env(parent = globalenv())
sys.source("script/scan_tissue_attempt/workhorse_slide_prior_em.R", worker_env)
root <- tempfile("slider_worker_test_"); dir.create(root)
set.seed(145)
n <- 160L
ids <- sprintf("GTEX-%05d", 1:n)
samples <- paste0(ids, "-SM-TEST")
X <- matrix(rbinom(n * 4, 2, .35), n, 4)
X[, 4] <- rep(c(0, 0, 1, 0), length.out = n)
colnames(X) <- paste0("chr1_", 1001:1004, "_A_C_b38_A")
write.table(data.frame(SUBJID = ids, SEX = rep(1:2, n / 2), AGE = "40-49"),
            file.path(root, "subjects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(SAMPID = samples, SMTS = rep(c("Brain", "Liver"), each = n / 2),
                       SMTSD = "test", SMGEBTCHT = "batch", SMAFRZE = "RNASEQ"),
            file.path(root, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
reads <- round(exp(5 + .7 * (X[, 1] + .5 * (X[, 1] == 1)) + rnorm(n, sd = .25)))
expression <- data.frame(Name = c("E1", "E2"), Description = c("GENE", "HOUSE"),
                         rbind(reads, rep(1e6, n)), check.names = FALSE)
names(expression)[-(1:2)] <- samples
writeLines(c("#1.2", paste(2, n, sep = "\t")), file.path(root, "expression.gct"))
suppressWarnings(write.table(expression, file.path(root, "expression.gct"), sep = "\t",
                             quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE))
raw <- data.frame(FID = ids, IID = ids, PAT = 0, MAT = 0, SEX = 1, PHENO = -9, X,
                   check.names = FALSE)
write.table(raw, file.path(root, "genotypes.raw"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines('get_gene_annotations <- function(path) data.frame(gene_name="GENE", chromosome="1", start=1000, end=2000, strand="+")',
           file.path(root, "annotations.R"))
worker_env$system <- function(command) {
  prefix <- sub(".* --out ", "", command)
  stopifnot(file.copy(file.path(root, "genotypes.raw"), paste0(prefix, ".raw")))
  0L
}
w <- matrix(1 / 17, 2, 17, dimnames = list(NULL, spe_columns()))
priors <- data.frame(tissue = c("Brain", "Liver"), w)
args <- list(target_gene = "GENE", tissue_priors = priors,
             project_dir = normalizePath(".", winslash = "/"), datadir = root,
             gene_annot_fun = file.path(root, "annotations.R"),
             subject_pheno_file = file.path(root, "subjects.tsv"),
             sample_attr_file = file.path(root, "samples.tsv"),
             expr_file = file.path(root, "expression.gct"),
             temp_dir = file.path(root, "plink"), L = 1L)
out <- do.call(worker_env$run_slide_prior_gene, args)
stopifnot(identical(names(out), c("Brain", "Liver")), !length(attr(out, "tissue_errors")))
for (tissue in names(out)) {
  f <- out[[tissue]]$fit_slide_prior
  counts <- spe_fit_counts(f, spe_prior(priors, tissue))
  stopifnot(out[[tissue]]$n_ind == 80, out[[tissue]]$n_SNP == 4,
            f$delta_forced[4], !out[[tissue]]$em_warm_started,
            is.finite(out[[tissue]]$min_pv))
  priors[priors$tissue == tissue, spe_columns()] <- spe_update_weights(colSums(counts), spe_prior(priors, tissue))
}
args$tissue_priors <- priors; args$previous_result <- out
next_out <- do.call(worker_env$run_slide_prior_gene, args)
for (tissue in names(out)) {
  f <- next_out[[tissue]]$fit_slide_prior
  spe_fit_counts(f, spe_prior(priors, tissue))
  stopifnot(next_out[[tissue]]$em_warm_started,
            identical(out[[tissue]]$sample_ids, next_out[[tissue]]$sample_ids))
}
stopifnot(!length(list.files(args$temp_dir)))
cat("PASS: synthetic GTEx import, donor alignment, QC, association metadata, real discrete fitting and warm-start refit.\n")

# Integrate the recovery wrapper with real fits at the frozen L=10 setting.
# The package clips L to p=4 here, which the recovery validator must accept.
for (file in c("slide_prior_recovery.R", "prepare_slide_prior_em_iteration.R"))
  source(file.path("script/scan_tissue_attempt", file))
project <- file.path(root, "r")
d <- file.path(project, "results_slide_prior_em/iteration_000")
dir.create(d, recursive = TRUE)
scripts <- file.path(project, "script/scan_tissue_attempt")
dir.create(scripts, recursive = TRUE)
invisible(file.copy("script/scan_tissue_attempt/slide_prior_recovery.R", scripts))
frozen <- data.frame(iteration = 0L, source_iteration = -1L, tissue = c("Brain", "Liver"), w,
                     n_fits = 0L, n_active = 0L, free_mass = 0, prior_max_change = 0,
                     m_step_q_gain = 0, source_elbo = NA_real_,
                     update_method = "slider_prior_active_counts_v1", created_at = "fixture")
em_atomic_write(frozen, file.path(d, "priors.csv"), csv = TRUE)
em_atomic_write(frozen, file.path(project, "results_slide_prior_em/prior_history.csv"), csv = TRUE)
em_atomic_write(data.frame(chunk = 1L, gene = "GENE"), file.path(d, "manifest.csv"), csv = TRUE)
em_atomic_write(data.frame(delta = spe_grid()), file.path(d, "grid.csv"), csv = TRUE)
em_atomic_write(c(spe_settings(), list(package_version = spe_check_packages())), file.path(d, "settings.rds"))
run_real <- function(...) {
  supplied <- list(...)
  paths <- args[c("datadir", "gene_annot_fun", "subject_pheno_file", "sample_attr_file", "expr_file")]
  supplied[names(paths)] <- paths
  supplied$project_dir <- normalizePath(".", winslash = "/")
  do.call(worker_env$run_slide_prior_gene, supplied)
}
ctx <- sre_context(project, d, check_package = TRUE)
sre_audit(ctx)
# Exercise the real GTF reader used to freeze chromosome scope on RCC.
invisible(file.copy("script/scan_tissue_attempt/get_gene_annotations.R", scripts))
scope_gtf <- file.path(root, "scope.gtf")
writeLines(c("# fixture GTF", paste("chr1", "ensembl", "gene", 1000, 2000, ".", "+", ".",
  'gene_id "E1"; gene_type "protein_coding"; gene_name "GENE";', sep = "\t")), scope_gtf)
ctx <- sre_exclude_sex_mt(ctx, gtf_file = scope_gtf)
stopifnot(!length(ctx$scope$excluded), identical(unname(ctx$scope$chromosomes), "1"))
sre_run_gene(ctx, "GENE", run_real)
sre_verify(ctx)
launch <- spe_prepare_iteration(project, "new")
ctx1 <- sre_context(project, launch$iteration_dir, check_package = TRUE)
sre_audit(ctx1)
sre_run_gene(ctx1, "GENE", run_real)
sre_verify(ctx1)
saved <- readRDS(sre_result(ctx1, "GENE"))
stopifnot(all(vapply(saved, function(x) isTRUE(x$em_warm_started), logical(1))))
cat("PASS: recovery wrapper, actual L clipping, completed-iteration gate, original M-step and real warm-start recovery.\n")
