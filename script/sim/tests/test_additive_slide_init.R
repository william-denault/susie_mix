# Small real fits, using the same installed packages as the simulation.
source("script/sim/sim_additive_slide_init.R")
source("script/sim/sim_workhorse.R")
check_init_packages()
near <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-8)))
expect_error <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))
test_dir <- file.path(init_project_dir, "tmp/additive_init_validation")
genotypes <- file.path(test_dir, "genotypes")
dir.create(genotypes, recursive = TRUE, showWarnings = FALSE)
# Fixed, small genotype fixture; no external data or live temp_plink required.
set.seed(9014)
X <- sapply(seq(.08, .45, length.out = 24), function(maf) rbinom(650, 2, maf))
colnames(X) <- paste0("chr1_", 100001:100024, "_A_G_b38_A")
raw <- data.frame(FID = 1:650, IID = paste0("donor", 1:650), PAT = 0, MAT = 0,
                  SEX = 0, PHENOTYPE = -9, X)
raw_file <- file.path(genotypes, "fixture.raw")
if (!file.exists(raw_file)) write.table(raw, raw_file, row.names = FALSE, quote = FALSE, sep = "\t")

data <- simulate_additive_init_data(1000001L, genotypes)
stopifnot(nrow(data$X) == 500, length(unique(data$true_pos)) == 2)
near(data$genetic_variance, .05)
near(data, simulate_additive_init_data(1000001L, genotypes))
# Execute the original generator up to its first fit, comparing X and y
# directly rather than relying on version-dependent optimizer defaults.
reference <- new.env()
for (nm in names(formals(sim_mix))) assign(nm, eval(formals(sim_mix)[[nm]]), reference)
reference$pve <- .05; reference$n <- 500; reference$L_add <- 2; reference$L_rec <- 0
reference$seed <- 1000001; reference$temp_dir <- genotypes
for (expr in as.list(body(sim_mix))[-1]) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("susie_res"))) break
  eval(expr, reference)
}
stopifnot(identical(data$causal_snps, reference$causal_snps),
          identical(data$true_pos, reference$true_pos))
near(data$X, reference$geno_all); near(data$y, reference$y)
near(data$beta_standardized, reference$beta)
fitted <- fit_additive_init(data)
# Independently reproduce the requested direct warm-start call.
args <- list(X = data$X, y = data$y, L = 10, standardize = TRUE,
  estimate_prior_method = "optim", coverage = .95, min_abs_corr = .5,
  max_iter = 1000, tol = .001, verbose = FALSE)
cold <- do.call(susieR::susie, args)
near(fitted$fits$SuSiE$pip, cold$pip)
args[[check_init_packages()$init_arg]] <- fitted$fits[["SuSiE-slide"]]
direct <- do.call(susieR::susie, args)
near(fitted$fits[["SuSiE-init-slide"]]$pip, direct$pip)
near(fitted$fits[["SuSiE-init-slide"]]$elbo, direct$elbo)
compact <- compact_additive_init(data, fitted, .05)
stopifnot(identical(compact$metrics$method, init_methods),
          all(compact$metrics$n_causal == 2), !is.null(compact$slide_delta_causal),
          !"full_fits" %in% names(compact))

out <- file.path(test_dir, "results")
path <- run_additive_initialization(1L, 2L, genotypes, out)
digest <- tools::md5sum(path)
run_additive_initialization(1L, 2L, genotypes, out)
stopifnot(identical(digest, tools::md5sum(path)))
expect_error(run_additive_initialization(1L, 2L, genotypes, out, pve = .1))
# A saved failure is retried without rerunning successful records.
saved <- readRDS(path)
first <- saved$results[[1]]
saved$results[[2]] <- list(seed = 1000002L, error = "simulated interruption")
saveRDS(saved, path)
run_additive_initialization(1L, 2L, genotypes, out)
recovered <- readRDS(path)
stopifnot(identical(first, recovered$results[[1]]), is.null(recovered$results[[2]]$error))
summary <- summarize_additive_initialization(out)
stopifnot(nrow(summary$metrics) == 6, nrow(summary$summary) == 3,
          nrow(summary$optimization) == 2, nrow(summary$errors) == 0,
          all(summary$summary$n_replicates == 2))
near(summary$optimization$additive_elbo_gain,
  summary$optimization$initialized_additive_elbo - summary$optimization$additive_elbo)
cat("PASS: paired additive generator, direct slide initialization, compact results, checkpoint resume/retry, and summaries.\n")
