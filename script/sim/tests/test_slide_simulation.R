# Run from the project root with Rscript. Uses installed fitting packages.
source("script/sim/write_jobs.R")
source("script/sim/run_job.R")
root <- normalizePath(".", winslash = "/")
validation <- file.path(root, "tmp/slide_simulation_validation")
dir.create(validation, recursive = TRUE, showWarnings = FALSE)

expect_error <- function(expression, pattern) {
  result <- tryCatch(force(expression), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}

# Independent combinatorial enumeration: choose up to three types, then assign
# positive counts summing to each K. No absent or four-type allocation is legal.
expected <- list()
for (types in 1:3) for (ids in combn(5, types, simplify = FALSE)) {
  allocations <- expand.grid(rep(list(1:5), types))
  allocations <- allocations[rowSums(allocations) <= 5, , drop = FALSE]
  for (i in seq_len(nrow(allocations))) {
    counts <- integer(5)
    counts[ids] <- as.integer(allocations[i, ])
    expected[[length(expected) + 1L]] <- sim_configuration(counts)
  }
}
conditions <- sim_conditions()
observed <- apply(conditions[sim_count_columns], 1, sim_configuration)
stopifnot(nrow(conditions) == 225L, !anyDuplicated(observed),
          setequal(observed, unlist(expected)), length(unique(conditions$name)) == 25L,
          identical(as.integer(table(conditions$K)), c(5L, 15L, 35L, 65L, 105L)))
legacy <- read.csv("script/sim/jobs/conditions.csv")
stopifnot(nrow(legacy) == 55L,
          all(apply(cbind(legacy[c("L_add", "L_rec", "L_dom")], L_prec = 0, L_pdom = 0),
                    1, sim_configuration) %in% observed))
expect_error(sim_scenario_name(c(1, 1, 1, 1, 0)), "sum")
expect_error(sim_scenario_name(c(0, 0, 0, 0, 0)), "sum")

manifest <- write_simulation_jobs(validation)
stopifnot(nrow(manifest) == 1125L, sum(manifest$reps_per_chunk) == 450000,
          identical(sort(unique(manifest$pve)), c(.05, .1, .2, .3, .4)),
          !anyDuplicated(manifest$output_file))
parsed <- sim_parse_checkpoints(manifest$output_file)
stopifnot(all(as.matrix(parsed[sim_count_columns]) == as.matrix(manifest[sim_count_columns])),
          all(parsed$pve == manifest$pve))
old_name <- "recessive_add0_rec1_dom0_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData"
stopifnot(all(sim_parse_checkpoints(old_name)[c("L_prec", "L_pdom")] == 0))

# Independent hard-call fixture, with enough common SNPs for every effect type.
genotypes <- file.path(validation, "genotypes")
dir.create(genotypes, recursive = TRUE, showWarnings = FALSE)
set.seed(9014)
X <- sapply(seq(.06, .45, length.out = 24), function(maf) rbinom(650, 2, maf))
colnames(X) <- paste0("chr1_", 100001:100024, "_A_G_b38_A")
raw <- data.frame(FID = 1:650, IID = paste0("donor", 1:650), PAT = 0, MAT = 0,
                  SEX = 0, PHENOTYPE = -9, X)
write.table(raw, file.path(genotypes, "fixture.raw"), row.names = FALSE, quote = FALSE, sep = "\t")

# One real three-method fit for all 25 supports, plus an unequal K=5 triple
# and a high-PVE repeat to check seed pairing and the noise scaling convention.
selected <- which(!duplicated(conditions$name))
selected <- c(selected, which(observed == "add0_rec2_dom0_prec1_pdom2"))
smoke_manifest <- manifest[manifest$pve == .05 &
                            manifest$name %in% conditions$name[selected], ]
smoke_manifest <- smoke_manifest[match(observed[selected],
    apply(smoke_manifest[sim_count_columns], 1, sim_configuration)), ]
smoke_manifest <- rbind(smoke_manifest,
  manifest[manifest$pve == .4 & manifest$L_add == 0 & manifest$L_rec == 2 &
             manifest$L_dom == 0 & manifest$L_prec == 1 & manifest$L_pdom == 2, ])
chunks <- file.path(validation, "smoke_chunks")
dir.create(chunks, recursive = TRUE, showWarnings = FALSE)
fitted <- list()
for (i in seq_len(nrow(smoke_manifest))) {
  cell <- smoke_manifest[i, ]
  args <- as.list(cell[c("pve", "n", "L", sim_count_columns)])
  args$seed <- 1000001
  args$temp_dir <- genotypes
  x <- suppressMessages(do.call(sim_mix, args))
  stopifnot(identical(x$metrics$method, c("SuSiE", "SuSiE-mix", "SuSiE-slide")),
            abs(x$genetic_variance - cell$pve) < 1e-12,
            length(unique(x$causal_snps)) == cell$K,
            identical(x$causal_delta, unname(sim_effect_delta[x$causal_coding])),
            identical(is.na(x$true_pos_mix), grepl("^partial_", x$causal_coding)),
            all(abs(x$susie_slide_delta_causal) <= 1),
            ncol(x$susie_slide_delta_causal) == cell$K,
            length(x$susie_slide_pip) == length(x$susie_pip),
            all(x$susie_slide_pip >= 0 & x$susie_slide_pip <= 1),
            !any(c("y", "g", "noise") %in% names(x)))
  full_indices <- !is.na(x$true_pos_mix)
  stopifnot(identical(x$mix_to_add[x$true_pos_mix[full_indices]], x$true_pos[full_indices]))
  for (j in seq_along(x$susie_slide_cs$cs)) {
    cs <- x$susie_slide_cs$cs[[j]]
    stopifnot(all(cs %in% seq_along(x$susie_pip)))
  }
  results <- list(x)
  save(results, file = file.path(chunks, basename(cell$output_file)))
  fitted[[i]] <- x
  cat("Validated support", i, "of", nrow(smoke_manifest), ":", cell$name, "\n")
}
low <- fitted[[length(fitted) - 1L]]
high <- fitted[[length(fitted)]]
stopifnot(identical(low$causal_snps, high$causal_snps),
          isTRUE(all.equal(high$beta_standardized / low$beta_standardized, rep(sqrt(8), 5))))

# Reconstruct partial predictors independently, including minor-allele orientation
# and donor sampling, and verify the original additive fit saw this phenotype.
set.seed(1000001)
invisible(sample.int(1, 1))
G <- qc_filter_geno(X)$X
G <- G[sample.int(nrow(G), 500), , drop = FALSE]
G <- G[, colSums(G) >= 5 & matrixStats::colSds(G) > 0, drop = FALSE]
counts <- sapply(0:2, function(g) colSums(G == g))
eligible <- colnames(G)[apply(counts, 1, min) >= 5]
snps <- eligible[sample.int(length(eligible), 5)]
stopifnot(identical(snps, low$causal_snps))
truth_X <- G[, snps, drop = FALSE]
# Effects are recessive, recessive, partial recessive, partial dominant, partial dominant.
truth_X[, 1:2] <- (truth_X[, 1:2] == 2) * 1
for (j in 3:5) truth_X[, j] <- truth_X[, j] + c(-.5, .5, .5)[j - 2] * (truth_X[, j] == 1)
signs <- sample(c(-1, 1), 5, replace = TRUE)
g <- drop(scale(truth_X) %*% signs)
g <- g * sqrt(.05 / var(g))
y <- g + rnorm(500, sd = sqrt(.95))
storage.mode(G) <- "double"
reference <- susieR::susie(G, y, L = 10, standardize = TRUE, estimate_prior_method = "optim",
                          coverage = .95, min_abs_corr = .5, max_iter = 1000)
# Algebraically equivalent rescaling can differ at roundoff level; the prior
# variance optimizer amplifies that slightly. Check the resulting PIPs to 1e-7.
stopifnot(max(abs(unname(reference$pip) - low$susie_pip)) < 1e-7)

# Check actual save/resume through the shared runner, including refusal to reuse
# an incompatible checkpoint. Use two reps only, in this isolated test project.
runner_root <- file.path(validation, "runner")
runner_manifest <- write_simulation_jobs(runner_root, pve_values = .05, reps_per_chunk = 2)
job_id <- runner_manifest$job_id[runner_manifest$L_prec == 1 & runner_manifest$K == 1]
checkpoint <- run_simulation_job(job_id, runner_root, genotypes)
saved <- new.env()
load(checkpoint, saved)
stopifnot(length(saved$results) == 2, all(vapply(saved$results, function(x) is.null(x$error), TRUE)))
original <- saved$results
# Truncate a completed checkpoint to emulate interruption after its first save.
saved$results <- saved$results[1]
save(list = c("results", "checkpoint_settings"), envir = saved, file = checkpoint)
run_simulation_job(job_id, runner_root, genotypes)
load(checkpoint, saved)
stopifnot(identical(saved$results, original))
before <- tools::md5sum(checkpoint)
run_simulation_job(job_id, runner_root, genotypes)
stopifnot(identical(before, tools::md5sum(checkpoint)))
saved$checkpoint_settings$design$delta_prec <- -.25
save(list = c("results", "checkpoint_settings"), envir = saved, file = checkpoint)
expect_error(run_simulation_job(job_id, runner_root, genotypes), "Checkpoint design")
saved$checkpoint_settings$design$delta_prec <- -.5
save(list = c("results", "checkpoint_settings"), envir = saved, file = checkpoint)

# Read all 25 scenario families through the plotting pipeline. The figure smoke
# exports one representative from each metric and all five scenario groups.
script <- parse("script/sim/plot_simulations.R")
preview <- file.path(validation, "figures")
overrides <- list(project_dir = root, chunk_dir = chunks, output_dir = preview,
                  bootstrap_reps = 50, write_png = FALSE)
e <- new.env()
for (expression in script) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) && is.symbol(expression[[2]])) {
    name <- as.character(expression[[2]])
    if (name %in% names(overrides)) expression[[3]] <- overrides[[name]]
  }
  # Loading/aggregation and plotting functions run; the large production export
  # loops are replaced below by a small representative selection.
  if (is.call(expression) && identical(expression[[1]], as.name("for")) &&
      identical(expression[[2]], as.name("metric"))) next
  if (is.call(expression) && identical(expression[[1]], as.name("save_roc_figures"))) next
  eval(expression, e)
}
stopifnot(sum(e$audit$included) == length(fitted), nrow(e$replicates) == 3L * length(fitted),
          length(unique(e$replicates$scenario)) == 25L,
          nrow(e$metric_analysis$differences) == 12L * length(fitted))
for (group in names(e$scenario_groups))
  e$save_figure("coverage", e$scenario_groups[[group]], paste0("coverage_", group))
for (metric in c("purity", "power", "cs_size", "roc", "power_fdr"))
  e$save_figure(metric, e$scenario_groups$triples_2, paste0(metric, "_triples_2"),
                only_K = if (metric %in% c("roc", "power_fdr")) 5 else NULL)
e$write_png <- TRUE
e$save_figure("power", e$scenario_groups$triples_2, "review_power_triples")

# All three paired contrasts use the same seed blocks; a fixed recovery gap
# must have the same point and interval even when configurations are repeated.
toy <- data.frame(scenario = "test", pve = .05, K = 5, configuration = "one", seed = 1:20,
                   method = "SuSiE", n_cs = 1, covered_cs = rep(0:1, 10),
                   purity_sum = .8, cs_size_sum = 2, recovered = rep(0:1, 10), n_causal = 5)
mix <- slide <- toy
mix$method <- "SuSiE-mix"; mix$recovered <- toy$recovered + 1
slide$method <- "SuSiE-slide"; slide$recovered <- toy$recovered + 2
paired <- rbind(toy, mix, slide)
duplicated_config <- paired; duplicated_config$configuration <- "two"
a <- e$summarize_metrics(paired, B = 100, seed = 42)
b <- e$summarize_metrics(rbind(paired, duplicated_config), B = 100, seed = 42)
contrasts <- a$differences[a$differences$metric == "power", ]
stopifnot(isTRUE(all.equal(contrasts$difference, c(.2, .4, .2))),
          max(abs(contrasts$lower - contrasts$difference)) < 1e-12,
          max(abs(contrasts$upper - contrasts$difference)) < 1e-12,
          isTRUE(all.equal(a$differences, b$differences)))
expect_error(e$summarize_metrics(paired[-1, ], B = 20), "same replicates")
empty <- e$cs_summary(list(cs = NULL, purity = NULL), list(), c(1, 2))
stopifnot(empty["n_cs"] == 0, empty["recovered"] == 0, empty["n_causal"] == 2)
t <- c(0, .1, .5, .9, 1)
p <- c(0, .1, .1, .9, 1)
pc <- e$pip_counts(p, c(2, 5), t)
stopifnot(all(pc[, "tp"] == sapply(c(t, Inf), function(cut) sum(p[c(2, 5)] >= cut))),
          all(pc[, "fp"] == sapply(c(t, Inf), function(cut) sum(p[-c(2, 5)] >= cut))))
cat("PASS: complete grid, 27 three-method datasets, generating effects/PVE, resume, plots and paired metrics.\n")
