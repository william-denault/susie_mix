# Resume a completed v2 run without CS metadata or changing its saved state.
# Uses synthetic files only; no GTEx input or cluster jobs are involved.
source("script/scan_tissue_attempt/em_utils.R")
source("script/scan_tissue_attempt/prepare_em_iteration.R")

run_unfiltered_resume_tests <- function() {
  equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))
  dir.create("tmp", showWarnings = FALSE)
  base <- normalizePath("tmp", winslash = "/", mustWork = TRUE)
  project <- tempfile("em_unfiltered_resume_", tmpdir = base)
  dir.create(file.path(project, "results_em"), recursive = TRUE)
  on.exit({
    resolved <- normalizePath(project, winslash = "/", mustWork = TRUE)
    stopifnot(startsWith(resolved, paste0(base, "/em_unfiltered_resume_")))
    unlink(resolved, recursive = TRUE)
  }, add = TRUE)

  coding <- c("additive", "recessive", "dominant")
  tissue <- function(alpha, V) {
    alpha <- matrix(alpha, ncol = 3, byrow = TRUE)
    active <- rep_len(V, nrow(alpha)) > 0
    list(weighted_fit_mix = list(alpha = alpha, V = V,
      pip = 1 - apply(1 - alpha[active, , drop = FALSE], 2, prod),
      converged = TRUE, elbo = -10, sets = list(cs = NULL)), mix_coding = coding,
      n_ind = 54L, min_pv = 1, mean_read = 1)
  }
  fits <- list(
    G1 = list(Tissue = tissue(c(.5, .3, .2, .2, .3, .5, .98, .01, .01), c(1, 1, 0))),
    G2 = list(Tissue = tissue(c(.1, .2, .7), 1e-12)),
    G3 = list(Tissue = tissue(rep(1/3, 18), 0)))
  manifest <- data.frame(chunk = c(1L, 1L, 2L), gene = names(fits))
  history <- NULL
  for (i in 1:12) {
    directory <- file.path(project, "results_em", sprintf("iteration_%03d", i))
    dir.create(directory)
    prior <- data.frame(iteration = i, source_iteration = i - 1L,
      source_fit = if (i == 1L) "susie_mix" else "weighted_fit_mix",
      tissue = "Tissue", pi_add = .6, pi_rec = .3, pi_dom = .1,
      update_method = "susie_active_component_alpha_v2")
    write.csv(prior, file.path(directory, "priors.csv"), row.names = FALSE)
    write.csv(manifest, file.path(directory, "manifest.csv"), row.names = FALSE)
    history <- rbind(history, prior)
  }
  history_file <- file.path(project, "results_em/prior_history.csv")
  write.csv(history, history_file, row.names = FALSE)
  previous_dir <- file.path(project, "results_em/iteration_012")
  dir.create(file.path(previous_dir, "results"))
  dir.create(file.path(previous_dir, "completed"))
  for (gene in names(fits)) saveRDS(fits[[gene]], file.path(previous_dir, "results", paste0(gene, ".rds")))
  for (chunk in 1:2) saveRDS(list(chunk = chunk),
    file.path(previous_dir, "completed", sprintf("chunk_%03d.done", chunk)))
  frozen <- list.files(file.path(project, "results_em"), recursive = TRUE, full.names = TRUE)
  frozen <- setdiff(frozen, history_file)
  hashes <- tools::md5sum(frozen)

  # All positive-variance rows count despite no reported CS, weak association,
  # low read count, small sample size, and a tiny positive prior variance.
  next_iteration <- em_prepare_iteration(project)
  stopifnot(basename(next_iteration$iteration_dir) == "iteration_013")
  equal(next_iteration$chunks, 1:3)
  prior <- read.csv(file.path(next_iteration$iteration_dir, "priors.csv"))
  equal(em_prior_for_tissue(prior, "Tissue"), c(.8, .8, 1.4)/3)
  stopifnot(prior$source_iteration == 12L, prior$alpha_total == 3,
            prior$n_components == 10L, prior$n_active_components == 3L,
            prior$n_zero_variance_components == 7L,
            prior$update_method == "susie_active_component_alpha_v2")
  equal(tools::md5sum(frozen), hashes)
  updated_history <- read.csv(history_file)
  equal(updated_history[updated_history$iteration <= 12L, names(history)], history)
  equal(read.csv(file.path(next_iteration$iteration_dir, "manifest.csv"))$gene, manifest$gene)

  # Old filtered metadata, including an empty eligibility set, cannot alter EM.
  source_files <- file.path(previous_dir, "results", paste0(names(fits), ".rds"))
  original <- em_estimate_coding_priors(source_files, "weighted_fit_mix")
  fits$G1$Tissue$weighted_fit_mix$coding_prior_cs <- list(components = integer())
  fits$G2$Tissue$weighted_fit_mix$sets <- list(cs = list(L1 = 1:3),
    cs_index = 1L, purity = data.frame(min.abs.corr = 0.1))
  for (gene in names(fits)) saveRDS(fits[[gene]], file.path(previous_dir, "results", paste0(gene, ".rds")))
  equal(em_estimate_coding_priors(source_files, "weighted_fit_mix"), original)

  # Retry the pending iteration without changing its priors or repartitioning.
  new_files <- list.files(next_iteration$iteration_dir, recursive = TRUE, full.names = TRUE)
  new_hashes <- tools::md5sum(new_files)
  retry <- em_prepare_iteration(project, "resume", n_chunks = 1L)
  equal(retry, next_iteration)
  equal(tools::md5sum(new_files), new_hashes)
  cat("Unfiltered continuation 12 -> 13 passed; old files and history rows preserved.\n")
}
run_unfiltered_resume_tests()
