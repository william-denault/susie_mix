# Rscript --vanilla script/analysis/tests/test_em_coloc_screen.R
# Synthetic saved fits, including unequal coding coverage. No GTEx/refitting needed.
source("script/analysis/em_analysis_utils.R")
source("script/analysis/em_coloc_screen_utils.R")
source("script/scan_tissue_attempt/em_utils.R")
root <- tempfile("em-coloc-test-")
dir.create(file.path(root, "script/scan_tissue_attempt"), recursive = TRUE)
dir.create(file.path(root, "script/analysis"), recursive = TRUE)
file.copy("script/scan_tissue_attempt/em_utils.R", file.path(root, "script/scan_tissue_attempt"))
file.copy(c("script/analysis/generate_summary_results.R", "script/analysis/em_analysis_utils.R",
            "script/analysis/em_coloc_screen_utils.R", "script/analysis/screen_em_coloc_candidates.R"),
          file.path(root, "script/analysis"))
dir.create(file.path(root, "results"))
snps <- paste0("chr1_", c(100000L, 105000L, 110000L, 130001L, 350000L, 450000L), "_A_G_b38_G")
add_map <- data.frame(predictor_index = 1:6, predictor_name = snps, snp = snps, coding = "additive")
mix_map <- do.call(rbind, lapply(c("additive", "recessive", "dominant"), function(coding)
  data.frame(predictor_name = paste(snps, coding, sep = "__"), snp = snps, coding = coding)))
mix_map <- mix_map[-c(8, 15), ] # Different numbers of predictors per coding/SNP.
mix_map$predictor_index <- seq_len(nrow(mix_map))
prior <- c(additive = .7, recessive = .1, dominant = .2)
weights <- em_predictor_weights(mix_map$coding, prior, mix_map$predictor_name)
make_fit <- function(map, weights, leads, members, components = seq_along(leads)) {
  alpha <- matrix(.01 / (nrow(map) - 1), max(components), nrow(map))
  for (i in seq_along(leads)) alpha[components[i], leads[i]] <- .99
  dimnames(alpha) <- list(NULL, map$predictor_name)
  # Each component's absolute BF has an explicit scale relative to its null.
  bf <- sweep(log(alpha), 2, log(weights), "-") + 20
  list(alpha = alpha, pip = 1 - apply(1 - alpha, 2, prod), lbf_variable = bf,
       pi = weights, converged = TRUE, elbo = c(-10, -5), KL = rep(1, nrow(alpha)),
       sets = list(cs = setNames(members, paste0("L", components)), cs_index = components))
}
add <- make_fit(add_map, rep(1/6, 6), c(1, 5), list(c(1L, 2L), 5L))
em_leads <- c(2, 3, which(mix_map$snp == snps[4] & mix_map$coding == "recessive"), 5, 6)
em <- make_fit(mix_map, weights, em_leads,
               list(2L, 3L, c(em_leads[3], 1L), 5L, 6L), c(3L, 1L, 5L, 2L, 4L))
baseline <- list(susie_add = add, susie_add_perm = add, susie_mix = em, susie_mix_perm = em,
                 add_predictor_map = add_map, mix_predictor_map = mix_map,
                 min_pv = 1e-10, n_ind = 30L, mean_read = 1, median_read = 1)
current <- list(weighted_fit_mix = em, mix_predictor_map = mix_map,
                weighted_mix_prior_weights = weights, mix_coding_prior = prior,
                n_ind = 30L, mean_read = 1, median_read = 1)
genes <- c("GOOD", "NOADD", "THRESHOLD", "BADBF", "BADLD", "NOEM", "MISSINGPOS")
manifest <- data.frame(chunk = 1L, gene = genes)
history <- NULL
for (iteration in 1:2) {
  d <- file.path(root, "results_em", sprintf("iteration_%03d", iteration))
  dir.create(file.path(d, "results"), recursive = TRUE)
  dir.create(file.path(d, "completed"))
  priors <- data.frame(iteration, tissue = "Tissue", pi_add = .7, pi_rec = .1, pi_dom = .2)
  history <- rbind(history, priors)
  write.csv(priors, file.path(d, "priors.csv"), row.names = FALSE)
  write.csv(manifest, file.path(d, "manifest.csv"), row.names = FALSE)
  if (iteration == 2) next
  for (gene in genes) {
    b <- baseline
    x <- current
    if (gene == "NOADD") b$susie_add$sets <- NULL
    if (gene == "THRESHOLD") b$min_pv <- 5e-8
    if (gene == "BADBF") x$weighted_fit_mix$lbf_variable <- NULL
    if (gene == "NOEM") x$weighted_fit_mix$sets <- NULL
    if (gene == "MISSINGPOS") {
      # Unparseable lead in every vanilla comparison prevents passing distance QC.
      b$add_predictor_map$snp[1] <- "rsUNKNOWN"
    }
    saveRDS(list(Tissue = b), file.path(root, "results", paste0(gene, ".rds")))
    out <- list(Tissue = x)
    attr(out, "em_iteration") <- 1L
    saveRDS(out, file.path(d, "results", paste0(gene, ".rds")))
  }
  saveRDS(list(iteration = 1L, chunk = 1L, n_genes = length(genes)), file.path(d, "completed/chunk_001.done"))
}
write.csv(history, file.path(root, "results_em/prior_history.csv"), row.names = FALSE)
set.seed(42)
g <- matrix(sample(0:2, 30*6, replace = TRUE), 30, dimnames = list(NULL, snps))
g[, 6] <- g[, 5] # Perfect LD across a 100-kb shift must not be called independent.
reads <- character()
reader <- function(gene, tissue, x) {
  reads <<- c(reads, gene)
  if (gene == "BADLD") stop("Synthetic LD failure")
  g
}
out <- suppressWarnings(em_screen_coloc(root, genotype_reader = reader))
stopifnot(all(out$all_em_cs$em_iteration == 1),
          !"THRESHOLD" %in% out$all_em_cs$gene, out$counts$n_no_em_cs == 1,
          all(out$candidates$min_distance_bp > 10000),
          !any(out$candidates$gene %in% c("NOADD", "MISSINGPOS")),
          nrow(out$no_vanilla_cs) == 5L,
          !any(reads %in% c("NOADD", "THRESHOLD", "NOEM", "MISSINGPOS")))
good <- subset(out$all_em_cs, gene == "GOOD")
stopifnot(identical(good$shifted, c(FALSE, FALSE, TRUE, FALSE, TRUE)),
          good$min_distance_bp[2] == 10000, good$min_distance_bp[4] == 0,
          abs(good$max_lead_r2[5] - 1) < 1e-12, good$max_shared_snps[3] == 1,
          subset(out$all_em_cs, gene == "BADBF")$bf_status[1] == "error",
          subset(out$all_em_cs, gene == "BADLD")$ld_status[1] == "error",
          all(c("bf_export", "ld") %in% out$audit$stage))
bundle <- readRDS(file.path(out$output_dir, good$bundle[1]))
stopifnot(identical(bundle$em$components, c(3L,1L,5L,2L,4L)),
          ncol(bundle$em$lbf) == 6L, nrow(bundle$vanilla$lbf) == 2L,
          all(abs(bundle$paired$em_retained_mass - 1) < 1e-10),
          colnames(bundle$em$lbf)[1] == "1:100000:A:G",
          length(unique(round(bundle$em$snp_prior, 6))) > 1)
# Sum over coding alternatives reproduces the joint evidence including unequal q_j.
for (i in seq_len(nrow(em$alpha))) {
  before <- sum(exp(em$lbf_variable[i, ]) * weights)
  after <- sum(exp(bundle$em$all_component_lbf[i, ]) * bundle$em$snp_prior)
  stopifnot(abs(before/after - 1) < 1e-12)
}
tab <- read.delim(gzfile(file.path(out$output_dir, sub("rds$", "tsv.gz", good$bundle[1]))))
stopifnot(nrow(tab) == (5+2)*6, all(is.finite(tab$log_bf)), !any(grepl("__", tab$snp)))
noadd <- readRDS(file.path(out$output_dir, out$no_vanilla_cs$bundle[1]))
stopifnot(nrow(noadd$vanilla$lbf) == 0L, nrow(noadd$vanilla$all_component_lbf) == 2L)
expect_error <- function(expr, pattern) {
  x <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(x, "error"), grepl(pattern, conditionMessage(x)))
}
bad <- em
bad$lbf_variable[1,1] <- bad$lbf_variable[1,1] + 3
expect_error(em_screen_bfs(bad, mix_map, weights), "inconsistent")
expect_error(em_screen_bfs(em, mix_map, rev(weights)), "names")
expect_error(em_screen_variants(c(snps[1], sub("_G$", "_A", snps[1]))), "Duplicate")
expect_error(em_screen_coloc(root, iteration = 2, compute_ld = FALSE), "unfinished")
expect_error(suppressWarnings(em_screen_coloc(root, compute_ld = FALSE, output_dir = out$output_dir)), "new or empty")
# The Source entry uses latest, does not need summaries, and marks unavailable LD.
saved_project <- Sys.getenv("SUSIE_MIX_PROJECT_DIR", unset = NA_character_)
Sys.setenv(SUSIE_MIX_PROJECT_DIR = root)
entry <- new.env(parent = globalenv())
suppressWarnings(source("script/analysis/screen_em_coloc_candidates.R", local = entry))
if (is.na(saved_project)) Sys.unsetenv("SUSIE_MIX_PROJECT_DIR") else Sys.setenv(SUSIE_MIX_PROJECT_DIR = saved_project)
stopifnot(nrow(entry$em_coloc_screen$candidates) == 6L,
          all(entry$em_coloc_screen$candidates$ld_status == "error"))
# No passing P values: empty CSVs retain usable headers, no genotype access.
empty <- suppressWarnings(em_screen_coloc(root, association_threshold = 1e-12, compute_ld = FALSE))
stopifnot(nrow(empty$candidates) == 0L, nrow(empty$bf_manifest) == 0L,
          "gene" %in% names(read.csv(file.path(empty$output_dir, "candidates.csv"))))
# Check the exporter against the installed SuSiE implementation, not just constructed BFs.
if (requireNamespace("susieR", quietly = TRUE)) {
  set.seed(451)
  gx <- matrix(rbinom(600*6, 2, .35), 600, dimnames = list(NULL, snps))
  mx <- cbind(gx, 1*(gx == 2), 1*(gx > 0))[, -c(8,15)]
  colnames(mx) <- mix_map$predictor_name
  yy <- 3*(gx[,4] == 2) + 2*gx[,6] + rnorm(600)
  real_fit <- susieR::susie(mx, yy, L = 2, prior_weights = weights, max_iter = 200)
  real_bf <- em_screen_bfs(real_fit, mix_map, weights)
  stopifnot(nrow(real_bf$lbf) == 2L, all(is.finite(real_bf$lbf)))
} else message("SKIP installed-SuSiE smoke test: susieR is unavailable")
cat("PASS: latest-complete selection, strict P/distance boundaries, all vanilla leads,\n",
    "CS overlap, tissue LD, BF marginalization, failures, no-CS cases, Source entry and empty output.\n", sep = "")
