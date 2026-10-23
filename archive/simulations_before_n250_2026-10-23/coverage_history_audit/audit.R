root <- normalizePath(getwd(), winslash = "/")
out <- file.path(root, "tmp/coverage_history_audit")
old <- read.csv("simulation results/figures/metric_summary.csv")
new <- read.csv("simulation results/slide_prior_v1/figures/metric_summary.csv")
fields <- c("pve", "K", "method", "coverage", "n_cs", "covered_cs", "power", "n_replicates")
a <- old[old$scenario == "Additive only", fields]
b <- new[new$scenario == "Additive only" & new$method %in% c("SuSiE", "SuSiE-mix"), fields]
comparison <- merge(a, b, by = c("pve", "K", "method"), suffixes = c("_old", "_new"), all = TRUE)
comparison <- comparison[order(comparison$pve, comparison$K, comparison$method), ]
write.csv(comparison, file.path(out, "additive_comparison.csv"), row.names = FALSE)
print(comparison, row.names = FALSE, digits = 5)
paired <- read.csv("simulation results/slide_prior_v1/figures/method_comparison.csv")
paired <- paired[paired$scenario == "Additive only" & paired$metric == "coverage" & paired$reference == "SuSiE" & paired$method == "SuSiE-mix", ]
cat("\nCURRENT PAIRED COVERAGE DIFFERENCES\n")
print(paired[, c("pve", "K", "difference", "lower", "upper")], row.names = FALSE, digits = 5)
write.csv(paired, file.path(out, "current_paired_coverage.csv"), row.names = FALSE)
load_one <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  cat("\nFILE:", path, "\nOBJECTS:", paste(ls(e), collapse = ", "), "\n")
  if (exists("checkpoint_settings", e, inherits = FALSE)) {
    s <- e$checkpoint_settings
    cat("CHECKPOINT SETTINGS (except full gene pool and input file inventory):\n")
    str(s[setdiff(names(s), "genotype_inputs")], max.level = 3)
    if (!is.null(s$genotype_inputs)) {
      print(s$genotype_inputs[setdiff(names(s$genotype_inputs), c("files", "loci"))])
      cat("GENE POOL SIZE:", nrow(s$genotype_inputs$loci), "\n")
    }
  }
  r <- e$results
  ok <- vapply(r, function(x) is.null(x$error), logical(1))
  cat("N RESULTS:", length(r), "SUCCESS:", sum(ok), "\n")
  z <- r[[which(ok)[1]]]
  cat("FIRST SUCCESS FIELDS:", paste(names(z), collapse = ", "), "\n")
  str(z[c("seed", "settings", "raw_file", "gene", "genotype_mode", "region", "causal_snps", "true_pos", "metrics")], max.level = 3)
  meta <- do.call(rbind, lapply(r[ok], function(x) data.frame(
    seed = x$seed, raw_file = x$raw_file,
    gene = if (is.null(x$gene)) sub("\\.raw$", "", basename(x$raw_file)) else x$gene,
    causal_snps = paste(x$causal_snps, collapse = ";"),
    p = length(x$susie_pip), mix_p = length(x$susie_mix_pip))))
  cat("DISTINCT REGIONS:", length(unique(meta$gene)), "\n")
  cat("FIRST FIVE REGIONS:\n"); print(head(meta, 5), row.names = FALSE)
  cat("P QUANTILES:\n"); print(quantile(meta$p))
  meta
}
old_file <- "simulation results/chunks/additive_add1_rec0_dom0_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData"
new_file <- "simulation results/slide_prior_v1/chunks/additive_add1_rec0_dom0_prec0_pdom0_n500_L10_pve0.1_seed1000000_reps400_chunk1.RData"
om <- load_one(old_file)
nm <- load_one(new_file)
mm <- merge(om, nm, by = "seed", suffixes = c("_old", "_new"))
cat("\nMATCHED SEEDS:", nrow(mm), "IDENTICAL CAUSAL SNPS:", sum(mm$causal_snps_old == mm$causal_snps_new), "IDENTICAL P:", sum(mm$p_old == mm$p_new), "\n")
write.csv(mm, file.path(out, "matched_seed_inputs_K1_PVE10.csv"), row.names = FALSE)
