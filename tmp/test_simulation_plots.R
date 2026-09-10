# Limited smoke test: at most 3 replicates from each of 9 selected files.
# The full 221-file analysis is deliberately not run.
script <- parse("script/sim/plot_simulations.R")
root <- normalizePath(".", winslash = "/")
preview <- file.path(root, "tmp/simulation_plot_smoke")
selected <- c(
  "additive_add1_rec0_dom0_n500_L10_pve0.2_seed1e+06_reps200_chunk1.RData",
  "additive_add1_rec0_dom0_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData",
  "dominant_add0_rec0_dom2_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData",
  "recessive_add0_rec3_dom0_n500_L10_pve0.3_seed1e+06_reps400_chunk1.RData",
  "additive_dominant_add1_rec0_dom1_n500_L10_pve0.2_seed1e+06_reps400_chunk1.RData",
  "additive_recessive_add1_rec1_dom0_n500_L10_pve0.3_seed1e+06_reps400_chunk1.RData",
  "recessive_dominant_add0_rec1_dom1_n500_L10_pve0.4_seed1e+06_reps400_chunk1.RData",
  "additive_recessive_dominant_add1_rec1_dom1_n500_L10_pve0.1_seed1e+06_reps400_chunk1.RData",
  "additive_add2_rec0_dom0_n500_L10_pve0.4_seed1e+06_reps400_chunk1.RData"
)
overrides <- list(project_dir = root, output_dir = preview,
                  max_reps_per_file = 3,
                  file_pattern = paste0("^(", paste(gsub(".", "[.]", selected, fixed = TRUE),
                                                    collapse = "|"), ")$"))
# Escape the literal '+' in the seed field for the file-selection regex.
overrides$file_pattern <- gsub("+", "[+]", overrides$file_pattern, fixed = TRUE)
for (i in seq_along(script)) {
  z <- script[[i]]
  if (is.call(z) && identical(z[[1]], as.name("<-")) && is.symbol(z[[2]])) {
    name <- as.character(z[[2]])
    if (name %in% names(overrides)) script[[i]][[3]] <- overrides[[name]]
  }
}
e <- new.env()
eval(script, envir = e)
stopifnot(nrow(e$audit) == 9, sum(e$audit$examined) == 27,
          sum(e$audit$included) == 24, sum(e$audit$duplicates) == 3,
          nrow(e$replicates) == 48,
          all(e$summary_table$coverage >= 0 & e$summary_table$coverage <= 1),
          all(e$summary_table$power >= 0 & e$summary_table$power <= 1),
          all(e$summary_table$purity >= .5 & e$summary_table$purity <= 1))

# Check thresholds, ties, and zero/one endpoints against direct counting.
pip <- c(0, .1, .1, .9, 1)
truth <- c(2, 5)
t <- c(0, .1, .5, .9, 1)
counts <- e$pip_counts(pip, truth, t)
expected_tp <- sapply(c(t, Inf), function(x) sum(pip[truth] >= x))
expected_fp <- sapply(c(t, Inf), function(x) sum(pip[-truth] >= x))
stopifnot(identical(as.numeric(counts[, 'tp']), as.numeric(expected_tp)),
          identical(as.numeric(counts[, 'fp']), as.numeric(expected_fp)))

# Empty CSs do not add a fake zero to the coverage/purity denominator.
empty <- e$cs_summary(list(cs = NULL, purity = NULL), list(), c(1, 2))
stopifnot(empty['n_cs'] == 0, empty['recovered'] == 0, empty['n_causal'] == 2)
known <- e$cs_summary(list(cs = list(c(1, 2), 3),
                           purity = data.frame(min.abs.corr = c(.7, 1))),
                      list(c(1, 2), 3), c(1, 4))
stopifnot(known['n_cs'] == 2, known['covered_cs'] == 1,
          known['purity_sum'] == 1.7, known['recovered'] == 1)

# Check one real replicate against the workhorse's already-saved metrics.
saved <- new.env()
load(file.path(e$chunk_dir, selected[2]), saved)
x <- saved$results[[1]]
d <- e$replicates[e$replicates$configuration == 'add1_rec0_dom0' &
                  e$replicates$pve == .2 & e$replicates$seed == x$seed, ]
d <- d[match(x$metrics$method, d$method), ]
stopifnot(isTRUE(all.equal(d$covered_cs / d$n_cs, x$metrics$cs_coverage)),
          isTRUE(all.equal(d$recovered / d$n_causal, x$metrics$causal_recall)))
cat('\nAll limited plotting, duplicate, CS-summary and tied-PIP checks passed.\n')
