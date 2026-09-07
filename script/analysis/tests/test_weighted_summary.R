# Run from the repository root:
# Rscript --vanilla script/analysis/tests/test_weighted_summary.R
# Uses base R and synthetic fits; no GTEx data or susieR installation required.

summary_path <- "script/analysis/generate_summary_results.R"
program <- parse(file = summary_path)

load_functions <- function(expressions, envir) {
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
      eval(expr, envir)
    }
  }
}

summary_env <- new.env(parent = globalenv())
load_functions(program, summary_env)
utils_env <- new.env(parent = globalenv())
load_functions(parse("script/scan_tissue_attempt/workhorse_utils.R"), utils_env)

make_fit <- function(p, cs, leads = integer(0), components = seq_along(cs)) {
  alpha <- matrix(1 / p, nrow = 4L, ncol = p)
  for (i in seq_along(cs)) {
    alpha[components[i], ] <- 0.001
    alpha[components[i], cs[[i]]] <- 0.1
    alpha[components[i], leads[i]] <- 0.8
    alpha[components[i], ] <- alpha[components[i], ] / sum(alpha[components[i], ])
  }
  list(
    alpha = alpha,
    pip = 1 - apply(1 - alpha, 2L, prod),
    elbo = c(-12, -10),
    KL = rep(0.5, 4L),
    converged = TRUE,
    sets = if (length(cs) == 0L) NULL else list(
      cs = cs, cs_index = components,
      coverage = rep(0.95, length(cs)),
      purity = data.frame(min.abs.corr = rep(0.8, length(cs)),
                          mean.abs.corr = 0.9, median.abs.corr = 0.9)
    )
  )
}

snps <- c("chr1_100_A_G_b38_A", "chr1_200_C_T_b38_C", "chr1_300_G_T_b38_G")
add_map <- data.frame(predictor_index = 1:3, predictor_name = snps,
                      snp = snps, coding = "additive")
mix_coding <- rep(c("additive", "recessive", "dominant"), each = 3L)
mix_map <- data.frame(
  predictor_index = 1:9,
  predictor_name = paste(rep(snps, 3L), mix_coding, sep = "__"),
  snp = rep(snps, 3L), coding = mix_coding
)
add_fit <- make_fit(3L, list(L1 = c(1L, 2L), L3 = 3L), c(1L, 3L), c(1L, 3L))
mix_fit <- make_fit(9L, list(L3 = c(1L, 2L), L1 = 6L), c(1L, 6L), c(3L, 1L))
weighted_fit <- make_fit(9L, list(L2 = c(7L, 2L), L4 = 6L), c(7L, 6L), c(2L, 4L))
weighted_fit$elbo <- c(-10, -8)
weighted_fit$converged <- FALSE

x <- list(
  susie_add = add_fit, susie_add_perm = make_fit(3L, list()),
  susie_mix = mix_fit, susie_mix_perm = make_fit(9L, list()),
  weighted_fit_mix = weighted_fit,
  add_predictor_map = add_map, mix_predictor_map = mix_map,
  mix_coding_prior = c(dominant = 0.15, additive = 0.8, recessive = 0.05),
  weighted_mix_prior_weights = setNames(
    c(rep(0.8 / 3, 3), rep(0.05 / 3, 3), rep(0.15 / 3, 3)),
    mix_map$predictor_name
  ),
  n_ind = 100L, n_SNP = 3L, n_mix_predictor = 9L,
  n_add_rm = 0L, n_rec_rm = 0L, n_dom_rm = 0L,
  mean_read = 150, median_read = 120, min_pv = 1e-10
)
# Match the workhorse's named numeric distance vectors (no separate TSS field).
for (key in c("susie_add", "susie_mix", "weighted_fit_mix")) {
  map <- if (key == "susie_add") add_map else mix_map
  x[[paste0(key, "_lead_snp_tss_distance")]] <-
    utils_env$get_cs_lead_tss_distance(x[[key]], map, tss = 150)
}

summarize <- function(value) {
  summary_env$summarize_tissue(value, "GENE", "TISSUE", "GENE.rds")
}
cs_summary <- function(value, row = summarize(value)) {
  summary_env$summarize_credible_sets(value, "GENE", "TISSUE", "GENE.rds", row)
}
expect_error <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}

row <- summarize(x)
stopifnot(
  row$has_weighted_fit_mix, row$ncs_weighted_fit_mix == 2L,
  row$elbo_weighted_mix == -8, row$dif_elbo_weighted_vs_add == 2,
  row$dif_elbo_weighted_vs_mix == 2, row$log_lik_weighted_mix == -6,
  identical(row$converged_weighted_fit_mix, FALSE),
  row$n_add_weighted == 0L, row$n_dom_weighted == 1L, row$n_rec_weighted == 1L,
  row$mix_weighted_overlap_snp == 3L, row$mix_weighted_jaccard_snp == 1,
  row$mix_weighted_same_cs_count, row$mix_weighted_same_cs_snp_sets,
  !row$mix_weighted_same_cs_predictor_sets,
  row$mix_weighted_same_lead_snp_set, !row$mix_weighted_same_lead_predictor_set,
  abs(row$weighted_prior_realized_additive - 0.8) < 1e-12,
  abs(row$weighted_prior_realized_dominant - 0.15) < 1e-12,
  abs(row$weighted_prior_realized_recessive - 0.05) < 1e-12
)

cs <- cs_summary(x, row)
wcs <- cs[cs$model_key == "weighted_fit_mix", ]
stopifnot(
  nrow(cs) == 6L, nrow(wcs) == row$ncs_weighted_fit_mix,
  all(wcs$model == "Weighted SuSiE-mix"),
  identical(wcs$lead_snp, snps[c(1, 3)]),
  identical(wcs$lead_coding, c("dominant", "recessive")),
  identical(wcs$component, c(2, 4)),
  identical(wcs$distance_to_tss_bp, c(-50, 150)),
  all(wcs$ncs_weighted_fit_mix == 2L), !anyDuplicated(cs$cs_id)
)

# Keep the generator and dedicated weighted-analysis contracts synchronized.
weighted_descriptive_program <- parse(
  file = "script/analysis/descriptive_results_weighted.R"
)
read_assignment <- function(program, name) {
  match <- vapply(program, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && identical(as.character(expr[[2]]), name)
  }, logical(1L))
  stopifnot(sum(match) == 1L)
  eval(program[[which(match)]], new.env(parent = globalenv()))
}
required_summary_columns <- read_assignment(
  weighted_descriptive_program,
  "required_columns"
)
required_cs_columns <- read_assignment(
  weighted_descriptive_program,
  "required_cs_columns"
)
stopifnot(length(setdiff(required_summary_columns, names(row))) == 0L,
          length(setdiff(required_cs_columns, names(cs))) == 0L)
# Prove that weighted distances are read from the weighted saved field.
different_distance <- x
different_distance$weighted_fit_mix_lead_snp_tss_distance[] <- c(123, 456)
changed_cs <- cs_summary(different_distance)
stopifnot(identical(
  changed_cs$workhorse_distance_to_tss_bp[changed_cs$model_key == "weighted_fit_mix"],
  c(123, 456)
))

# Legacy metadata must not partially match the missing weighted fit name.
legacy <- x
legacy$weighted_fit_mix <- NULL
legacy_row <- summarize(legacy)
new_columns <- names(summary_env$summarize_weighted_mix(legacy))
old_columns <- setdiff(names(row), new_columns)
stopifnot(
  !legacy_row$has_weighted_fit_mix,
  all(is.na(legacy_row[setdiff(new_columns, "has_weighted_fit_mix")])),
  identical(row[old_columns], legacy_row[old_columns]),
  !anyDuplicated(names(row)),
  !any(grepl("weighted.*perm|perm.*weighted", names(row)))
)
legacy_cs <- cs_summary(legacy)
original_cs <- cs[cs$model_key != "weighted_fit_mix", setdiff(names(cs), new_columns)]
stopifnot(identical(original_cs, legacy_cs[setdiff(names(legacy_cs), new_columns)]))

# Optionally compare every pre-existing value with an earlier script supplied
# as the first command-line argument. The default test is self-contained.
baseline_paths <- commandArgs(trailingOnly = TRUE)
if (length(baseline_paths) > 0L) {
  baseline <- new.env(parent = globalenv())
  load_functions(parse(file = baseline_paths[1]), baseline)
  baseline_row <- baseline$summarize_tissue(x, "GENE", "TISSUE", "GENE.rds")
  baseline_cs <- baseline$summarize_credible_sets(
    x, "GENE", "TISSUE", "GENE.rds", baseline_row
  )
  stopifnot(identical(row[names(baseline_row)], baseline_row),
            identical(original_cs[names(baseline_cs)], baseline_cs))
  cat("PASS: all original summary values match the supplied baseline script.\n")
}

zero <- x
zero$weighted_fit_mix <- make_fit(9L, list())
zero$weighted_fit_mix_lead_snp_tss_distance <- setNames(numeric(0), character(0))
zero_row <- summarize(zero)
stopifnot(
  zero_row$has_weighted_fit_mix, zero_row$ncs_weighted_fit_mix == 0L,
  zero_row$n_add_weighted == 0L, zero_row$mix_weighted_overlap_snp == 0,
  is.na(zero_row$mix_weighted_same_lead_snp_set),
  !any(cs_summary(zero)$model_key == "weighted_fit_mix"),
  nrow(rbind(row, legacy_row, zero_row)) == 3L
)
only_weighted <- x
only_weighted$susie_add <- make_fit(3L, list())
only_weighted$susie_mix <- make_fit(9L, list())
stopifnot(all(cs_summary(only_weighted)$model_key == "weighted_fit_mix"))

# Renaming/reordering CS components must not affect set agreement.
reordered <- mix_fit
reordered$sets$cs <- setNames(rev(mix_fit$sets$cs), c("renamed1", "renamed2"))
reordered$sets$cs_index <- rev(mix_fit$sets$cs_index)
comparison <- summary_env$compare_mixed_fits(mix_fit, reordered, mix_map)
stopifnot(comparison$mix_weighted_same_cs_snp_sets,
          comparison$mix_weighted_same_cs_predictor_sets)
partial <- summary_env$compare_mixed_fits(
  make_fit(9L, list(A = c(1L, 2L)), 1L),
  make_fit(9L, list(B = c(2L, 3L)), 2L), mix_map
)
stopifnot(abs(partial$mix_weighted_jaccard_snp - 1 / 3) < 1e-12,
          partial$mix_weighted_overlap_snp == 1L)

shuffled <- x
shuffled$weighted_mix_prior_weights <- rev(x$weighted_mix_prior_weights)
stopifnot(identical(summarize(shuffled), row))
# If a whole coding class was filtered out, report its realized mass as zero.
no_recessive <- x
no_recessive$mix_predictor_map <- mix_map[c(1:3, 7:9), ]
no_recessive$mix_predictor_map$predictor_index <- 1:6
no_recessive$susie_mix <- no_recessive$weighted_fit_mix <-
  make_fit(6L, list(L1 = c(1L, 2L)), 1L)
no_recessive$susie_mix_perm <- make_fit(6L, list())
no_recessive$weighted_mix_prior_weights <- c(rep(0.8 / 0.95 / 3, 3),
                                           rep(0.15 / 0.95 / 3, 3))
no_rec_row <- summarize(no_recessive)
stopifnot(no_rec_row$weighted_prior_realized_recessive == 0,
          no_rec_row$weighted_prior_target_recessive == 0.05,
          abs(no_rec_row$weighted_prior_realized_additive - 0.8 / 0.95) < 1e-12)
no_metadata <- x
no_metadata$mix_coding_prior <- no_metadata$weighted_mix_prior_weights <- NULL
stopifnot(summarize(no_metadata)$has_weighted_fit_mix,
          is.na(summarize(no_metadata)$weighted_prior_realized_additive))
bad_weights <- x
names(bad_weights$weighted_mix_prior_weights)[1] <- "wrong_predictor"
expect_error(summarize(bad_weights), "names do not match")
bad_fit <- x
bad_fit$weighted_fit_mix$alpha <- matrix(1, 4, 2)
expect_error(summarize(bad_fit), "weighted_fit_mix_predictor_map")

# Execute the complete summarizer, overriding only its path assignments.
# All RDS/CSV/RData fixtures live in this R session's temporary directory.
run_pipeline <- function(fixtures) {
  root <- tempfile("weighted-summary-test-")
  dir.create(root)
  input <- file.path(root, "input")
  dir.create(input)
  for (gene in names(fixtures)) {
    saveRDS(list(TISSUE = fixtures[[gene]]), file.path(input, paste0(gene, ".rds")))
  }
  outputs <- c("summary_file", "summary_csv", "cs_summary_file", "cs_summary_csv",
               "cs_error_csv", "error_csv", "failed_gene_file")
  extensions <- c("RData", "csv", "RData", "csv", "csv", "csv", "txt")
  paths <- c(list(path_res = input), setNames(as.list(
    file.path(root, paste0(outputs, ".", extensions))
  ), outputs))
  run_env <- new.env(parent = globalenv())
  warnings <- character(0)
  log <- capture.output(withCallingHandlers({
    for (expr in program) {
      name <- if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
                  is.symbol(expr[[2]])) as.character(expr[[2]]) else ""
      if (name %in% names(paths)) {
        assign(name, paths[[name]], envir = run_env)
      } else {
        eval(expr, envir = run_env)
      }
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }))
  stopifnot(all(file.exists(unlist(paths[outputs]))))
  saved <- new.env()
  load(paths$summary_file, saved)
  load(paths$cs_summary_file, saved)
  stopifnot(identical(saved$res_summary, run_env$res_summary),
            identical(saved$res_cs_summary, run_env$res_cs_summary))
  csv <- read.csv(paths$summary_csv, check.names = FALSE)
  # The existing pipeline writes a zero-column CSV when no model has a CS.
  cs_csv <- if (ncol(saved$res_cs_summary) > 0L) {
    read.csv(paths$cs_summary_csv, check.names = FALSE)
  } else {
    data.frame()
  }
  stopifnot(identical(names(csv), names(saved$res_summary)))
  list(saved = saved, csv = csv, cs_csv = cs_csv, warnings = warnings, log = log)
}

pipeline <- run_pipeline(list(NEW = x, LEGACY = legacy, ZERO = zero, BAD = bad_fit))
stopifnot(
  nrow(pipeline$csv) == 3L, nrow(pipeline$saved$res_errors) == 1L,
  pipeline$saved$res_errors$gene == "BAD", nrow(pipeline$saved$cs_errors) == 0L,
  length(pipeline$warnings) == 0L,
  sum(pipeline$csv$has_weighted_fit_mix) == 2L,
  sum(pipeline$cs_csv$model_key == "weighted_fit_mix") == 2L,
  pipeline$csv$ncs_weighted_fit_mix[pipeline$csv$gene == "ZERO"] == 0L,
  is.na(pipeline$csv$ncs_weighted_fit_mix[pipeline$csv$gene == "LEGACY"])
)
legacy_pipeline <- run_pipeline(list(OLD = legacy))
stopifnot(any(grepl("No weighted_fit_mix was found", legacy_pipeline$warnings)))
all_zero <- zero
all_zero$susie_add <- make_fit(3L, list())
all_zero$susie_mix <- make_fit(9L, list())
empty_cs_pipeline <- run_pipeline(list(NO_CS = all_zero))
stopifnot(nrow(empty_cs_pipeline$saved$res_cs_summary) == 0L,
          empty_cs_pipeline$csv$has_weighted_fit_mix,
          empty_cs_pipeline$csv$ncs_weighted_fit_mix == 0L,
          length(empty_cs_pipeline$warnings) == 0L)

cat("PASS: weighted metrics, coding-aware agreement, component order, saved TSS distances,\n",
    "legacy/zero-CS cases, prior alignment, unchanged original summaries,\n",
    "and full RDS-to-CSV/RData pipeline including invalid-fit error reporting.\n")
