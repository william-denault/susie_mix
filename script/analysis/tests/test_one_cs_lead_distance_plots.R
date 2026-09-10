# Run from the repository root; only base R is required:
# Rscript --vanilla script/analysis/tests/test_one_cs_lead_distance_plots.R
source("script/analysis/one_cs_lead_distance_plot_utils.R")
output_dir <- file.path("tmp", "one_cs_distance_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use different counted alleles and alternating upstream/downstream shifts
# to check that genomic coordinates, not predictor indices or alleles, are used.
fixture <- data.frame(
  gene = sprintf("Gene%02d", 1:25), tissue = "Test tissue", status = "plotted",
  additive_lead_snp = "chr1_1000000_A_G_b38_A",
  lead_snp = paste0("chr1_", 1000000 + (1:25) * 1000 * rep(c(-1, 1), length.out = 25),
                     "_C_T_b38_C"),
  additive_n_cs = rep(c(0L, 1L, 2L, NA_integer_), length.out = 25),
  n_cs = 1L, lead_coding = "dominant", lead_pip = seq(.1, .9, length.out = 25)
)
failed <- fixture[1, ]
failed$gene <- "Failed"
failed$status <- "error"
failed$lead_snp <- "chr1_999999999_A_G_b38_A"
same <- fixture[1, ]
same$gene <- "Same SNP"
same$lead_snp <- same$additive_lead_snp
input <- rbind(fixture, failed, same)
result <- plot_one_cs_lead_distance_overview(
  input, "dominant", file.path(output_dir, "top20.png"), top_n = 20L
)
stopifnot(
  result$n_completed == 26L, result$n_changed == 25L,
  result$n_scatter == 25L, result$n_invalid == 0L,
  identical(result$top_cases$gene, sprintf("Gene%02d", 25:6)),
  identical(result$distance_data$distance_bp, as.numeric((25:1) * 1000)),
  !"distance_bp" %in% names(input),
  file.info(result$output_file)$size > 10000,
  grDevices::dev.cur() == 1L
)

# Smaller successful cohorts should render up to the available count.
small <- plot_one_cs_lead_distance_overview(
  fixture[1:2, ], "dominant", file.path(output_dir, "two_cases.png")
)
stopifnot(nrow(small$top_cases) == 2L)

# Same-position different alleles have distance zero, not a fabricated log
# coordinate. Missing PIP excludes only the scatter point, not the ranking.
zero <- fixture[1, ]
zero$gene <- "Same position"
zero$lead_snp <- "chr1_1000000_C_T_b38_C"
missing_pip <- fixture[2, ]
missing_pip$lead_pip <- NA_real_
edge <- plot_one_cs_lead_distance_overview(
  rbind(zero, missing_pip), "dominant", file.path(output_dir, "zero_missing_pip.png")
)
stopifnot(nrow(edge$top_cases) == 2L, edge$n_changed == 2L,
          edge$n_scatter == 0L, identical(edge$distance_data$distance_bp, c(2000, 0)))

# Do not compute inter-chromosome, mixed-build or malformed-ID distances.
invalid <- fixture[1:3, ]
invalid$lead_snp <- c("chr2_1001000_C_T_b38_C", "chr1_1001000_C_T_b37_C", "invalid")
warnings <- character()
bad <- withCallingHandlers(
  plot_one_cs_lead_distance_overview(invalid, "dominant", file.path(output_dir, "invalid.png")),
  warning = function(e) { warnings <<- c(warnings, conditionMessage(e)); invokeRestart("muffleWarning") }
)
stopifnot(bad$n_invalid == 3L, nrow(bad$top_cases) == 0L,
          length(warnings) == 1L, grepl("comparable GRCh38", warnings))

# Empty and all-failed runs have only audit columns in the runner's output.
for (status in list(character(), "error")) {
  empty <- plot_one_cs_lead_distance_overview(
    data.frame(status = status), "recessive", file.path(output_dir, "empty.png")
  )
  stopifnot(empty$n_completed == 0L, empty$n_changed == 0L,
            nrow(empty$top_cases) == 0L, file.info(empty$output_file)$size > 10000)
}
unchanged <- plot_one_cs_lead_distance_overview(
  same, "dominant", file.path(output_dir, "same_lead.png")
)
stopifnot(unchanged$n_completed == 1L, unchanged$n_changed == 0L)

# Execute each entry script's actual final two expressions with its expensive
# four-panel runner stubbed. The real overview must consume that run's rows
# after the runner finishes, without requiring a ranking CSV or GTEx inputs.
for (coding in c("dominant", "recessive")) {
  env <- new.env(parent = globalenv())
  env$project_dir <- file.path(output_dir, paste0("entry_", coding))
  env$runner_finished <- FALSE
  env$fixture <- fixture
  env$fixture$lead_coding <- coding
  assign(paste0(coding, "_cases"), data.frame(gene = "unused"), envir = env)
  env$run_one_cs_fit_mix_plots <- function(...) {
    env$runner_finished <- TRUE
    env$fixture
  }
  env$plot_one_cs_lead_distance_overview <- function(...) {
    stopifnot(env$runner_finished)
    env$overview <- plot_one_cs_lead_distance_overview(...)
  }
  expressions <- parse(file.path("script", "analysis", paste0("finding_interesting_", coding, "_1cs.R")))
  for (expr in tail(expressions, 2L)) eval(expr, env)
  stopifnot(nrow(env$overview$top_cases) == 20L,
            env$overview$n_completed == 25L, file.exists(env$overview$output_file))
}
stopifnot(grDevices::dev.cur() == 1L)
cat("One-CS lead-distance overview tests passed.\n")
