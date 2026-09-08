# Illustrate gene-tissue pairs for which the saved unweighted SuSiE-mix
# (fit_mix / susie_mix) has exactly one credible set and its lead predictor
# uses dominant coding.

library(data.table)

project_dir <- "/project2/mstephens/wdenault/susie_mix"

source(
  file.path(
    project_dir,
    "script",
    "analysis",
    "fit_mix_expression_plot_utils.R"
  )
)

load(file.path(project_dir, "res_summary.RData"))

if (!exists("res_summary")) {
  stop("res_summary.RData does not contain res_summary.")
}

association_threshold <- 5e-8
minimum_mean_reads <- 100

required_columns <- c(
  "gene",
  "tissue",
  "min_pv",
  "mean_count",
  "ncs_susie_mix",
  "n_add",
  "n_rec",
  "n_dom"
)

missing_columns <- setdiff(required_columns, names(res_summary))

if (length(missing_columns) > 0L) {
  stop(
    "res_summary is missing: ",
    paste(missing_columns, collapse = ", "),
    ". Rerun generate_summary_results.R."
  )
}

res <- as.data.table(res_summary)

# The selection concerns fit_mix itself. Requiring ncs_susie_mix == 1 and
# n_dom == 1 means the sole mixed-model CS has a dominant lead predictor.
dominant_cases <- res[
  is.finite(min_pv) &
    min_pv < association_threshold &
    is.finite(mean_count) &
    mean_count >= minimum_mean_reads &
    ncs_susie_mix == 1L &
    n_dom == 1L,
  .(
    gene,
    tissue,
    min_pv,
    mean_count,
    ncs_susie_mix,
    n_add,
    n_rec,
    n_dom
  )
]

setorder(dominant_cases, min_pv, gene, tissue)

cat(
  "One-CS dominant fit_mix candidates:",
  nrow(dominant_cases),
  "\n"
)

print(dominant_cases)

run_one_cs_fit_mix_plots(
  cases = dominant_cases,
  expected_coding = "dominant",
  plot_dir = file.path(
    project_dir,
    "plot",
    "one_cs_dominant"
  ),
  summary_filename = "fit_mix_one_cs_dominant_plot_summary.csv",
  project_dir = project_dir,
  genotype_axis_ticks = FALSE
)
