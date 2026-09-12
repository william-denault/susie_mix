options(susie_mix.lead_pip_plots.run = FALSE)
source("script/sim/plot_lead_pip_coding.R")
config <- readRDS(file.path(lead_pip_settings$output_dir, "plot_settings.rds"))
a <- readRDS(file.path(config$output_dir, "lead_pip_analysis.rds"))
if ("--redraw" %in% commandArgs(TRUE)) lp_render_plots(a, config)
r <- a$replicates
stopifnot(nrow(r) == 87988, sum(a$audit$included) == nrow(r),
  sum(a$audit$duplicates) == 9, sum(a$audit$errors) == 12, sum(a$audit$invalid) == 0,
  all(r$outcome %in% lp_outcomes), all(r$lead_pip >= 0 & r$lead_pip <= 1),
  all(is.na(r$called_coding) == (r$outcome %in% c("ambiguous", "no_support"))),
  all(r$n_top[r$outcome == "ambiguous"] > 1), all(r$n_top[!is.na(r$called_coding)] == 1),
  all(r$n_top_exact[r$outcome == "exact"] == 1),
  all(r$n_top_exact[r$outcome %in% c("wrong_coding", "wrong_snp")] == 0),
  all(r$called_coding[r$scenario == "Additive only" & r$outcome == "exact"] == "additive"))
for (name in c("bins", "summary", "conditions", "additive")) {
  x <- a[[name]]
  stopifnot(all(rowSums(x[paste0("n_", lp_outcomes)]) == x$n_runs),
    all(rowSums(x[paste0("n_lead_", cm_classes)]) == x$n_unique),
    all(abs(rowSums(x[paste0("fraction_", lp_outcomes)]) - 1) < 1e-12))
}
for (scope in unique(a$summary$scope)) {
  raw <- if (scope == "Additive only") r[r$scenario == "Additive only", ] else r
  bins <- a$bins[a$bins$scope == scope, ]
  acc <- a$accuracy[a$accuracy$scope == scope, ]
  stopifnot(sum(bins$n_runs) == nrow(raw), sum(acc$n_leads) == sum(!is.na(raw$called_coding)),
    sum(acc$n_exact) == sum(raw$outcome == "exact"))
}
stopifnot(length(list.files(config$output_dir, "[.]pdf$")) == 4,
  length(list.files(config$output_dir, "[.]png$")) == 4)
print(table(r$outcome))
print(a$additive[c("pve", "K", "n_runs", "n_unique", "n_ambiguous", "share_unique_additive", "share_unique_recessive", "share_unique_dominant")], row.names = FALSE)
cat("All full-data lead-PIP consistency checks passed.\n")
