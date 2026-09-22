# Analyze the last completed EM iteration

Run from the project root, in this order:

```bash
Rscript --vanilla script/analysis/generate_summary_results_em.R
Rscript --vanilla script/analysis/descriptive_results_em.R
Rscript --vanilla script/analysis/finding_interesting_dominant_1cs_em.R
```

Each script also supports RStudio **Source**: edit the settings list at its top.
The project defaults to the current project root, then the RCC path, with
`SUSIE_MIX_PROJECT_DIR` as an override. To pin all steps to iteration 12:

```bash
Rscript --vanilla script/analysis/generate_summary_results_em.R /project2/mstephens/wdenault/susie_mix 12
Rscript --vanilla script/analysis/descriptive_results_em.R /project2/mstephens/wdenault/susie_mix 12
Rscript --vanilla script/analysis/finding_interesting_dominant_1cs_em.R /project2/mstephens/wdenault/susie_mix 12 --expression-plots
```

`latest` selects the highest iteration in `prior_history.csv` with a valid
manifest, every chunk completion marker, all gene RDS files, and matching saved
priors. Prepared or partially copied iterations are reported and skipped.
An explicitly requested iteration must be complete. Completion includes saved
gene/tissue error records; it does not imply every fit succeeded. The number 12
is not hard-coded, and no EM fitting or additional prior update is performed.

## Inputs and interpretation

The EM workhorse saves only `weighted_fit_mix`. The generator joins each EM
gene/tissue with `results/GENE.rds` for additive and original unweighted fits,
permutations, marginal association P-values, and metadata. It checks predictor
order, sample/phenotype metadata where available, convergence, iteration tags,
and coding priors. The original stored weighted fit is replaced by the chosen
EM fit; it is never used as a fallback for missing EM results. Missing tissues,
failed genes, and invalid comparisons are recorded in the error CSVs.

The summary retains the established object names `res_summary` and
`res_cs_summary`, adds `em_iteration`, and uses the existing `*_weighted*`
columns and `model_key == "weighted_fit_mix"` for the final EM model. The plain
`susie_mix` columns still describe the original unweighted model. Original
permutation fits are retained as baseline information: there is no EM
permutation fit. Use `descriptive_results_em.R` for the final EM analysis;
pointing the original unweighted `descriptive results.R` at this summary would
still analyze the original unweighted columns.

## Outputs

All new products live under `results_em/iteration_012/` (or the selected iteration):

- `summary/`: `res_summary.RData`, `res_cs_summary.RData`, CSV versions,
  `res_errors.csv`, `res_cs_errors.csv`, `failed_genes.txt`,
  `em_summary_metadata.rds`, and `em_priors_used.csv`.
- `descriptive_results/`: the existing weighted descriptive CSV/PDF outputs,
  now describing the final EM fit versus additive and original unweighted fits.
- `plot/one_cs_dominant/`: `lead_snp_comparisons.csv`, `lead_distance_summary.csv`,
  `largest_lead_shifts.csv`, and `lead_snp_distance_overview.png`.
- With `expression_plots = TRUE` or `--expression-plots`, the last folder also
  contains `expression/` with four-panel PIP/expression plots and an error audit.

The descriptive analysis keeps the existing P < 1e-8 and mean reads >= 100
thresholds. The one-CS analysis keeps the original P < 5e-8 and mean reads >= 100
thresholds and requires one **EM** CS with a dominant lead. It does not require
one additive CS by default. For multiple additive CSs it selects the highest-PIP
CS lead; with no additive CS it uses the highest-PIP SNP and labels that fallback.
Set `both_one_cs = TRUE` for strictly one CS in both models. Set
`expected_coding = "recessive"` to run the same analysis for recessive leads.

Distances use biological SNP coordinates, retain unchanged leads and zero
distances in the comparison table, and audit unavailable coordinates. The
overview ranks changed leads and distinguishes additive CS counts. Overlap and
Jaccard columns compare the selected additive CS with the sole EM CS. A large
physical shift does not establish low LD. The ELBO + KL difference is a
descriptive metric, not a calibrated likelihood-ratio test. Distances and the
overview use all eligible saved fits, independent of expression-plot success.

Summary generation and distance plots use base R. Descriptive tables require
`data.table`; expression plots also require `matrixStats`, `susieR`, PLINK, and
the original GTEx data at `datadir`. The local checkout currently contains EM
metadata but no EM gene result directories, so actual iteration-12 outputs must
be generated on the cluster or after copying those RDS files locally.

## Reusable functions and tests

The shared `generate_summary_results.R` and `descriptive_results_weighted.R`
now define reusable functions when sourced. Their Rscript entry points still
run the original analyses. To use them interactively outside the EM wrappers,
call `generate_summary_results(path_res, output_dir)` or
`run_weighted_descriptive_results(summary_file, cs_summary_file, output_dir)`.

```bash
Rscript --vanilla script/analysis/tests/test_em_analysis.R
Rscript --vanilla script/analysis/tests/test_weighted_summary.R
Rscript --vanilla script/analysis/tests/test_one_cs_plots.R
Rscript --vanilla script/analysis/tests/test_one_cs_lead_distance_plots.R
```
