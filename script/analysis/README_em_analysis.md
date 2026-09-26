# Analyze the last completed EM iteration

Run from the project root, in this order:

```bash
Rscript --vanilla script/analysis/generate_summary_results_em.R
Rscript --vanilla script/analysis/descriptive_results_em.R
Rscript --vanilla script/analysis/finding_interesting_dominant_1cs_em.R
Rscript --vanilla script/analysis/finding_interesting_recessive_1cs_em.R
Rscript --vanilla script/analysis/plot_additive_1cs_vs_mix_2cs_em.R
```

Each script also supports RStudio **Source**: edit the settings list at its top.
The project defaults to the current project root, then the RCC path, with
`SUSIE_MIX_PROJECT_DIR` as an override. To pin all steps to iteration 12:

```bash
Rscript --vanilla script/analysis/generate_summary_results_em.R /project2/mstephens/wdenault/susie_mix 12
Rscript --vanilla script/analysis/descriptive_results_em.R /project2/mstephens/wdenault/susie_mix 12
Rscript --vanilla script/analysis/finding_interesting_dominant_1cs_em.R /project2/mstephens/wdenault/susie_mix 12
Rscript --vanilla script/analysis/finding_interesting_recessive_1cs_em.R /project2/mstephens/wdenault/susie_mix 12
Rscript --vanilla script/analysis/plot_additive_1cs_vs_mix_2cs_em.R /project2/mstephens/wdenault/susie_mix 12
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
- `plot/one_cs_dominant/` and `plot/one_cs_recessive/`: `lead_snp_comparisons.csv`, `lead_distance_summary.csv`,
  `largest_lead_shifts.csv`, and `lead_snp_distance_overview.png`.
- By default, both coding folders also contain `expression/` with individual
  `GENE_TISSUE_add_vs_mix.png` figures: additive/EM SNP PIPs above expression
  grouped by each model's lead-SNP genotype. `expression/plot_summary.csv`
  records successful plots and errors. These are the same four-panel plots as
  in the original unweighted scripts, using the selected EM fit.
  When the two lead SNPs differ, the subtitle includes their physical GRCh38
  distance (bp below 1 kb, otherwise kb). The same text is saved in the
  expression audit's `lead_comparison` column. Non-comparable coordinates are
  labeled as distance unavailable. This subtitle helper is shared with the
  original one-CS plots.
- `plot/additive_1cs_vs_mix_2cs/`: individual `GENE_TISSUE_add1_mix2.png`
  figures, `selected_cases.csv`, and `additive_1cs_vs_mix_2cs_plot_summary.csv`.
  Each figure shows additive PIPs, EM PIPs, and expression at each of the two
  EM CS leads adjusted for the other lead's predictor coding. The title reports
  the closest additive-to-EM lead distance and the ELBO + KL comparison; the
  audit retains all three pairwise lead distances, codings, PIPs, genotype
  counts, and the EM iteration. Errors are recorded as rows in the audit.

Individual plots are enabled by default (`expression_plots = TRUE`). For a
summary-only run without GTEx data, set this to `FALSE` or pass `--summary-only`.
The older `--expression-plots` flag remains accepted but is no longer needed.

The descriptive analysis keeps the existing P < 1e-8 and mean reads >= 100
thresholds. The one-CS analysis keeps the original P < 5e-8 and mean reads >= 100
thresholds and requires one **EM** CS with a dominant or recessive lead, respectively. It does not require
one additive CS by default. For multiple additive CSs it selects the highest-PIP
CS lead; with no additive CS it uses the highest-PIP SNP and labels that fallback.
Set `both_one_cs = TRUE` for strictly one CS in both models. Both overview plots
show the 20 largest physical shifts beside a distance-versus-PIP scatter. The
six largest shifts with usable PIPs are labeled using the current EM results.

`plot_additive_1cs_vs_mix_2cs_em.R` keeps the original comparison's P < 1e-8
and mean reads >= 100 thresholds. It selects **one original additive CS and
two final EM CSs**, across all coding combinations, using
`ncs_weighted_fit_mix` and `overlap_snp_add_weighted` from the selected EM
summary. The original `ncs_susie_mix` count is not used for EM selection.
Run `generate_summary_results_em.R` for the same iteration first. Individual
four-panel figures are always enabled and require the GTEx inputs. Edit
`datadir` in the entry script for a different data location. The original and
EM scripts share `add1_mix2_plot_utils.R` and `fit_mix_expression_plot_utils.R`.
The EM entry also uses `em_add1_mix2_utils.R` and the existing EM summary helpers;
sync the updated `script/analysis/` directory when updating the cluster.

### TSS plots: only disagreeing CS leads

The original `descriptive results.R`, the weighted analysis, and the EM wrapper
now exclude individual credible sets whose **biological lead SNP** is also a
CS lead in the other model for the same gene/tissue. Coding is ignored for this
decision. A shared lead is excluded even if another CS in the same gene/tissue
disagrees. Different leads remain eligible even when their CSs overlap; CSs
reported by only one model remain eligible too. Other descriptive/coding
statistics still use the full primary analysis set.

Overall and tissue-specific TSS distributions use only those selected CSs.
Each model is normalized by its number of selected, finite distances within
the existing plotted window (10-kb bins centered from -200 to +200 kb).
`tss_cs_agreement_audit.csv` records inclusion for each CS, and
`tss_cs_selection_summary.csv` counts shared leads excluded, selected CSs,
missing distances, and distances outside the window. Weighted/EM outputs use
the same filenames prefixed with `weighted_`. TSS plot filenames are unchanged;
rerun the descriptive script to replace the old all-CS distributions.

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
Rscript --vanilla script/analysis/tests/test_tss_disagreement.R
Rscript --vanilla script/analysis/tests/test_add1_mix2_comparison.R
Rscript --vanilla script/analysis/tests/test_em_add1_mix2_plots.R
```
