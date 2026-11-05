# Analyze the last completed EM iteration

For the separate first-stage GWAS colocalization screen, Source
`screen_em_coloc_candidates.R` or run it with Rscript. It uses `latest` completed
fits, P < 5e-8, shifts >10 kb from every vanilla CS lead, LD/CS-overlap diagnostics,
and SNP-level log-BF exports. It does not require pre-generated summaries.
See [the screening guide](README_em_coloc_screen.md) and [trait shortlist](coloc_trait_panel.csv).

For the next stage, see [GWAS downloads and paired EM/vanilla colocalization](README_em_gwas_coloc.md).
Source `download_em_coloc_gwas.R`, then `run_em_gwas_coloc.R` with the printed frozen job path.
New jobs include 50 fixed GWAS/biomarker traits plus gene-matched UKB-PPP proteins
and external eQTL Catalogue signals. The screen can also export unshifted controls;
see the coloc guide for Synapse access and the prespecified tissue matches.

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
thresholds. The interesting one-CS plots use only P < 5e-8 for eligibility,
with no sample-size or mean-read-count cutoff, and require one **EM** CS with
a dominant or recessive lead, respectively. The original unweighted one-CS
plots use the same P-value eligibility rule. Expression reconstruction still checks
that the data match the saved fit. Both EM one-CS scripts default to
`both_one_cs = TRUE`: exactly one additive CS and exactly one EM CS per
gene/tissue pair. This does not require a CS to contain only one SNP.
Set `both_one_cs = FALSE` only to include other additive CS counts; in that
mode, multiple additive CSs use their highest-PIP CS lead, and no additive CS
uses the highest-PIP SNP with an explicit fallback label. Both overview plots
show the 20 largest physical shifts beside a distance-versus-PIP scatter. The
six largest shifts with usable PIPs are labeled using the current EM results.
The overview subtitle counts shifts strictly greater than 10 kb and 100 kb.
These counts also appear as `n_over_10kb` and `n_over_100kb` in
`lead_distance_summary.csv`, using all comparable lead pairs in the selection.
The overview is always written before optional expression plots, including
with `--summary-only`, to `plot/one_cs_dominant/lead_snp_distance_overview.png`
or `plot/one_cs_recessive/lead_snp_distance_overview.png` within the selected
EM iteration. Sync the entry scripts, `em_one_cs_utils.R`, and
`one_cs_lead_distance_plot_utils.R` together when updating the cluster.

`plot_additive_1cs_vs_mix_2cs_em.R` uses P < 5e-8 with no mean-read-count or
sample-size cutoff, matching the one-CS EM scripts. It selects **one original additive CS and
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
The original, weighted/EM, and slide descriptive scripts now plot Gaussian
kernel density estimates of the raw signed distances, with one common
`tss_bandwidth_kb = 10` for both models and all tissues. Change that setting
to 5 or 20 to check sensitivity. The y-axis is density per kb.

Each KDE uses **all finite selected distances**, including leads outside the
displayed -200 to +200 kb window. Each model's density integrates to one over
the full domain; the visible portion need not integrate to one. Legend `n`
therefore includes all finite selected leads and can exceed the old in-window
count. No observation is moved to zero or discarded to smooth the curve.

The existing `*tss_distance_distribution.csv` filenames now contain
`distance_to_tss_kb`, `density`, `bandwidth_kb`, `kernel`, `normalization`,
`n_cs_selected`, `n_cs_total` (finite distances used), `n_cs_in_window`,
`n_cs_missing_distance`, and `n_cs_outside_window`. The old histogram
`count`, `proportion`, `percentage`, and bin-boundary columns are replaced.
Empty groups have NA densities and zero usable observations. A single usable
lead is displayed with the same fixed Gaussian bandwidth, with n = 1.

`tss_cs_agreement_audit.csv` records inclusion for each CS, and
`tss_cs_selection_summary.csv` counts shared leads excluded, selected CSs,
missing distances, and distances outside the exact +/-200 kb display window
(there are no outer half bins). Weighted/EM outputs use
the same filenames prefixed with `weighted_`. TSS plot filenames are unchanged;
rerun the descriptive script to replace the binned curves. Copy the updated
`tss_disagreement_utils.R` along with the descriptive scripts to the cluster.

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

## Four-panel CS agreement figure

After generating the EM summaries for the converged iteration, Source
`plot_cs_agreement_em.R` in RStudio, or run from the project root:

```bash
Rscript --vanilla script/analysis/plot_cs_agreement_em.R /project2/mstephens/wdenault/susie_mix 18
```

Use `latest` instead of `18` to select the latest complete iteration using the
same checks as the other EM scripts. This does not itself test convergence of
the EM prior updates; pin the iteration you have identified as converged.
The script uses `weighted_fit_mix` from that iteration, never SuSiE-slide or
the original fixed-prior fit. It needs the region and CS summaries produced by
`generate_summary_results_em.R`; no expression data or refitting is required.

Outputs go to `results_em/iteration_018/descriptive_results/cs_agreement/`
(with the selected iteration number). The main figure is
`em_cs_agreement_4panel.png`, with a vector PDF of the same name:

- **A: four tissue CS-count tables.** Rows are additive SuSiE CS counts,
  columns are EM CS counts, and cells are numbers of gene/tissue regions.
  Green diagonal cells have equal counts; orange off-diagonal cells differ.
  Select the four tissues with the lowest equal-count rate (highest
  disagreement). Ties use larger sample size, then tissue name.
  `minimum_tissue_regions = 100` avoids ranking very small cohorts; edit the
  setting as needed. With fewer than four qualifying tissues, the script warns
  and displays the available tissues without duplication.
- **B: lead-distance empirical cumulative distributions across all tissues.**
  One curve describes non-identical paired leads, including distinct variants
  at the same position. Two separate curves describe extra EM or extra SuSiE
  CSs by their distance to the nearest lead in the other model. The transform
  is `log10(1 + distance_bp)`, with physical-distance tick labels and zero
  retained. Identical paired leads are counted separately. Extra CSs in
  regions with no other-model CS have no distance and are counted separately.
- **C: extra-lead TSS-distance empirical cumulative distributions.** Use the
  exact extra EM and extra SuSiE CS groups in B, with matching colors and line
  styles. Each distance is from that extra CS's lead SNP to its own gene's TSS,
  using `abs(distance_to_tss_kb) * 1000` and the same logarithmic spacing as B.
  Zero distances are retained. Model-only regions remain excluded in both
  panels; missing TSS distances are counted separately. Legend counts show
  the number with TSS information over the corresponding group size in B.
- **D: recessive and dominant EM coding priors by tissue.** Separate points
  show `pi_rec` and `pi_dom` from the selected iteration's validated
  `priors.csv`; tissues are sorted by `pi_rec + pi_dom`. These probabilities
  retain additive prior mass in the denominator and are not renormalized
  within nonadditive codings. They summarize the full EM fitting cohort,
  independently of the A-C filters, and are not lead PIPs or significance
  tests. The caption gives the source iteration when available (e.g., priors
  used in iteration 18 were estimated from iteration 17 fits). All prior-table
  tissues are shown, including any without regions passing the A-C filters.

The default cohort for A-C matches `descriptive_results_em.R`: P < 1e-8 and mean reads
>= 100, then at least one CS in either model. Zero/zero regions are excluded
from tissue ranking and retained in `region_audit.csv`. One-sided zero-CS
regions count as disagreement. CS size is unrestricted. Each tissue's rate
uses regions, each distance curve uses CS pairs or unpaired CSs; multi-CS
regions therefore contribute several distances. This is a descriptive pooled
distribution, without independence assumptions or inferential intervals.

Matching reserves identical biological lead SNP IDs first, ignoring coding,
and pairs each CS at most once. Remaining leads use a one-to-one assignment
minimizing total absolute distance, with unmatched CSs allowed on the larger
side. For 1 SuSiE CS and 2 EM CSs with one shared lead, the shared lead is paired
at zero distance, and the second EM CS is extra. Its nearest-lead distance is
a reference distance, not a second match to the same SuSiE CS. Matching is
deterministic under input reordering. It is positional and does not establish
shared signals, causality, CS equivalence, or low LD.

Incomplete CS summaries and invalid coordinates are audited. Valid count
comparisons remain in panel A; valid lead identities remain in the exported
retention summary even when coordinates cannot be used in B/C. Invalid regions are never partially
matched or treated as evidence of extra CSs. Both summaries must have a single
matching EM iteration, and the entry script checks saved summary provenance.

The script prints the selected contingency tables using R's table display and
leaves `em_agreement_results` in the RStudio session. It also saves:
`tissue_ranking.csv`, `selected_tissues.csv`, `tissue_cs_count_tables.csv`,
`matched_cs_pairs.csv`, `unpaired_cs.csv`, `lead_distance_summary.csv`,
`tss_distance_summary.csv`, `tissue_coding_priors.csv`,
`lead_retention_summary.csv`, `extra_em_coding_summary.csv`,
`region_audit.csv`, `orphan_cs_rows.csv`, and `analysis_settings.csv`.
`unpaired_cs.csv` includes signed and absolute TSS distances and the panel-C
inclusion/exclusion status. The prior export retains estimation metadata and
adds the number of regions passing the A-C filters. The former panels C/D
remain available as retention and extra-CS coding CSV summaries.
Sync `plot_cs_agreement_em.R` and `em_cs_agreement_utils.R` along with the
existing `em_analysis_utils.R` when updating the cluster.

## Reusable functions and tests

The shared `generate_summary_results.R` and `descriptive_results_weighted.R`
now define reusable functions when sourced. Their Rscript entry points still
run the original analyses. To use them interactively outside the EM wrappers,
call `generate_summary_results(path_res, output_dir)` or
`run_weighted_descriptive_results(summary_file, cs_summary_file, output_dir)`.

```bash
Rscript --vanilla script/analysis/tests/test_em_analysis.R
Rscript --vanilla script/analysis/tests/test_em_cs_agreement.R
Rscript --vanilla script/analysis/tests/test_weighted_summary.R
Rscript --vanilla script/analysis/tests/test_one_cs_plots.R
Rscript --vanilla script/analysis/tests/test_one_cs_lead_distance_plots.R
Rscript --vanilla script/analysis/tests/test_tss_disagreement.R
Rscript --vanilla script/analysis/tests/test_add1_mix2_comparison.R
Rscript --vanilla script/analysis/tests/test_em_add1_mix2_plots.R
```
