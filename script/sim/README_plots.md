# Simulation figures

Run this in R or RStudio; no additional packages are needed:

```r
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/sim/plot_simulations.R")
```

The script reads `simulation results/chunks` and writes to `simulation results/figures`. It only reads saved results: it does not source the workhorse, refit models, or modify the simulation checkpoints. Change `project_dir` at the top if you move the project.

## Figures

Each main figure is saved as both PDF and PNG. Blue is SuSiE; pink is SuSiE-mix. Columns are PVE 10%, 20%, 30%, and 40%.

Both coverage figures share a y-axis ending at 1. The lower limit is calculated from the lowest coverage across the pure and mixed scenarios, with a small margin, rounded down to a multiple of 0.05. For the current full results, this gives 0.75 to 1.00. The 0.95 reference remains visible.

| Output | Rows | Within each panel |
|---|---|---|
| `coverage_pure`, `purity_pure`, `power_pure` | Additive only; dominant only; recessive only | Number of true causal SNPs from 1 to 5 |
| `coverage_mixed`, `purity_mixed`, `power_mixed` | Additive + dominant; additive + recessive; recessive + dominant; all three | Number of true causal SNPs from 2 to 5, or 3 to 5 for all three |
| `roc_pure_L1` through `roc_pure_L5` | Additive only; dominant only; recessive only | One causal SNP count per figure, with two solid curves per panel |
| `roc_mixed_L2` through `roc_mixed_L5` | The four mixed-scenario rows | One causal SNP count per figure, with two solid curves per panel |
| `roc_pure_by_L.pdf`, `roc_mixed_by_L.pdf` (optional) | The same respective scenario rows | Collections with one page per causal SNP count; enable `write_roc_pages` |

The true number of causal SNPs is `L_add + L_rec + L_dom`. It is different from the fitted upper bound `L = 10` in your jobs. Each ROC title specifies the true count. Each panel contains only the blue SuSiE curve and pink SuSiE-mix curve, with identical axis ranges across L. The script no longer generates the old overlaid `roc_pure` and `roc_mixed` figures; existing copies of those older files are not replaced.

For example, additive + recessive with three causal SNPs pools the `(L_add=1, L_rec=2)` and `(L_add=2, L_rec=1)` configurations. Every included simulation contributes; this is not an equal-weight average of configuration-level averages. If jobs have unequal completion or failure rates, their contributions differ. The saved configuration counts make that visible.

## Metric definitions

- **Coverage:** total number of reported credible sets containing at least one generating causal biological SNP, divided by the total number of reported sets. The dashed reference is 0.95. Runs with no CS contribute no sets to this denominator; they are not assigned zero coverage. This is not the saved posterior probability `sets$coverage`.
- **Purity:** mean of the saved `min.abs.corr` values across reported sets. Singleton purity is 1. For SuSiE-mix, these correlations are between its fitted predictor columns, which may use different codings; they are not recomputed as additive SNP LD after collapsing codings. No genotype matrices are available in the compact results to recompute an alternative purity measure.
- **Power:** fraction of generating causal SNPs appearing in the union of reported credible sets. Each causal SNP counts once even if multiple sets contain it. No-CS runs have zero power.
- **ROC:** causal SNP detection using PIPs, with TPR = TP / number of causal SNPs and FPR = FP / number of noncausal SNPs, pooled over replicates in each scenario/PVE/causal-count group. A discovery means PIP at least the threshold. FPR is different from false discovery rate. Counts are calculated at a dense common threshold grid; the displayed lines join these threshold points. Tied PIPs enter together. The full range is saved; the display shows FPR from 0 to 0.25, as in your reference.

Mixed credible sets are evaluated using the saved biological-SNP mapping. ROC uses `susie_mix_pip_snp`, which your workhorse already computes by combining codings within each single effect before calculating the SNP PIP. It does not sum final coding-specific PIPs or treat the three codings as independent SNPs.

The plotted points are pooled estimates, without confidence intervals. The scripts compare the saved standard additive and stacked mixed SuSiE fits; they do not change fitting settings or evaluate coding-identification accuracy.

## Checkpoints and incomplete results

The directory currently has 221 checkpoint files, including an older file for one condition. The script prefers the larger advertised checkpoint and removes repeated condition/seed pairs, so overlapping saves are not counted twice. Existing successful saved replicates are used even when fewer than the filename's requested number are present. Error entries are logged and omitted from both methods, never treated as null/no-discovery fits.

By default, successful fits with `converged=FALSE` remain in the comparison and their counts are reported. Set `exclude_nonconverged <- TRUE` to remove a replicate from both methods whenever either method did not converge. This preserves the paired comparison.

The script checks the saved settings against filenames, uses n = 500 and fitted L = 10 by default, and refuses to silently combine unexpected QC settings, all-additive paired controls, or missing SNP-level mixed PIPs. Blank panels are labelled when no matching saved results are available. The all-three scenario at true L = 2 is labelled not applicable in its ROC figure.

## Tables for checking the figures

- `metric_summary.csv`: pooled coverage, purity and power, together with their counts and the number of contributing replicates/configurations.
- `configuration_counts.csv`: number of included replicates per exact causal allocation, PVE and method.
- `file_audit.csv`: saved, examined, included, duplicate, error and convergence counts per checkpoint.
- `replicate_metrics.rds`: compact metrics per included replicate and method.
- `roc_counts.rds`: threshold-specific TP, FP, TPR and FPR for every plotted curve, including thresholds outside the displayed FPR range.

For a quick local preview, change `max_reps_per_file <- Inf` to 5, and preferably give previews a different `output_dir`. Restore `Inf` for the final figures. `write_png <- FALSE` writes only PDFs; `roc_max_fpr <- 1` displays the full ROC range.

Only a small smoke test was run during initial development: 3 records from each of 9 selected files, leaving 24 distinct simulations after duplicate removal. All seven scenario types were represented. The test also checked no-CS handling, tied PIPs, threshold endpoints and agreement with saved workhorse metrics. Test output is under `tmp/simulation_plot_smoke`, not the final figures directory. The development test did not run the full simulation collection.
