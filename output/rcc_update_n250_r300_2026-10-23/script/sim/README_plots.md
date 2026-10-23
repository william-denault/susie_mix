# Plots for the n=250, three-method benchmark

From the project root, run `source("generate_sim_figures.R")` in R or:

```sh
Rscript --vanilla generate_sim_figures.R
```

The entry point explicitly selects `simulation results/n250_r300_v1/chunks`,
n=250 and SuSiE, SuSiE-mix, SuSiE-slide-prior, then restores session options.
It never fits models. The old generate_sim_figures_slide_prior.R name forwards
to this entry point. Archived results have their own frozen plotting scripts.

## September-style panels

- Four PVE columns: **10%, 20%, 30%, 40%**.
- Causal count **1–5**; pairs start at 2 and triples at 3.
- Blue SuSiE, pink SuSiE-mix, green SuSiE-slide-prior, with offset point markers,
  right-side scenario labels and a dashed 95% coverage reference.
- Point-only coverage, purity, power and CS-size panels by default, as in the
  September screenshots. Analytic intervals and paired comparisons remain in
  the tables. Set `options(susie.sim.show_intervals=TRUE)` to display metric bars.
- Five pure rows, ten pair rows and ten triple rows, in groups `pure`,
  `pairs_1`, `pairs_2`, `triples_1`, `triples_2` with at most five rows each.
- Coverage starts at .75, or lower when necessary to show all estimates at
  n=250. Every panel within a figure shares limits. The new runs are not
  restricted to the K=1–4 display of the archived publication figures.

PDF and PNG outputs go to `simulation results/n250_r300_v1/figures/`.
There are coverage, purity, power and CS-size panels; ROC and power-FDR curves
for individual K and all K pooled; and pooled-K PIP calibration by method.
PIP calibration retains its +/-2 empirical-SE bars across seed blocks.

## Metrics and uncertainty

Coverage = reported CSs containing any true biological causal SNP / reported
CSs. No-CS runs add no sets to that denominator. Power = causal SNPs found in
the union of reported sets / all causal SNPs; no-CS runs contribute zero.
Purity is mean saved minimum absolute correlation under the fitted coding.
CS size is the mean number of unique biological SNPs per reported set.

SuSiE-mix PIPs aggregate coding probabilities within each single effect before
combining active effects. Mixed sets are mapped to biological SNPs. All methods
use the same included replicate identities. Error records are excluded;
nonconverged fits remain by default and are counted in the audit.

Counts are pooled before ratios. Marginal coverage/power/purity intervals use
the existing normal proportion approximation; CS-size intervals use the
empirical size variance. Three paired method contrasts use seed-block
delta-method intervals. ROC and empirical-FDR curves pool threshold counts,
retain threshold order, and have no uncertainty bands. See the helper code
for formulas; these intervals do not imply guaranteed FDR control.

Tables include file_audit.csv, configuration_counts.csv, metric_summary.csv,
method_comparison.csv, replicate_metrics.rds, ROC count caches and calibration
CSV/RDS caches. Check the audits for errors and incomplete conditions before
interpreting figures. After adding chunks, rerun the full entry point. For
an unchanged selection, `source("script/sim/refresh_simulation_figures.R")`
reuses compact saved summaries to redraw more quickly.
