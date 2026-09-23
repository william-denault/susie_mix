# Simulation figures

Run `source("script/sim/plot_simulations.R")` from R at the project root, or use the script's absolute path. Only base R is needed. Set `SUSIE_MIX_PROJECT_DIR` or edit `project_dir` if the project moves.

The script reads `simulation results/slide_v1/chunks` and writes to `simulation results/slide_v1/figures`. It never refits models. See [README_simulations.md](README_simulations.md) for the generating design and cluster commands.

To redraw an existing analysis after changing figure settings, use
`source("script/sim/refresh_simulation_figures.R")`. This reuses the saved CS/ROC
summaries. The first refresh of older summaries reads the checkpoints once to
calculate exact PIP calibration bin totals; later refreshes reuse those totals.
Use the full `plot_simulations.R` after adding results or changing inclusion
settings, so all summaries use the new selection.

## Figures

Blue is SuSiE, pink is SuSiE-mix, and green is SuSiE-slide. Columns are PVE **5%, 10%, 20%, 30%, 40%**. Figures contain at most five scenario rows:

| Suffix | Scenarios |
|---|---|
| `pure` | All five pure effect types |
| `pairs_1`, `pairs_2` | All ten pairs of distinct effect types |
| `triples_1`, `triples_2` | All ten triples of distinct effect types |

Effect order is additive, recessive, dominant, partial recessive, partial dominant. Pairs and triples follow combination order. Every allocation with the indicated total causal count K contributes to its scenario. For example, recessive + partial recessive at K = 3 pools allocations (1, 2) and (2, 1). K ranges from 1–5 for pure effects, 2–5 for pairs, and 3–5 for triples; fitted L remains 10.

Each figure is saved as PDF and PNG by default:

- `coverage_<group>`, `purity_<group>`, `power_<group>`, `cs_size_<group>`: points and 95% bootstrap intervals against the number of causal SNPs.
- `roc_<group>_L<K>` and `power_fdr_<group>_L<K>`: three method curves at a fixed true causal count.
- `roc_<group>_all_L` and `power_fdr_<group>_all_L`: all applicable causal counts pooled into one curve per method. Set `write_roc_all_L <- FALSE` to omit.
- Optional `roc_<group>_by_L.pdf` and `power_fdr_<group>_by_L.pdf`: multipage collections when `write_roc_pages <- TRUE`.
- `pip_calibration_<group>_all_L` and `pip_calibration_<group>_L<K>`: mean PIP
  against the fraction of SNPs that are causal, for each method and PVE. The
  same `write_roc_all_L` and `write_roc_pages` settings control pooled and
  optional multipage calibration exports.

Each scenario/PVE panel has its own labelled y-axis. Coverage and purity run
from the lowest point estimate to 1, ignoring confidence bounds when choosing
the limits. CS size runs from 0 to the highest point estimate. Power-FDR runs
from 0 to the highest curve value within the displayed FDR window, including
line intersections with the window boundary. Confidence intervals can be
clipped by these limits. Power, ROC and calibration retain a 0-1 y-axis;
calibration also uses a 0-1 x-axis. The 0.95 coverage reference is shown when
inside the panel. An all-one coverage/purity panel uses 0.95-1 to avoid a
zero-height plot; empty/all-zero panels use 0-1. Missing cells say "No saved
results". PNGs use Cairo when available to preserve labels on Windows.

## Metric definitions

- **Coverage:** reported credible sets containing a true biological causal SNP, divided by the number of reported sets. No-CS runs contribute no sets to this denominator. This is not posterior `sets$coverage`.
- **Purity:** mean saved minimum absolute correlation across reported CSs. SuSiE-mix uses its fitted coding columns; SuSiE-slide uses each component's fitted transformed genotypes. These are each model's native purity measures, not a common additive-LD measure. Singleton purity is 1.
- **Power:** fraction of true causal SNPs found in the union of reported CSs. Each causal SNP counts once; no-CS runs have zero power.
- **CS size:** mean number of unique biological SNPs per reported CS. Multiple coding columns of one SNP count once within a mixed CS. No-CS runs contribute no fictitious size of zero.
- **ROC:** pooled SNP-level PIP threshold counts, TPR = TP / causal SNPs and FPR = FP / noncausal SNPs. Tied PIPs enter together. The full threshold range is saved; the display defaults to FPR up to 0.25.
- **Power–FDR:** pooled PIP-based TPR against FP / (TP + FP). This is pooled empirical FDP, labelled empirical FDR; it differs from FPR and the mean of replicate FDPs. No discoveries give (0, 0). Curves retain threshold order, including reversals, rather than a sorted or optimized envelope.

SuSiE-mix SNP PIPs combine coding probabilities **within each single effect** before computing the union across active effects. Final coding PIPs are not summed. All three methods are evaluated at biological SNPs, including partial-effect scenarios without an exact generating predictor in SuSiE-mix. Slider PIPs and CSs condition on fitted deltas.

Pooling sums counts before calculating ratios, including pooled-K curves. Larger loci contribute more noncausal SNPs, and allocations with more completed replicates contribute more observations. The analysis does not equally weight configuration averages; `configuration_counts.csv` exposes imbalance.

## Intervals and paired comparisons

The bootstrap resamples **simulation seeds**, carrying all their methods, configurations and PVE values together. Defaults are 1,000 resamples and 95% percentile intervals. Every resample recomputes pooled ratios. The reader and summary code require identical replicate identities for all selected methods.

`method_comparison.csv` reports all three paired contrasts: SuSiE-mix minus SuSiE, SuSiE-slide minus SuSiE, and SuSiE-slide minus SuSiE-mix. The `method` and `reference` columns specify subtraction direction. Positive CS-size differences mean larger sets, not improved resolution. Overlapping marginal intervals are not a paired comparison.

Intervals describe Monte Carlo uncertainty within this design. Fewer than two seed blocks gives unavailable bounds. Undefined bootstrap ratios are omitted for that metric and valid draw counts are saved. At boundaries intervals can collapse. ROC and power–FDR figures have no uncertainty bands.

## PIP calibration

The calibration definition follows Supplementary Figure 1 of the mvSuSiE
supplement (DOI `10.1038/s41588-025-02486-7`). Biological SNP PIPs from all
included replicates are assigned to ten equal-width bins: `[0,0.1)`, ...,
`[0.9,1]`. The x-coordinate is the **actual mean PIP**, not the bin midpoint.
The y-coordinate is the number of causal SNPs divided by the number of SNPs
in that bin. Zero and one are included; empty bins are omitted from the plot
and recorded with zero counts and unavailable estimates in the tables.

All-L figures pool SNP counts and PIP sums over the available causal counts
and allocations before calculating ratios; they do not average the separate
L-specific frequencies. Thus incomplete allocations contribute fewer SNPs.
SuSiE-mix uses its SNP-level PIPs, with each biological SNP counted once per
replicate. Points below the diagonal indicate PIPs larger than the observed
causal frequency; points above it indicate smaller PIPs.

Error bars are the observed frequency plus/minus **two empirical standard
errors**, clipped to 0-1, matching the reference's display convention. Here the
SE of the pooled ratio is calculated across seed blocks, keeping SNPs within
a replicate and repeated uses of a seed together. If `C_s` and `N_s` are a
seed's causal and total bin counts, `q = sum(C_s)/sum(N_s)` and
`SE = sqrt(S/(S-1) * sum((C_s-q*N_s)^2)) / sum(N_s)`. Seeds with no SNPs in
that bin have zero counts and remain in the calculation. Fewer than two
seeds gives unavailable error bars. These are approximate empirical error
bars, not independent-SNP binomial intervals or posterior credible intervals.

## Checkpoints, audits and tables

The reader checks all five causal counts, true effect labels/deltas, PVE, sample size, fitted L, QC settings and slider settings. It requires three-method saves by default. Partial effects are never coerced into endpoint truth categories. Paired all-additive controls need separate plots.

Larger advertised checkpoints are preferred; repeated configuration/seed pairs are removed. Successful saved replicates are included even when jobs are incomplete. Error records are audited and excluded from every method. By default, nonconverged fits remain in the comparison and are counted. Set `exclude_nonconverged <- TRUE` to exclude a replicate from every method whenever any fit did not converge. All-error input stops after writing the audit.

- `file_audit.csv`: saved, examined, included, duplicate, error and convergence counts.
- `configuration_counts.csv`: included replicates per exact allocation/PVE/method.
- `metric_summary.csv`: pooled metrics, interval bounds and contributing counts.
- `method_comparison.csv`: paired method-minus-reference differences and intervals.
- `replicate_metrics.rds`: per-replicate compact metric counts.
- `roc_counts.rds` and `roc_counts_all_L.rds`: threshold counts/rates by K and pooled K.
- `pip_calibration.csv` and `pip_calibration_all_L.csv`: bin counts, causal
  counts, PIP sums/means, observed frequencies, empirical SEs and error bars.
- `pip_calibration_seed_counts.rds`: exact per-seed bin totals, source file
  signatures and the included replicate/method selection for rapid redraws.

Run `Rscript script/sim/tests/test_plot_calibration.R` for base-R checks of
bin boundaries, pooled means/frequencies, seed-level uncertainty and panel
limits; this does not fit models.

For a small preview, set `max_reps_per_file <- 5`, use another output directory and optionally disable PNG output. Restore `Inf` for final figures. `tests/test_slide_simulation.R` validates all 25 families and exports representative figures to `tmp/slide_simulation_validation/figures`.

Legacy two-method saves remain in their original directory. To inspect them, explicitly set `chunk_dir` to `simulation results/chunks`, use a separate output directory, set `method_names <- c("SuSiE", "SuSiE-mix")`, and use the first two colors. The filename reader accepts the old format. Legacy saves cannot supply a three-method comparison. The exact-coding diagnostic scripts still target the original three-coding data; see README_simulations.md.
