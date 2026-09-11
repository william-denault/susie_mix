# Coding classification and mixture recovery

Run the new script in R or RStudio:

```r
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/sim/plot_coding_mixture.R")
```

Or, from the project root:

```bash
Rscript --vanilla script/sim/plot_coding_mixture.R
```

It reads the existing `simulation results/chunks` checkpoints and writes eight
figures as PDF and PNG, plus CSV summaries, into
`simulation results/coding_figures`. It uses only base R. It does not run
simulations, refit SuSiE, or modify the existing workhorse, plotting script,
or checkpoints. The default project path matches this project on Windows and
the cluster; `SUSIE_MIX_PROJECT_DIR` overrides it.

## What your simulations already provide

`sim_workhorse.R` uses real genotype loci, simulates known causal SNPs and
codings, and saves both coding-specific PIPs and biological-SNP PIPs. The
current job generator includes every allocation of 1–5 causal effects among
additive, recessive, and dominant coding, at PVE 10%, 20%, 30%, and 40%, with
sample size 500. In particular, the additive-only conditions are already a
control for apparent nonadditive support under additive truth.

`plot_simulations.R` assesses biological-SNP localization after combining
codings. This new script assesses the coding itself and the pooled mixture
estimate proposed for EM.

## Figures and how to read them

| Figure | Question and interpretation |
|---|---|
| `additive_only_pip_allocation` | How much additive, recessive, and dominant PIP mass appears when every generating effect is additive? The truth is additive share 1 and the other shares 0. Nonzero nonadditive mass is the simulation baseline for coding ambiguity/false support. |
| `mixture_recovery` | Does the pooled PIP mixture recover the known fraction of causal effects of each coding? Each point is one causal allocation; columns show PVE, rows show coding, symbols show the true number of causal SNPs K. The dashed line is exact recovery. Bars show 95% seed-block bootstrap intervals. |
| `coding_confusion` | At each known causal SNP, which coding has the highest PIP? Rows are true coding; columns are called coding, ties, and undetected SNPs. Every causal SNP stays in the denominator. This is classification given the true SNP location. |
| `coding_classification` | How often is that coding call correct as K increases? The upper row uses all causal SNPs; the lower row conditions on detected causal SNPs. Colors indicate the generating coding. |
| `coding_pip_calibration` | For predictors assigned a given PIP, how often is that exact SNP/coding pair a generating effect? The diagonal is empirical agreement. Points below it indicate overprediction in the selected simulation conditions. Larger points have more predictors. |
| `coding_pip_calibration_additive_only` | The same diagnostic under additive truth alone. Recessive/dominant predictors are all negative labels; positive PIPs on these codings make the ambiguity/false-support baseline explicit. |
| `coding_discovery_accuracy` | Among all coding-specific calls with PIP at least 0.9, what fraction have the exact correct SNP AND coding (precision)? What fraction of the true effects of each coding are recovered (recall)? This includes errors at other SNPs as well as wrong codings at causal SNPs. |
| `additive_only_false_calls` | What fraction of purely additive simulations have at least one recessive or dominant call with PIP at least 0.9 anywhere in the locus? This measures confident false nonadditive calls separately from diffuse PIP mass. |

The default detection cutoff is 0.9. Classification uses the saved
`susie_mix_pip_snp` to decide whether the **causal biological SNP** is detected,
then selects the largest coding-specific PIP at that SNP. A tie within
`1e-12` is called ambiguous; an undetected SNP or all-zero coding scores are
called undetected. Missing coding predictors receive zero for classification.
The largest PIP is used only as a ranking score: normalized coding PIPs are
**not** treated as mutually exclusive posterior coding probabilities.

Classification at the known causal SNP is easier than locating the SNP and
its coding. Read the confusion/classification plots together with the
discovery-accuracy plot and the original SNP-localization figures.

## Mixture estimate and uncertainty

For each exact simulation condition, the estimated coding share is:

```text
sum of PIPs for this coding over all included replicates and predictors
--------------------------------------------------------------------
sum of PIPs over all three codings, replicates, and predictors
```

This is the same pooling rule used in the EM preparation. It includes every
predictor, without a PIP or credible-set filter. It does not average
replicate-normalized shares. The true share is the number of generating
effects of that coding divided by the total number of generating effects.
The script reports bias relative to that count-based target; this is a
diagnostic of the proposed estimator, not an assumption that PIP pooling is
unbiased for effect proportions.

Mixture intervals resample whole seeds. Your generator reuses seeds across
PVE and coding allocations, so all comparisons share the same bootstrap
weights. The default is 500 resamples and percentile 95% intervals. These
describe simulation sampling uncertainty, not uncertainty for the real-data
EM estimates. Conditions with fewer than two distinct seeds have no interval.
Zero-total-PIP replicates are counted; if an entire condition has zero mass,
its estimated shares are undefined rather than replaced with an arbitrary prior.

Other figures show pooled empirical rates without uncertainty bars. Their
counts and denominators are exported. Pooling is weighted by available
replicates, true effects, or predictors as appropriate; it is not an
equal-weight average across causal allocations. Larger loci contribute more
predictors to PIP calibration. The simulation grid fixes causal counts and
effect sizes; these calibration plots describe that grid, rather than a
sample drawn from the model's full prior distribution.

## Saved tables

- `mixture_recovery.csv`: true/estimated shares, bias, PIP mass, bootstrap
  bounds, replicate and seed counts for every exact condition and coding.
- `coding_confusion.csv`: counts by scenario, PVE, K, true coding, and called
  coding. Zero-count cells are retained.
- `coding_classification.csv`: all-causal and detected-causal denominators
  and accuracy, pooled across allocations by PVE/K/generating coding.
- `coding_pip_calibration.csv`: ten PIP bins per coding and scenario/PVE/K,
  counts, summed PIP, true-positive count, and mean Brier score. Bins are
  `[0,.1)`, ..., `[.9,1]`, with exact PIP 1 in the final bin.
- `coding_discovery_accuracy.csv`: precision, recall, and pooled false
  discovery proportion at PIP cutoffs 0, .1, .5, .8, .9, .95, .99, and 1.
  False discoveries are split into wrong coding at a true causal SNP versus
  a noncausal biological SNP. No-call precision/FDP and no-truth recall are
  undefined (NA), not assigned perfect scores or zero.
- `additive_only_false_call_rates.csv`: fractions of additive-only
  replicates with at least one incorrect call of each coding at the chosen
  cutoff. For additive coding this means a call at a noncausal SNP.
- `configuration_counts.csv`: included replicates for each condition.
- `file_audit.csv` and `excluded_replicates.csv`: inclusion, duplicate,
  error, validation, and convergence records.
- `replicate_coding_metrics.rds`: compact replicate-level PIP masses, truth
  counts, classification counts, and false-call indicators for further work.
- `plot_settings.rds`: the options used for the output.

## Compact checkpoint compatibility and exclusions

The loader prefers explicit `mix_coding` labels or coding suffixes in PIP
names when present. Current compact saves omit both. In that case, it
reconstructs the additive/recessive/dominant blocks from `mix_to_add`, using
the block order in `sim_workhorse.R` and the fact that SNP indices increase
within each block. The generating coding and predictor indices must agree
with the reconstruction. Different block sizes after QC are handled.

If more than one recessive/dominant boundary fits the saved information,
the replicate is excluded and reported; the script never assumes equal
block lengths or guesses the boundary. Explicit coding labels can be saved
in future simulations to handle such cases. Existing checkpoints are not changed.

Saved PIPs outside `[0,1]` by at most `1e-10` are clipped to the boundary for
analysis and counted in `roundoff_clamped_replicates`. This handles harmless
floating-point error when grouping SNP probabilities. Larger excursions or
nonfinite PIPs are excluded as invalid; the checkpoint itself is not edited.

Overlapping checkpoints are deduplicated by the original causal allocation,
`all_additive` flag, PVE, sample size, fitted L, and seed, preferring the file
requesting more replicates and then the newer copy. Failed/invalid copies do
not prevent a usable duplicate from being included. All-additive controls
use their actual generating truth, and retain their original SNP-selection
design in the mixture table so they are not mistaken for an independent copy
of the ordinary additive-only condition.

Nonconverged mixed fits are included and counted by default. Set
`exclude_nonconverged = TRUE` to exclude them. Additive-only fit convergence
does not affect these mixed-fit diagnostics. Review audit and configuration
counts when checkpoint coverage is incomplete or exclusions are substantial.

## Scope relative to the genome-wide EM analysis

These are simulations of the existing **unweighted** SuSiE-mix fit. They do
not validate convergence of the new EM loop or its final weighted estimates.
The simulation uses `standardize=TRUE`, prior-variance method `optim`, and
Gaussian phenotypes; the scan uses `standardize=FALSE`, prior-variance method
`EM`, and processed expression phenotypes. Simulations also enforce eligible
causal genotypes, equal individual contribution variances, and the selected
PVE/sample-size grid. Use them to assess coding ambiguity and estimator
behavior, while keeping these differences in mind when comparing numerical
shares to real tissues.

For a quick preview or a different cutoff, load the functions without running:

```r
options(susie_mix.coding_plots.run = FALSE)
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/sim/plot_coding_mixture.R")
config <- coding_plot_settings
config$max_reps_per_file <- 5
config$output_dir <- file.path(project_dir, "tmp/coding_preview")
preview <- run_coding_plots(config)
```

Run the synthetic validation checks with:

```bash
Rscript --vanilla script/sim/tests/test_coding_mixture.R
```
