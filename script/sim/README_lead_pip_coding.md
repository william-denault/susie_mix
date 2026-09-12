# Exact discovery by lead PIP

Run from the project root:

```bash
Rscript --vanilla script/sim/plot_lead_pip_coding.R
```

Or in R/RStudio:

```r
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/sim/plot_lead_pip_coding.R")
```

The script reads existing `simulation results/chunks` checkpoints and writes
four figures as PDF and PNG, summary CSV files, and compact RDS files to
`simulation results/lead_pip_figures`. It requires only base R and does not
run simulations or fine mapping. `SUSIE_MIX_PROJECT_DIR` overrides the default
project location on Windows or the cluster.

## What is selected and what counts as correct?

Each simulation contributes **one global lead**, selected by the largest
coding-specific PIP across every SNP and all three codings. Selection does
not use the causal truth, credible sets, biological-SNP PIPs, or a detection
threshold. Truth is used to score the selected predictor:

- **Exact SNP + coding:** the lead is one of the generating SNP/coding pairs.
- **Wrong coding at causal SNP:** the SNP is causal, but its selected coding
  is not its generating coding.
- **Noncausal SNP:** the selected SNP is not causal, regardless of its coding.
- **Tied leads:** more than one predictor is within `tie_tolerance` of the
  maximum PIP (default absolute tolerance `1e-12`). Ties may involve different
  codings, different SNPs, or several true effects; they are marked ambiguous
  even if all tied candidates are true. The script never breaks a tie using
  column order or simulation truth.
- **All PIPs zero:** no positive lead exists.

Very high PIPs can round to the same saved value. These leads are also tied;
the script cannot recover an ordering lost from the saved PIPs. The replicate
table retains the number of tied candidates, number of exact true pairs among
them, and the PIP gap between the two highest-scoring predictors.

For multiple causal effects, an exact lead can match **any** generating pair.
This measures the strongest reported hit, not recovery of every signal.
One lead per credible set is a different analysis and is not performed here.

## Figures

| File | Interpretation |
|---|---|
| `lead_pip_outcomes` | Outcome proportions against lead PIP, across all simulated coding allocations. Rows are causal count K; columns are PVE. Each occupied PIP bin is a stacked bar whose denominator includes every simulation in that bin, including ties and zero support. Counts above bars show the denominator. |
| `lead_pip_outcomes_additive_only` | The same analysis restricted to simulations whose generating effects are all additive. Correct leads must be additive at a causal SNP. |
| `lead_pip_exact_accuracy` | Among positive, unambiguous leads of each **called** coding, how often is the exact SNP/coding pair true? The x-coordinate is the mean lead PIP in a bin; the y-coordinate is the observed exact-match fraction. Point size increases with bin count. The diagonal indicates equality. |
| `additive_only_lead_coding` | For each PVE and K, the count of unambiguous leads using each coding divided by the number of positive, unambiguous leads under additive truth. All SNPs compete. Panel annotations report unique leads / all simulations; detailed exclusions are in the tables. |

The default bins are `[0,0.1)`, ..., `[0.9,1]`; exact PIP 1 belongs in the
last bin. Empty bins are blank and are omitted from the bin tables. Neither
empty bins nor zero-denominator accuracy are assigned perfect performance.

In the additive-only coding-share plot, a dominant lead is counted whether
it is at a causal SNP or another SNP. This answers how often dominance wins
the ranking. It differs from the pooled dominant PIP mass / total PIP mass
used for the proposed EM update. Both quantities can be examined using these
figures and `plot_coding_mixture.R`.

## Tables and reproducibility

- `lead_pip_outcomes_by_bin.csv`: counts and fractions for all five outcomes,
  mean lead PIP, called-coding counts, and bin boundaries, by scope/PVE/K/bin.
- `lead_pip_accuracy_by_coding.csv`: exact-match counts, number of unique
  leads, mean lead PIP, and exact accuracy by scope/PVE/K/called coding/bin.
  It includes both all-simulation and additive-only scopes.
- `lead_pip_summary.csv`: outcomes and lead-coding shares by scope/PVE/K.
- `lead_pip_outcomes_by_condition.csv`: the same quantities separately for
  every original causal allocation and all-additive control design.
- `additive_only_lead_coding.csv`: additive-only rows of the summary,
  including `n_runs`, `n_unique`, `n_ambiguous`, `n_no_support`, and
  `share_unique_additive`, `share_unique_recessive`, `share_unique_dominant`.
- `lead_pip_replicates.rds`: one row per included simulation with seed,
  original and actual generating configuration, lead index, biological-SNP
  index, called coding, lead PIP, outcome, and tie diagnostics. Predictor/SNP
  indices refer to the corresponding saved replicate; ties have no selected
  index or called coding.
- `lead_pip_analysis.rds`: replicate data, summaries, and audit for redrawing
  figures without rereading all raw checkpoints.
- `plot_settings.rds`: settings used for the analysis.
- `file_audit.csv` and `excluded_replicates.csv`: checkpoint inclusion,
  duplicates, failed or invalid simulations, convergence, and PIP roundoff.

The shared loader in `coding_mixture_utils.R` validates coding maps and truth,
checks simulation/QC settings, reconstructs unequal coding blocks when labels
are missing, and deduplicates overlapping checkpoints. Failed/invalid copies
do not suppress a usable duplicate. Actual generating labels determine
whether a run is additive-only, including explicit all-additive controls.
Roundoff-sized PIP excursions outside `[0,1]` (at most `1e-10`) are clipped
for analysis and audited; larger errors are excluded. Mixed-fit nonconvergence
is included and counted by default; it can be excluded in the settings.

Within each PVE/K panel, observations are pooled by simulation, not by the
number of predictors or true causal effects. Different causal allocations
and any explicit additive-control designs contribute according to their
available replicates; their separate counts are retained in the condition
table. Figures show descriptive rates and denominators, without confidence
intervals. Seeds are reused across conditions, so plotted observations across
panels should not be treated as independent draws for inference.

These results describe the existing unweighted simulation fits. They do not
validate the weighted EM loop, and their numerical rates are not a direct
calibration of the tissue scan, which uses different fitting/preprocessing
settings. See `README_coding_mixture.md` for those differences.

## Customization and fast redraw

```r
options(susie_mix.lead_pip_plots.run = FALSE)
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/script/sim/plot_lead_pip_coding.R")
config <- lead_pip_settings
config$max_reps_per_file <- 5    # Optional preview; whole checkpoint files still load.
config$output_dir <- file.path(project_dir, "simulation results/lead_pip_preview")
analysis <- run_lead_pip_plots(config)
```

After a full run, reload the saved analysis and settings to redraw:

```r
config <- readRDS(file.path(lead_pip_settings$output_dir, "plot_settings.rds"))
analysis <- readRDS(file.path(config$output_dir, "lead_pip_analysis.rds"))
lp_render_plots(analysis, config)
```

Metric, checkpoint, and plotting checks run with:

```bash
Rscript --vanilla script/sim/tests/test_lead_pip_coding.R
Rscript --vanilla script/sim/tests/test_coding_mixture.R
```
