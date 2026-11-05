# Simulation figures for publication

Open `generate_publication_plots.R` in RStudio and click **Source** to refresh
the selected publication panels and Figure 1. It finds its own folder, so RStudio's
working directory does not need to match. From the repository root, this is:

```r
source("publication_plots/generate_publication_plots.R")
```

To load the functions first and run everything explicitly from any working directory:

```r
options(susie.publication.autorun = FALSE)
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/publication_plots/generate_publication_plots.R")
figures <- generate_publication_figures()
```

`figures$figure1$figure` is the editable Figure 1 patchwork object.

Figure 1 requires `ggplot2` and `patchwork`. If needed, install them once:

```r
install.packages(c("ggplot2", "patchwork"))
```

This folder is separate from the original analysis. All publication scripts,
including both Figure 1 entry scripts, are copies; the original scripts and
figures remain available for discussions with collaborators.

The entry point reads the **original n=500 publication snapshots** in
`input_data/`, restoring the old comparison with **causal counts 1–4**.
PVE is **10%, 20%, 30%, 40%**; 5% remains excluded as requested.
It reads the three saved summaries (`metric_summary.csv`, `roc_counts.rds`, and
`pip_calibration_seed_counts.rds`), without rebuilding from raw checkpoints.

These snapshots include the finite-prior results saved as `SuSiE-slide-prior`.
The older raw `../simulation results/slide_v1/chunks/` files contain the
continuous-slide method, so they cannot replace the finite-prior snapshots.
The copied simulation design refers to the older `slide_prior_v1` run.
The n=250 results and summaries remain in their existing folder for later use;
they are no longer inputs to the publication entry point. The n=500 input check
prevents silently using them with old sample-size labels.

Outputs are PDF and PNG files in **`figures/`**, with Figure 1
in **`figures/figure1/`**. Run the script again to refresh them.
The original analysis figures are kept in their existing folder.

| Saved method identifier | Figure label | Color | Point |
|---|---|---|---|
| SuSiE | SuSiE | Blue `#2489FF` | Circle |
| SuSiE-mix | SuSiE-MiX | Red/pink `#D81B60` | Triangle |
| SuSiE-slide-prior | SuSiE-slide | Bluish green `#009E73` | Diamond |

The green comes from the Okabe–Ito palette. Distinct shapes and line styles also
identify methods. The continuous-delta method originally saved as `SuSiE-slide`
is excluded: the publication label **SuSiE-slide always refers to SuSiE-slide-prior**.
Data exports retain the original identifiers to avoid ambiguity.

True causal counts **K = 1–4** are included (pairs start at 2; triples at 3).
Both pooled calibration and pooled power–FDR curves are recomputed after
filtering K. Calibration retains its existing seed-block empirical error bars.
The fitted SuSiE upper bound L remains 10; it is not the true causal count K.

- `coverage_<group>`, `purity_<group>`, `cs_size_<group>`: value versus causal count.
- `power_<group>`: current power scale, 0–1.
- `power_<group>_zoom_075_1`: the same points/bars with the y-axis restricted to 0.75–1.
- `power_fdr_<group>_all_L`: curves pooled across applicable K up to 4.
  The title omits the causal-count pooling suffix; the manuscript captions
  describe the pooling. Per-K power–FDR figures and all ROC figures are omitted.
- `pip_calibration_pure_all_L_<method>`: one method per calibration figure,
  restricted to architectures with a single coding type. Mixed-architecture
  calibration figures are omitted.
- `figure1/main_simulations`: the updated six-column manuscript figure. Column F
  shows mean credible-set size (unique biological SNPs per reported CS) at 30% PVE,
  against 1–4 causal SNPs, with the existing CS-size error bars and one shared
  y-axis range across the three scenario rows.
  The PIP threshold markers use `>=0.95` and `>=0.80`, matching the saved counts.

Groups are `pure`, `pairs_1`, `pairs_2`, `triples_1`, and `triples_2`, retaining the
existing scenario arrangement and the four PVE columns. Missing results are labeled.

The printed normal-CI formula is removed. Coverage, power and CS-size error bars
retain their current definitions. Purity has no error bars and its y-axis label
is exactly `purity`. Points are offset by −0.24, 0, +0.24 around each causal count;
Figure 1 uses the equivalent dodge width of 0.72. A zoomed power plot can omit
estimates below 0.75; the full-scale version includes them.
Coverage uses the original 0.70–1.00 scale.

For editable settings, use `script/sim/plot_publication_simulations.R` and
`script/sim/main_manuscript_fig1.R`. The alternate Figure 1 entry
`script/sim/plot_manuscript_simulations.R` has the same changes. Shared method
labels and colors live in `script/sim/publication_plot_helpers.R`.

The saved CSV/RDS tables and `publication_settings.rds` record the exact plotted
data, method mapping, causal-count selection, offsets, axis ranges, and source hashes.

The manuscript supplement contains 33 selected simulation figures plus the
slider illustration (S1–S34). Older unselected outputs remain available for
collaborator discussions but are not regenerated or included in the manuscript.
