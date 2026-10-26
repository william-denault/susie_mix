# Simulation figures for publication

Open `generate_publication_plots.R` in RStudio and click **Source** to refresh
all publication panels and Figure 1. It finds its own folder, so RStudio's
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

`figures$figure1$figure` is the editable Figure 1 patchwork object;
`figures$audit` contains the checkpoint and replicate counts.

Figure 1 requires `ggplot2` and `patchwork`. If needed, install them once:

```r
install.packages(c("ggplot2", "patchwork"))
```

This folder is separate from the original analysis. All publication scripts,
including both Figure 1 entry scripts, are copies; the original scripts and
figures remain available for discussions with collaborators.

The entry point reads the synced checkpoints in
`../simulation results/n250_r300_v1/chunks/`: **n=250**, **300 requested replicates
per condition**, all **225 effect allocations**, three methods, and **PVE 10%,
20%, 30%, 40%**. Temporary `.tmp` and `.filepart` files are not read.
It rebuilds the summaries in `../simulation results/n250_r300_v1/figures/`
before drawing the publication plots. It does not run simulations or refit models.
The old `input_data` snapshots are no longer the default input.

The first run reads all checkpoints (about 23 GB for the complete run).
Subsequent calls reuse summaries only when the checkpoint names, sizes and
modification times, summary-script hashes, and saved summary hashes match the
manifest. Syncing added or updated checkpoints automatically triggers a rebuild.
To force a rebuild, use `generate_publication_figures(rebuild_summaries = TRUE)`.
Finish syncing before running: a changed input signature during summarization
stops the run instead of marking inconsistent summaries as current.

The function prints the file count, included replicate count, error count, and
number of checkpoints with fewer than 300 recorded attempts. Error records are
excluded from all three methods together; nonconverged fits remain included.
Review `file_audit.csv` and `configuration_counts.csv` in the summary folder.
Three hundred recorded attempts do not necessarily mean 300 successful replicates.

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

True causal counts **K = 1–5** are included (pairs start at 2; triples at 3).
Both pooled calibration and pooled ROC/power–FDR curves are recomputed after
filtering K. Calibration retains its existing seed-block empirical error bars.
The fitted SuSiE upper bound L remains 10; it is not the true causal count K.

- `coverage_<group>`, `purity_<group>`, `cs_size_<group>`: value versus causal count.
- `power_<group>`: current power scale, 0–1.
- `power_<group>_zoom_075_1`: the same points/bars with the y-axis restricted to 0.75–1.
- `roc_<group>_L<K>`, `power_fdr_<group>_L<K>`: curves for each applicable K up to 5.
- `roc_<group>_all_L`, `power_fdr_<group>_all_L`: curves pooled across applicable K up to 5.
- `pip_calibration_<group>_all_L_<method>`: one method per calibration figure.
- `figure1/main_simulations`: the updated six-column manuscript figure. Column F
  shows mean credible-set size (unique biological SNPs per reported CS) at 30% PVE,
  against 1–5 causal SNPs, with the existing CS-size error bars and one shared
  y-axis range across the three scenario rows.
  The PIP threshold markers use `>=0.95` and `>=0.80`, matching the saved counts.

Groups are `pure`, `pairs_1`, `pairs_2`, `triples_1`, and `triples_2`, retaining the
existing scenario arrangement and the four PVE columns. Missing results are labeled.

The printed normal-CI formula is removed. Coverage, power and CS-size error bars
retain their current definitions. Purity has no error bars and its y-axis label
is exactly `purity`. Points are offset by −0.24, 0, +0.24 around each causal count;
Figure 1 uses the equivalent dodge width of 0.72. A zoomed power plot can omit
estimates below 0.75; the full-scale version includes them.
Coverage keeps its usual 0.70 lower limit when possible, and expands downward
in 0.05 steps when needed to show lower estimates in the new n=250 results.

For editable settings, use `script/sim/plot_publication_simulations.R` and
`script/sim/main_manuscript_fig1.R`. The alternate Figure 1 entry
`script/sim/plot_manuscript_simulations.R` has the same changes. Shared method
labels and colors live in `script/sim/publication_plot_helpers.R`.

The saved CSV/RDS tables and `publication_settings.rds` record the exact plotted
data, method mapping, causal-count selection, offsets, axis ranges, and source hashes.
