# Simulation figures for publication

Open `generate_publication_plots.R` in RStudio and click **Source**. It finds
its own folder, so RStudio's working directory does not need to match. Alternatively,
from the folder containing the script, run:

```r
source("generate_publication_plots.R")
```

Figure 1 requires `ggplot2` and `patchwork`. If needed, install them once:

```r
install.packages(c("ggplot2", "patchwork"))
```

This folder is separate from the original analysis. All publication scripts,
including both Figure 1 entry scripts, are copies; the original scripts and
figures remain available for discussions with collaborators.

The script reads snapshots of the existing summaries in `input_data`.
It does not run simulations or refit models. After adding simulation results,
refresh the original summaries and replace the three files in `input_data`
(`metric_summary.csv`, `roc_counts.rds`, and `pip_calibration_seed_counts.rds`),
or pass the refreshed summary folder as `input_dir` to `generate_publication_figures()`.

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

Only true causal counts **K = 1–4** are included (pairs start at 2; triples at 3).
Both pooled calibration and pooled ROC/power–FDR curves are recomputed after
filtering K. Calibration retains its existing seed-block empirical error bars.
The fitted SuSiE upper bound L remains 10; it is not the true causal count K.

- `coverage_<group>`, `purity_<group>`, `cs_size_<group>`: value versus causal count.
- `power_<group>`: current power scale, 0–1.
- `power_<group>_zoom_075_1`: the same points/bars with the y-axis restricted to 0.75–1.
- `roc_<group>_L<K>`, `power_fdr_<group>_L<K>`: curves for each applicable K up to 4.
- `roc_<group>_all_L`, `power_fdr_<group>_all_L`: curves pooled across applicable K up to 4.
- `pip_calibration_<group>_all_L_<method>`: one method per calibration figure.
- `figure1/main_simulations`: the updated five-column manuscript figure.

Groups are `pure`, `pairs_1`, `pairs_2`, `triples_1`, and `triples_2`, retaining the
existing scenario arrangement and PVE columns. Missing results are labeled.

The printed normal-CI formula is removed. Coverage, power and CS-size error bars
retain their current definitions. Purity has no error bars and its y-axis label
is exactly `purity`. Points are offset by −0.24, 0, +0.24 around each causal count;
Figure 1 uses the equivalent dodge width of 0.72. A zoomed power plot can omit
estimates below 0.75; the full-scale version includes them.

For editable settings, use `script/sim/plot_publication_simulations.R` and
`script/sim/main_manuscript_fig1.R`. The alternate Figure 1 entry
`script/sim/plot_manuscript_simulations.R` has the same changes. Shared method
labels and colors live in `script/sim/publication_plot_helpers.R`.

The saved CSV/RDS tables and `publication_settings.rds` record the exact plotted
data, method mapping, causal-count selection, offsets, axis ranges, and source hashes.
