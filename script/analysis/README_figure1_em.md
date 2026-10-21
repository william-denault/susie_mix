# Figure 1: count tables, ECDFs and posterior coding support

Open `script/analysis/plot_figure1_em.R` in RStudio on RCC and click **Source**:

```r
source("/project2/mstephens/wdenault/susie_mix/script/analysis/plot_figure1_em.R")
```

Sync the updated entry script and `figure1_em_utils.R` into the existing project's
`script/analysis/` directory together. They reuse `em_analysis_utils.R`,
`em_cs_agreement_utils.R`, `generate_summary_results.R` and
`script/scan_tissue_attempt/em_utils.R`. The existing `plot_cs_agreement_em.R`
remains the original prior-based figure; run **plot_figure1_em.R** for this revision.

This version uses **base R only**. It reads saved EM fits and summaries; it does
not refit models, run PLINK or load GTEx expression/genotype files.

## Settings and inputs

The editable settings block defaults to:

- Project: `/project2/mstephens/wdenault/susie_mix`. Edit `project_dir` or set
  `SUSIE_MIX_PROJECT_DIR` if needed.
- `iteration = "latest"`: resolved once to the latest completed EM iteration.
  Use `18L` to pin iteration 18.
- Minimum marginal **additive P < 1e-8**, **mean reads >=100**, **>=1 CS in either
  model**. Missing filter values and missing EM fits do not enter the cohort.
- `minimum_tissue_regions = 100L`: only the threshold for choosing the four
  tissue tables. It does **not** restrict the global headlines, ECDFs or panel D.
- `show_plot = interactive()`: PNG/PDF are saved, then the RStudio preview is
  displayed. Use Zoom or open the PDF to read the complete figure.

Required files are the selected iteration's completed manifest, markers, priors,
raw `results/*.rds` files, and matching `summary/res_summary.RData`,
`summary/res_cs_summary.RData`, `summary/em_summary_metadata.rds`.
If summaries are missing or stale, run `generate_summary_results_em.R` first
with the same iteration. The CS summary must include `component`, `cs_snps`,
`cs_size_snps` and lead fields.

The script reads each cohort gene's EM file once for panel D and prints progress
every 100 genes. Saved CS counts, components, membership, lead identity/coding and
PIP must agree with the summaries. Malformed probabilities, missing selected
fits or mismatches stop the run; there is no fallback to older or baseline
weighted fits. Completion checks do not independently establish EM convergence.

## Panels and denominators

**A — Four CS-count tables, plus overall agreement.** Keep the four tissues with
the highest CS-count disagreement rate among tissues meeting the table threshold.
Ties prefer larger region counts, then tissue name. Rows are additive SuSiE CS
counts; columns are EM counts; cells are exact numbers of gene-tissue regions.
The 0/0 cell is outside the cohort. The headline is:

```text
number of cohort regions with n_add == n_em / all cohort regions
```

It is not an average of the four displayed tissue rates.

**B — Lead-distance ECDFs, plus overall lead retention.** The headline is the
fraction of all matched CS pairs with exactly the same biological variant ID,
ignoring coding. Matching first reserves identical leads, then minimizes total
one-to-one physical distance among remaining CSs in each region.

The three curves retain the original definitions:
1. Changed leads in matched pairs.
2. Extra EM CSs to the nearest SuSiE lead.
3. Extra SuSiE CSs to the nearest EM lead.

Nearest-lead references are not additional one-to-one matches. Model-only
regions have no other-model lead, so their extra CSs are excluded from these
curves and counted in the note/audit. Identical matched leads enter the headline,
but not the changed-pair ECDF. Different alleles at the same coordinate remain
different leads with zero distance. The `log10(1 + bp)` axis retains zero.

**C — TSS-distance ECDFs for exactly the extra-CS groups in B.** Plot absolute
lead-to-TSS distances, excluding missing values, and show both the valid TSS
sample size and corresponding B denominator. Model-only regions remain excluded.
The two groups are not paired observations, and TSS proximity does not establish
causal accuracy. This preserves the preferred figure's C population.

**D — Posterior coding support in the same cohort.** For each EM component
associated with a reported credible set, calculate:

```text
q(component, coding) = sum of alpha[component, predictor]
                      across ALL predictors with that coding
```

Use the full component posterior, including mass outside the reported credible
set. Additive, recessive and dominant probabilities sum to one (within numerical
tolerance). Do not replace them with lead-coding counts, sums of PIPs, or the
learned coding priors; do not truncate/renormalize to CS members. Components
without reported CSs are excluded.

Average these probabilities with **equal weight per reported EM CS** within each
tissue. The plot shows mean recessive and dominant support; additive support
remains in the denominator and is reported in the pooled headline and tables.
The pooled mean is across all selected EM CSs, not an unweighted tissue average.
Tissues are ordered by mean recessive + dominant support, with their EM CS counts
at the right. The number of CSs is a descriptive denominator, not a claim of
independent replications.

All cohort regions remain in the audit. Regions with no EM CS have no component
to average; they are not assigned zero nonadditive probability. A tissue with no
EM CS remains displayed with undefined support and an explicit label. Panel D is
posterior support **among reported CS components**, not the fraction of regions
that truly have nonadditive effects.

## Outputs

Each run creates a separate directory:

```text
results_em/iteration_NNN/descriptive_results/figure1/run_YYYYMMDD_HHMMSS_PID/
```

The console prints the directory, also available as:

```r
figure1_result$output_dir
```

Main outputs:

- `figure1_em.pdf` (vector, 16 x 12.5 inches) and `figure1_em.png` (4000 x 3125).
- `figure1_caption.txt`, `figure1_plot_data.rds`, `sessionInfo.txt`.
- `metrics.csv`: headline numerators, denominators and rates.
- `cohort.csv`, `region_audit.csv`: included regions and all selection reasons.
- `tissue_ranking.csv`, `selected_tissues.csv`, `tissue_cs_count_tables.csv`
  (the four tables), `tissue_count_table.csv` (counts across all cohort tissues).
- `matched_cs_pairs.csv`, `unpaired_cs.csv`, `lead_retention_summary.csv`,
  `lead_distance_summary.csv`, `tss_distance_summary.csv`, `orphan_cs_rows.csv`.
- `coding_cs_support.csv`: each reported EM CS's full additive/recessive/dominant
  posterior support, component index and posterior mass inside the CS.
- `coding_region_audit.csv`: validates how many reported components each cohort
  region contributes, including zero-CS regions.
- `coding_tissue_support.csv`, `coding_overall_support.csv`: the exact D data.
- `priors_used.csv`: provenance only; these are **not** the D values.

For a fast redraw, set `plot_data_file` to the exported
`figure1_plot_data.rds` in the entry settings and Source again. This uses the
export's frozen iteration and cohort, without raw fits. A conflicting explicit
iteration is rejected.

The previous CYP4F11/scatter version's `figure1_em_v1` bundles do not contain the
required posterior coding data and are intentionally rejected. Run once on RCC
with `plot_data_file = NULL` to generate the new v2 bundle. Existing exported
files are not modified by the default new-run output path.

## Validation

```text
Rscript --vanilla script/analysis/tests/test_figure1_em.R
```

The test uses synthetic saved EM fits and summaries, with 27 tissues and a newer
unfinished iteration. It runs the actual RStudio Source entry, checks the shared
cohort and both headline denominators, compares A-C data directly with the
original table/ECDF implementation, verifies full-alpha coding probabilities
(including mass outside the CS and excluding non-reported components), rejects
stale/mismatched inputs, and tests zero-CS tissues, empty ECDFs and frozen redraws.
It renders PNG/PDF and a small-device preview without GTEx inputs or extra
packages. Previews under `tmp/figure1_em_v2_validation/` are explicitly synthetic.
Real results still require the EM files on RCC.

