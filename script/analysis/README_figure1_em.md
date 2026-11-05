# Figure 1: final EM SuSiE-mix versus additive SuSiE

Open `script/analysis/plot_figure1_em.R` in RStudio on RCC and click **Source**.
The script uses the existing project helpers; copy/sync the new entry script and
`figure1_em_utils.R` into the project's `script/analysis/` directory together.
It does not run SuSiE or restart EM.

```r
source("/project2/mstephens/wdenault/susie_mix/script/analysis/plot_figure1_em.R")
```

The settings at the top default to project
`/project2/mstephens/wdenault/susie_mix`, GTEx directory
`/project2/mstephens/gtex`, and `iteration = "latest"`. If the project is elsewhere,
edit `project_dir`, or set `SUSIE_MIX_PROJECT_DIR`. Set `iteration = 18L` for a
permanently pinned iteration-18 figure.

“Latest” is resolved **once**, using `prior_history.csv`, the iteration manifest,
completion markers, result files and coding priors. An incomplete newer iteration
is skipped with a warning. Completion is not itself a claim of EM convergence.
All panels use that one iteration's `weighted_fit_mix` fits. Old baseline weighted
fits, unweighted SuSiE-mix and SuSiE-slide are not used as EM results.

## Required inputs

- Matching EM summaries under `results_em/iteration_NNN/summary/`, including
  `em_summary_metadata.rds`, `res_summary.RData`, and `res_cs_summary.RData`.
  If missing or stale, run `script/analysis/generate_summary_results_em.R` first,
  with the same iteration setting. CS summaries must include `cs_snps` and
  `cs_size_snps`; regenerate summaries if using an older summary schema.
- Original baseline fits under `results/`, and selected EM fits under
  `results_em/iteration_NNN/results/`. CYP4F11's saved fits and summary are checked
  for matching CS counts, lead identities/codings, memberships and PIPs.
- The GTEx v8 subject/sample annotations, gene-read GCT, GRCh38 gene GTF, genotype
  PLINK dataset and `plink2` executable already used by the existing expression
  scripts. Filenames are defined in `fig1_reconstruct_example()`.
- R packages **data.table** and **matrixStats** for the existing GTEx reconstruction
  helper. Graphics, tables and replotting use base R; no ggplot2 or susieR required.

Use an RCC RStudio session with sufficient memory for the original GTEx expression
reconstruction pipeline. This reads the gene-read matrix to reconstruct library
sizes, even though panel D displays just one gene. Temporary PLINK extraction uses
R's `tempdir()` and the existing unique-prefix cleanup helper.

## One cohort, four panels

The parent cohort is gene-tissue regions with **minimum marginal additive
P < 1e-8, mean reads >= 100, and at least one CS in either model**. Missing EM fits
and missing filter values are excluded. No tissue-size cutoff is imposed. These
are gene-tissue regions, not unique genes or independent replications.

| Panel | Display and denominator |
| --- | --- |
| A | Pooled and tissue-level region outcome bars. Same count/same biological lead multiset; same count/changed leads; fewer EM CSs with both models positive; SuSiE only; more EM CSs with both positive; EM only; unavailable lead identity. The colors are mutually exclusive, and all regions contribute. |
| B | Same-lead retention among all matched CS pairs, then lead separation versus coding-specific EM lead PIP for changed pairs. Color = EM coding; shape = CS membership overlap. Invalid PIPs are audited and omitted from the scatter only. |
| C | SuSiE versus EM absolute lead-to-TSS distance for the same changed-pair population, with both usable distances. Below the diagonal = closer in EM. No PIP threshold; missing, inconsistent or stale TSS records are audited. |
| D | CYP4F11 / Adrenal Gland, checked to be in the same parent cohort. Saved regional predictor PIPs and genotype–expression plot at its recessive EM lead, preferring a lead shared with additive SuSiE. |

B and C are defined subsets of the parent cohort, rather than new region cohorts.
Their valid-data denominators can differ and are printed. Matching first reserves
exact biological lead IDs (ignoring coding), then minimizes total one-to-one
physical distance for the remaining CSs. Extra CSs are exported, not forced into
pairs. Different alleles at the same coordinate remain different leads with zero
physical separation. The log1p axes retain zero. CS overlap uses unique biological
SNP membership and Jaccard, also ignoring coding.

D must pass the same filters; the script stops rather than relaxing selection if
CYP4F11 / Adrenal Gland no longer qualifies or has no recessive EM lead. Its
phenotype and genotypes are reconstructed with the existing GTEx pipeline and
aligned to the saved predictors and sample count. Conditional expression subtracts
the other EM CS leads' fitted contributions from a joint OLS regression. The
additive/recessive lines are descriptive OLS curves on that adjusted expression,
**not fitted SuSiE posterior predictions or independent validation**. The saved
tables include the unadjusted transformed phenotype, group counts, means and
conditioning coefficients (in the RDS bundle).

Lead identity, CS membership, coding and CS counts answer different questions.
Neither a physical lead shift nor TSS proximity alone establishes improved causal
localization. Coding-specific lead PIP is not a CS-existence probability. A selected
example does not establish genome-wide nonadditivity. Simulations can provide the
calibration and ground-truth assessment that the real-data comparisons cannot.

## Outputs and redraws

By default each run creates:

```text
results_em/iteration_NNN/descriptive_results/figure1/run_YYYYMMDD_HHMMSS_PID/
    figure1_em.pdf                  # Vector figure, 15 x 12 inches
    figure1_em.png                  # 3750 x 3000 pixels
    figure1_caption.txt             # Definitions and calculated denominators
    figure1_plot_data.rds           # Frozen plotting data and provenance
    metrics.csv
    overall.csv
    tissue_outcomes.csv
    tissue_count_table.csv          # Exact CS-count contingency counts
    cohort.csv                     # Included regions and outcomes
    region_audit.csv               # Included/excluded region reasons
    matched_cs_pairs.csv            # Matching, PIPs, Jaccard, TSS and panel flags
    unpaired_cs.csv                 # Extra CSs; nearest-lead distances and TSS
    membership_audit.csv
    priors.csv
    example_leads.csv
    example_regional.csv
    example_groups.csv
    example_expression.csv         # Anonymous individual plotting values
    sessionInfo.txt
```

The path is printed in the console and available as `figure1_result$output_dir`.
The PNG/PDF are written before the optional RStudio Plots display. Use **Zoom**
or open the PDF for the full layout. If the Plots pane is too small, set
`show_plot = FALSE` and view the saved files.

For a fast redraw without GTEx data or installed reconstruction packages, set
`plot_data_file` to the exported `figure1_plot_data.rds` in the entry settings and
Source again. Replotting keeps the export's iteration, cohort, CYP4F11 data and
any synthetic-data label. A conflicting explicit iteration is rejected. The
`latest` setting in replot mode means the frozen export, not a new scan of results.

## Local validation

```text
Rscript --vanilla script/analysis/tests/test_figure1_em.R
```

This base-R test builds synthetic summaries and saved fits, exercises the actual EM
discovery/summary/join pipeline (including a newer incomplete iteration), verifies
shared-cohort counts, matching, membership, TSS exclusions and expression adjustment,
and renders all four panels. Only external GTEx reconstruction is replaced. It also
checks the Source entry's replot route and the zero-changed-leads case. Synthetic
previews in `tmp/figure1_em_validation/` are explicitly labelled and are not GTEx
results. Real-data execution still needs the inputs on RCC.
