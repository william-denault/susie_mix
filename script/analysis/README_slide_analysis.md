# Analyze the saved SuSiE-slide fits against additive SuSiE

On RCC, run these from `/project2/mstephens/wdenault/susie_mix`:

```bash
Rscript --vanilla script/analysis/generate_summary_results_slide.R
Rscript --vanilla script/analysis/descriptive_results_slide.R
Rscript --vanilla script/analysis/finding_interesting_dominant_1cs_slide.R
Rscript --vanilla script/analysis/finding_interesting_recessive_1cs_slide.R
```

Each script accepts the project directory as its first argument, or reads
`SUSIE_MIX_PROJECT_DIR`. RStudio **Source** also works: edit the settings list
at the top. No fitting, EM iteration, PLINK, or expression data is needed.

## Inputs and outputs

Inputs are `results/GENE.rds`, as written by the current `workhorse.R`:
`susie_add`, `susie_add_perm`, `fit_slide`, `fit_slide_perm`, and
`add_predictor_map`. Slide fits use the additive genotype matrix, not the
expanded mixed-coding map. Existing mixed-model fits are not substituted for
missing slide fits. Missing/invalid fits and saved workhorse failures appear
in `res_errors.csv`; an absent permutation is NA, not a zero-CS result.

Outputs are isolated under `results_slide/`:

- `summary/`: `res_summary.RData`, `res_cs_summary.RData`, CSV copies,
  `res_errors.csv`, and `slide_summary_metadata.rds`. The RData object names
  remain `res_summary` and `res_cs_summary`. CS rows retain component indices,
  lead SNP/PIP/delta, count-forced status, coverage, purity, SNP membership,
  and TSS-distance provenance. Tissue rows include real/permutation CS counts,
  coding counts, convergence, biological SNP overlap, lead agreement, and
  ELBO differences.
- `descriptive_results/`: overall/tissue counts, coding categories, lead
  agreement (all pairs and one CS in both), CS-count comparisons, paired
  permutation summaries, primary pairs, CS deltas, TSS audits/distributions,
  `coding_and_tss_summary_2x2.pdf` (also PNG), and
  `tissue_specific_coding_and_tss_4x2.pdf`.
- `plot/one_cs_dominant/` and `plot/one_cs_recessive/`:
  `lead_snp_comparisons.csv`, `lead_distance_summary.csv`,
  `largest_lead_shifts.csv`, and `lead_snp_distance_overview.png`.

The descriptive filter remains P < 1e-8 and mean reads >= 100; the one-CS
filter remains P < 5e-8 and mean reads >= 100. Convergence flags are retained;
they are not an additional filter. Regenerate summaries after replacing RDS
files. ELBO differences describe the fitted objectives, not a calibrated test.

## Continuous coding and one-CS comparisons

Slide uses `x + delta * I(x == 1)`. Delta = -1 is recessive, 0 additive,
and +1 dominant for the workhorse's minor-allele-oriented genotype. Delta is
component-by-SNP: the code uses `sets$cs_index`, not the CS list position,
to select its row. The lead is the largest component alpha within that CS,
matching the existing agreement summaries. SNP PIP is not a coding-specific
posterior probability.

Summaries preserve the exact delta and distinguish additive, partial
recessive, recessive, partial dominant, and dominant, using a numerical
tolerance of 1e-6 at -1, 0, and +1. Count-forced additive leads are recorded.
The one-CS scripts default to endpoint dominant/recessive leads. To include
partial effects, use the second argument `direction`:

```bash
Rscript --vanilla script/analysis/finding_interesting_dominant_1cs_slide.R /project2/mstephens/wdenault/susie_mix direction
Rscript --vanilla script/analysis/finding_interesting_recessive_1cs_slide.R /project2/mstephens/wdenault/susie_mix direction
```

The comparison table retains `slide_coding_class` and `lead_delta` in either
mode. Set `both_one_cs = TRUE` in the script for strictly one CS in both
models. Otherwise, the additive comparator is its highest-PIP CS lead; when
it has no CS, it is the highest-PIP SNP, explicitly labeled in the table and
plot. All eligible pairs, including equal leads, remain in the CSV. The plot
ranks changed leads and shows the 20 largest shifts beside distance versus
slide SNP PIP. Physical distance alone does not establish low LD.

## TSS disagreement rule for slide, original, and EM results

Agreement means the **same biological CS lead SNP in the same gene/tissue**,
ignoring coding. Every shared lead is excluded individually, even in a region
with other disagreeing CSs. Different leads remain eligible even if CS members
overlap. A CS found only by one model remains eligible. Unknown leads are
audited and excluded. Coding and general summaries retain the full primary set.

The original `descriptive results.R` and weighted/EM descriptive scripts already
use this rule via `tss_disagreement_utils.R`. Rerun those scripts to replace
previous all-CS plots. The slide analysis uses the same helper. Each model's
TSS curve is normalized over its selected finite distances inside the existing
window (10-kb bins centered at -200 through +200 kb). The agreement audit and
selection-count CSVs report exclusions, missing distances, and out-of-window CSs.

The current workhorse does not save slide-to-TSS distances. The slide summary
uses explicit saved gene metadata first, then the same GTF/first gene row as the
workhorse when available. It can recover TSS as position minus the named signed
offsets saved by the current `workhorse_utils.R`, requiring all usable estimates
to agree. This fallback assumes signed genomic offsets, not legacy absolute
distances. Without a usable TSS, distance is NA and is audited. GTF-derived
distances are oriented along transcription; recovered offsets use genomic
orientation when strand is unavailable. Every CS records its orientation/source.

All slide analysis and plots use base R. Loading the RCC GTF additionally needs
`data.table`. The slide package itself is not needed to read saved fits.

## Validation

```bash
Rscript --vanilla script/analysis/tests/test_slide_analysis.R
Rscript --vanilla script/analysis/tests/test_tss_disagreement.R
Rscript --vanilla script/analysis/tests/test_em_analysis.R
Rscript --vanilla script/analysis/tests/test_weighted_summary.R
```

The synthetic slide test covers reordered components, endpoint/partial coding,
permutations, missing fits, zero CSs, signed TSS recovery, shared/partly shared
leads, distances, no-CS additive fallback, plots, and empty selections.
