# EM screening for downstream GWAS colocalization

Run in RStudio by sourcing `script/analysis/screen_em_coloc_candidates.R`.
The settings list at its top controls the project, GTEx directory, iteration,
P-value threshold, and distance threshold. On RCC:

```r
project_dir <- "/project2/mstephens/wdenault/susie_mix"
source(file.path(project_dir, "script/analysis/screen_em_coloc_candidates.R"))
em_coloc_screen$output_dir
```

Or from the project root:

```bash
Rscript --vanilla script/analysis/screen_em_coloc_candidates.R /project2/mstephens/wdenault/susie_mix latest /project2/mstephens/gtex
```

Copy both `screen_em_coloc_candidates.R` and `em_coloc_screen_utils.R` to RCC.
They use the existing `em_analysis_utils.R`, `generate_summary_results.R`,
`fit_mix_expression_plot_utils.R`, and GTEx annotation helper in this repository.
No EM refitting, new iteration, summary generation, GWAS download, or actual
colocalization is performed. The screen reads saved fits directly; an iteration
does not need pre-generated summaries. The latest complete iteration is fixed
once at startup using the existing manifest/completion/priors checks. A prepared
newer iteration is skipped with a warning. Completed iterations can still contain
saved gene/tissue errors; these are audited.

## Selection and diagnostics

1. Require saved gene-tissue `min_pv < 5e-8`, strictly. This is the minimum
   **marginal additive** cis-SNP P value saved by the baseline workhorse. It is
   neither a CS-specific P value nor evidence for nonadditivity. This initial
   cohort can miss pure nonadditive associations without an additive signal.
2. Check convergence, valid predictor maps, CS component indices, and the
   existing EM/baseline provenance checks. No extra sample-size or read-count
   eligibility filter is applied. Genotype reconstruction reuses the original
   variant QC (MAF/HWE/predictor support) and verifies the saved sample count.
3. Keep **all EM CSs**, including those with an additive lead. A CS lead is the
   highest component-alpha predictor within the saved CS, as in the existing
   one-CS scripts. `em_nonadditive_mass` sums the component's posterior mass on
   recessive and dominant predictors across the region; a lead label alone is
   not a posterior probability that the effect is nonadditive.
4. Compare each EM lead with **every vanilla CS lead**. The minimum of these
   distances must be strictly greater than 10,000 bp. Exactly 10,000 bp fails.
   Missing or incomparable coordinates fail this screen. A vanilla fit with no
   CS belongs to `no_vanilla_cs.csv`; it does not automatically pass the shift
   test. There is no one-CS restriction on either fit.
5. For every pair of EM and vanilla CSs, report overlap in biological SNPs,
   Jaccard overlap, and the fraction of each CS that overlaps. Coding duplicates
   count once. For gene-tissue pairs with a shifted CS, compute every lead-pair
   LD r-squared from raw additive genotypes in the reconstructed tissue samples.
   Also report squared correlation between the EM lead's actual coding and the
   vanilla additive lead. The latter is a predictor correlation, not standard
   LD. Record genotype class counts for EM leads. No LD/overlap cutoff is imposed.

Distance is a way of finding illustrative cases, not a test for distinct causal
signals. A 100-kb shift can have r-squared = 1. Rank cases using distance, overlap,
LD, posterior uncertainty, and genotype support together. LD here describes the
eQTL sample; it is not the ancestry-matched LD required to fine-map a GWAS.

The entry script also defaults to `unshifted_control_fraction = .01`: a reproducible
1% sample of otherwise eligible gene–tissue pairs whose EM leads are all unshifted.
The hash-based selection is independent of external associations/H4. Controls get
BF and LD exports, with `selection_group = "unshifted_control"` in the BF manifest,
and do not enter `candidates.csv`. Set 0 to disable or 1 to export every eligible
unshifted pair. Compare groups separately; these are not power/MAF-matched controls.

The default LD reader needs GTEx files, PLINK2, and the packages used by the
existing expression plotting helper. `compute_ld = FALSE` explicitly skips LD.
Missing genotypes/packages or reconstruction failures leave LD as NA with an
error/disabled status; they never become zero LD. For another genotype source,
call `em_screen_coloc(..., genotype_reader = function(gene, tissue, x) ...)` with
a complete 0/1/2 matrix from the fitted samples, named with saved SNP IDs.

## Outputs

Each run creates a separate directory under
`results_em/iteration_NNN/coloc_screen/screen_TIMESTAMP_ID/`. The path is printed
and returned in `em_coloc_screen$output_dir`. An optional explicit `output_dir`
must be empty. This prevents stale exports from a previous screen.

| File | Content |
| --- | --- |
| `all_em_cs.csv` | Every CS passing gene-tissue P/convergence/map checks, with shift classification |
| `candidates.csv` | Strict shifts from all vanilla leads; consult BF and LD status columns |
| `no_vanilla_cs.csv` | Significant gene-tissue pairs with EM CSs and no vanilla CS |
| `cs_pair_comparisons.csv` | Every EM-CS/vanilla-CS distance, overlap, and available LD |
| `bf_manifest.csv` | Gene/tissue to export path, status, common-SNP count and retained posterior mass |
| `bf/case_NNNNNN.rds` | Native and common-SNP BF matrices, variant metadata, priors and CS definitions |
| `bf/case_NNNNNN.tsv.gz` | Full-region SNP-level natural-log BFs for every exported CS (optional) |
| `audit.csv` | Gene/tissue errors, invalid inputs, missing BFs and LD failures |
| `screen_counts.csv` | Counts through the screen, including failures/no-CS outcomes |
| `screen_metadata.rds`, `em_priors_used.csv` | Iteration, settings, priors, definitions and R session information |

BF export is attempted for shifted gene-tissue pairs and the separate no-vanilla-CS
group. All CS components of both methods are exported for each such pair, not
only the shifted component. All SNPs with positive prior support in each fitted
region are included, not just CS members. Failure to export BFs does not erase a
shift: it remains in the candidate table with `bf_status = "error"`.

## BF definition and coloc handoff

For component l, biological SNP j and coding c, let w[j,c] be the saved prior
probability of that predictor and q[j] = sum_c w[j,c]. Export

```
log BF[l,j] = logsumexp_c(log BF_predictor[l,j,c] + log(w[j,c] / q[j])).
```

This marginalizes coding uncertainty using **prior** weights. It does not take
the best coding BF, average log BFs, add posterior weights, or multiply q[j]
into the BF. The workhorse divides each coding's mass among its retained
predictors; if coverage differs by coding, q[j] need not be uniform. Therefore
using `pi_add/pi_rec/pi_dom` alone is not generally correct. The exporter uses
`weighted_mix_prior_weights`, checks it against `fit$pi` when present, and
validates that the saved alpha is proportional to BF times predictor prior.
Missing/inconsistent BFs are audited rather than inferred from PIPs. Zero-prior
SNPs are explicitly recorded and excluded from supported matrices.
As in SuSiE colocalization, these are component-level approximate BFs conditional
on the fitted other effects. Learned EM coding priors are treated as plug-in
estimates; uncertainty in those priors is not integrated by this exporter.

`bundle$em$lbf` and `bundle$vanilla$lbf` have rows in `sets$cs_index` order and
columns named `CHR:POS:REF:ALT` on GRCh38. These are natural logarithms, ready for
the log-BF interface of [coloc.bf_bf](https://chr1swallace.github.io/coloc/reference/coloc.bf_bf.html).
The API does not require the GWAS and eQTL to use the same genetic coding.
The biological sharing hypothesis is at the variant level. It does require
meaningful comparable evidence under each model and correctly identified variants.

`variants` preserves source SNP IDs, REF/ALT, build, and the PLINK counted allele.
The parser accepts GTEx v8 IDs; unknown IDs are audited for BF export, not guessed
from position. `snp_prior` is separate. `alpha` is coding-collapsed component
posterior mass. `all_component_lbf` also preserves components without CSs for
diagnostics; do not automatically treat these as discovered signals. No vanilla
CS means a zero-row vanilla CS BF matrix, not evidence for H0.

Before a paired coloc analysis:

- Obtain full regional GWAS evidence and signal-specific log BFs, generally by
  GWAS SuSiE fine-mapping with appropriate LD. Catalog lead associations alone
  are insufficient ([Catalog summary-statistics documentation](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/)).
- Harmonize genome build and variant alleles, and use the **same three-way SNP
  intersection** for EM, vanilla, and GWAS. The bundle's `paired` matrices only
  apply the first intersection, between the two eQTL methods.
- Recompute retained posterior mass in the full native signals after the final
  intersection. Large losses make a comparison unreliable. Do this before
  subsetting; coloc cannot detect mass already discarded from its input.
- Set the same p1, p2, p12 and the same SNP-location prior policy for both method
  comparisons. A uniform per-SNP coloc prior is a clear primary choice. The
  saved model-specific `snp_prior` can be a sensitivity analysis using coloc's
  variant prior-weight arguments in versions that support them. Declare this
  change explicitly: it tests spatial-prior changes as well as coding changes.
- Compare every relevant signal against the same GWAS signals and report losses
  as well as gains. Component numbers need not match across methods.

Minimal handoff, after GWAS harmonization and coverage checks:

```r
b <- readRDS(file.path(em_coloc_screen$output_dir,
                      subset(em_coloc_screen$bf_manifest,
                             bf_status == "exported")$bundle[1]))
# gwas_lbf: matrix of natural-log BFs, rows = GWAS signals,
# columns = verified GRCh38 CHR:POS:REF:ALT identities.
common <- Reduce(intersect, list(colnames(b$em$lbf),
                                colnames(b$vanilla$lbf), colnames(gwas_lbf)))
stopifnot(length(common) > 1L, nrow(b$vanilla$lbf) > 0L)
em_retained <- rowSums(b$em$alpha[, common, drop = FALSE])
vanilla_retained <- rowSums(b$vanilla$alpha[, common, drop = FALSE])
# Also check GWAS retained mass against its native signal posterior.
# Choose and justify a common coverage rule before looking at coloc changes.
coloc_priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6) # illustrative; assess sensitivity
em_coloc <- do.call(coloc::coloc.bf_bf, c(list(
  bf1 = b$em$lbf[, common, drop = FALSE],
  bf2 = gwas_lbf[, common, drop = FALSE]), coloc_priors))
vanilla_coloc <- do.call(coloc::coloc.bf_bf, c(list(
  bf1 = b$vanilla$lbf[, common, drop = FALSE],
  bf2 = gwas_lbf[, common, drop = FALSE]), coloc_priors))
```

The BF export was tested against installed `susieR`; actual GWAS coloc is a
separate stage. This implementation only handles EM mix. The same variant-level
interface can later accept slider BFs that properly marginalize its coding prior.

## Trait panel

`coloc_trait_panel.csv` contains the original 32 disease, growth and physiological
traits plus 18 biomarkers, with source links and analysis notes. Exact releases,
schemas and checksums are pinned in `coloc_gwas_sources.json` and
`coloc_biomarker_sources.json`; these registries govern downloads. The panel is
not a claim that every row is the newest/largest accessible study. Public headers
were checked on 2026-10-11/12; complete RCC downloads have not been run locally.
Prefer dense, compatible summary statistics to a larger study with inadequate
regional coverage. Gene-specific protein/eQTL sources are selected separately
when preparing a new job, before examining external associations.

The panel includes AD; schizophrenia, bipolar disorder and depression; seven
early-growth traits; adult height, BMI and waist-to-hip ratio; four lipid traits;
four glycemic traits; T2D and CAD; three kidney traits; IBD and its two main
subtypes; breast and prostate cancer. Priority 1 is a suggested starting subset,
not a statistical weight. Correlated traits and repeated tissues are not
independent replications. Look up GWAS by the screened genomic regions; requiring
a prior gene-to-trait annotation would bias the analysis against new links.

For AD, the [2026 consortium consensus study](https://www.nature.com/articles/s41588-026-02583-1)
provides main, no-proxy and no-biobank summary releases (GCST90704646/7/8).
Use the latter two as phenotype sensitivity analyses, not additional traits.
For growth, [EGG](https://egg-consortium.org/) supplies dedicated releases; adult
GTEx findings are mechanistic follow-up hypotheses rather than demonstrations
of an effect in developing tissue. The 2024 pubertal-growth download is
trans-ancestry, so single-population reference LD is not automatically suitable.

## Interpreting changed colocalization

With expression as dataset 1 and GWAS as dataset 2, H0 means neither is associated,
H1 expression only, H2 GWAS only, H3 both associated at different variants, and
H4 shared variant. High H3 changing to high H4 is a particularly informative
change in localization. H0/H1/H2 or an unresolved posterior changing to H4 is
more naturally described as evidence gained. H0 is not the generic label for
failure to colocalize. At a strongly associated eQTL and GWAS, inspect a large
H0 carefully for signal strength, missing SNPs, scaling and prior effects.

Report all five posterior probabilities, delta PP.H4, and PP.H4/(PP.H3+PP.H4)
when the denominator is substantial. Report the posterior distribution of the
shared SNP under H4 to establish whether *where* the traits colocalize changes.
Neither a high conditional H4 ratio with little evidence for association nor
an arbitrary H4 threshold crossing alone is strong evidence of improvement.
See [coloc's SuSiE method](https://pmc.ncbi.nlm.nih.gov/articles/PMC8504726/) and
[prior sensitivity guidance](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html).

Multiple robust changes at independent loci would support practical value.
The shifted subset is deliberately selected for changed localization, so it
cannot estimate the overall improvement rate. A later benchmark should use a
fixed trait panel and a broader eQTL cohort, count gains and losses, and calibrate
the entire selection/analysis procedure in simulations with additive truth,
nonadditive truth, shared variants and distinct variants. Check prior sensitivity,
SNP/LD coverage and independent replication. Colocalization supports a shared
causal variant; it does not on its own establish expression-mediated causality.

## Validation

The next stage is implemented in [README_em_gwas_coloc.md](README_em_gwas_coloc.md):
pinned downloads for the original panel plus 18 biomarkers, gene-matched protein
and external eQTL comparisons, and paired EM/vanilla H0–H4 comparisons,
with gains, losses, shared-SNP locations and explicit unassessed cases.

```bash
Rscript --vanilla script/analysis/tests/test_em_coloc_screen.R
```

Fixtures cover the latest completed iteration, strict P and 10-kb boundaries,
all vanilla leads, CS overlap, perfect LD despite a shift, missing coordinates,
no vanilla/EM CSs, failed BF/LD exports, unequal coding coverage, evidence-preserving
BF marginalization, empty outputs and the RStudio Source entry. If `susieR` is
installed, a two-signal fitted-data smoke test also checks the exporter.
