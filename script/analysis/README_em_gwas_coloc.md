# GWAS follow-up: EM mix versus vanilla

This follows `screen_em_coloc_candidates.R`. It freezes the latest **completed**
EM iteration and its most recent completed screen, downloads a fixed panel of 32
full-summary GWAS releases, and compares the saved expression log BFs with GWAS
log BFs. It does not refit the expression models or change screening thresholds.
The screened gene–tissue P < 5e-8 is the baseline marginal additive association
test, not a test for a specific CS or for nonadditivity.

## Run on RCC

Requirements: Python >=3.9 (standard library only), R packages `data.table`,
`jsonlite`, `digest`, `R.utils`, `coloc` (tested with 5.2.3), and `susieR` for the
GWAS multi-signal analysis. PLINK 2 is required when using reference genotypes.
Full GWAS downloads require substantial disk space; use a cluster filesystem.

1. Source `script/analysis/download_em_coloc_gwas.R`. Defaults select the latest
   completed iteration and latest completed screening run for that iteration.
   Run the screening script first if there is no such screen. `traits = NULL`
   requests all 32; `c("ad", "ldl")` is a useful initial subset. No GWAS P-value
   threshold is applied to the extracted regions.
2. The script prints a **frozen job directory** and sets `em_gwas_job_dir` in
   RStudio. To resume a download, set `job_dir` to that path in the download script.
   Completed cached files are checksum-verified and reused; interrupted individual
   downloads restart. Do not start simultaneous download processes for the same
   job/cache. Editing a source definition requires a new job.
3. Source `script/analysis/run_em_gwas_coloc.R` with that `job_dir`. Its default
   `modes = c("abf", "susie")` runs the explicitly labelled single-signal GWAS
   sensitivity analysis and records `pending_ld` for GWAS SuSiE until LD is set.
   This is not a substitution of ABF results for the primary analysis.
4. After identifying the LD reference, set `ld_bfile`, `ld_ancestry = "EUR"`,
   `ld_build = "GRCh38"`, `plink`, and (where needed) `ld_keep_file`, then rerun.
   Each analysis gets a new output directory; downloads and expression BFs stay
   frozen even if EM subsequently advances.

```bash
Rscript script/analysis/download_em_coloc_gwas.R /project2/mstephens/wdenault/susie_mix latest
Rscript script/analysis/run_em_gwas_coloc.R /project2/mstephens/wdenault/susie_mix /path/to/frozen/job
# Resume/extract only selected traits directly:
python3 script/analysis/download_coloc_gwas.py --job-dir /path/to/frozen/job \
  --sources /path/to/frozen/job/sources.json --cache /path/to/gwas_cache --traits ad,ldl
```

GTEx access alone does not establish a suitable LD reference. Check build,
ancestry, unrelatedness and donor count. A PLINK keep file has FID and IID columns.
All donors in the prefix are used if no keep file is supplied; the prefix must
then already contain a suitable subset. GRCh37 reference genotypes must be
properly harmonized before supplying them here. Prefer a large matched or
in-study GWAS reference when available: a small external GTEx subset can have
substantial LD estimation error. No individual genotypes are uploaded.

## Sources and deliberate exceptions

`coloc_gwas_sources.json` pins URLs, schemas, builds, ancestry, published checksums
where available, and effect-scale assumptions. Headers and initial records were
checked on 2026-10-11; tests include all 32 schemas. The job copies this registry.
The downloader retains original files and records SHA256, MD5 and chain-file
provenance. A checksum-less consortium URL is frozen by its first cached content;
it is not guaranteed to be historically immutable before the first download.

| Traits | Pinned release |
|---|---|
| Alzheimer disease | EADB/IGAP/PGC 2026, GCST90704646 |
| Schizophrenia; bipolar | PGC 2022 European autosomes; PGC 2021 all bipolar, European |
| Depression | PGC 2025 European, full results excluding 23andMe |
| Fetal birth weight; birth length | EGG 2019 maternal-adjusted fetal effect; EGG 2014 discovery |
| Childhood BMI; childhood obesity | EGG 2020 GCST90002409; European columns of EGG 2019 |
| Infant head circumference | EGG 2022, corrected file dated 2026-09-17, age 6–30 months |
| Pubertal growth; adiposity rebound age | EGG 2024 PG_II combined sex; EGG_AR_age |
| Adult height; BMI; WHRadjBMI | GIANT 2022 GCST90245992 public subset; Yengo 2018; Pulit 2019 |
| LDL; HDL; triglycerides; total cholesterol | GLGC 2021 European GCST90239658/652/664/676 |
| Fasting glucose; insulin; HbA1c; 2h glucose | MAGIC 2021 European GCST90002232/238/244/227 |
| T2D; CAD | Suzuki 2024 European; Aragam 2022 GCST90132314 |
| eGFR; UACR; CKD | CKDGen 2019 European GCST008059/8794/8065 |
| IBD; Crohn; ulcerative colitis | de Lange 2017 GCST004131/132/133 |
| Breast; prostate cancer | BCAC 2017 GCST004988; PRACTICAL 2018 GCST006085 |

The height file is a public subset, not the entire 5-million-person study. The
downloadable schizophrenia file also has a specific European sample definition.
These are not automatically replaced by newer/larger studies. The depression
top-10,000-hit download is deliberately excluded. Consult the pinned source
pages for acknowledgements and data terms. The T2D download form requires not
attempting participant identification; the consortium also prohibits reposting
its summary data to third-party websites. Cached data are local analysis inputs.

Two traits require special treatment:

- **Maternal-adjusted fetal birth weight:** beta and SE support ABF sensitivity,
  but `n_ownBW + n_offBW` is not a valid effective N for the adjusted estimate.
  GWAS RSS is withheld until a defensible sampling model/effective N is available.
- **Pubertal growth:** the public MR-MEGA file contains a multi-ancestry
  meta-regression. Its `beta_0/se_0` is not used as an ordinary GWAS Z score. The
  provided natural `lnBF` supports only a separately labelled single-signal
  sensitivity analysis with the source's own effect prior. A European LD matrix
  does not turn that analysis into an ancestry-matched SuSiE fit.

MAGIC's raw-trait beta/SE and inverse-normalized-trait P values have different
meanings. BFs use beta/SE together, never Z reconstructed from that P value.
For quantitative ABFs, a documented standardized scale uses sdY=1; otherwise
sdY is estimated from N, MAF and SE using coloc's variance relationship. This
regional estimate is approximate, especially for meta-analysis or reference
frequencies. It requires >=100 usable SNPs. Supply a study-specific raw-effect
`effect_prior_sd` when a reliable estimate is unavailable (often 2h glucose).
Case/control ABF prior SD is 0.2 on the log-odds scale; quantitative default is
0.15 times sdY. Each evaluated locus saves the actual effect prior used.

## Comparability and QC

- Download/extract the full expression SNP span plus 100 kb on either side;
  do not restrict GWAS to CS members, leads, significant SNPs or gene annotations.
- Convert source GRCh36/37 SNV coordinates to GRCh38 with UCSC chains; complement
  alleles on reverse mappings. Reject missing/ambiguous mapped bases. Parse
  harmonized Catalog `hm_*` columns together, including `hm_beta`.
- Match chromosome, position **and both alleles**. Retain autosomal biallelic
  SNVs only. Exclude strand-ambiguous A/T and C/G SNPs, all duplicate allele keys,
  and INFO <0.8 when supplied. There is no expression read-count/sample-size filter.
- Fit GWAS SuSiE once per trait/region using reference dosages oriented to the
  GWAS effect allele. Reference QC requires >=100 donors, MAF >0.01, missingness
  <=5%, frequency disagreement <=0.2 when GWAS EAF is supplied, >=80% regional
  SNP coverage, and retention of the strongest regional GWAS Z-score SNP.
- Use the median supplied per-SNP N, recording whether the source supplies total
  or effective N; withhold RSS if N varies by >20% or the `estimate_s_rss`
  diagnostic exceeds 0.1. These configurable cutoffs are conservative screening
  diagnostics, not guarantees of calibrated fine-mapping. Metadata-only N is an
  approximation. Require convergence and GWAS CSs; no CS is not an H0 estimate.
- Intersect EM, vanilla and GWAS SNPs once for the paired comparison. Require
  >=100 common SNPs and >=90% posterior support retained in **each tested
  component**, measured before restriction. QTL coverage uses both saved native
  alpha and the uniform-SNP-prior BF posterior; both must pass. A missing causal
  SNP should cause an unassessed comparison, not a spurious loss.
- Call `coloc::coloc.bf_bf` on **absolute natural-log BFs**, with uniform per-SNP
  coloc priors and the same GWAS fit for both methods. EM coding proportions have
  already been marginalized into expression BFs and are not multiplied in again.
  Defaults: p1=p2=1e-4 and p12=1e-6, 5e-6, 1e-5. Never renormalize BFs to sum one.

## Outputs

`all_signal_pairs.csv` contains every expression-CS × GWAS-signal × prior pair,
its status, retained mass, **H0, H1, H2, H3, H4**, H3+H4, and shared-SNP location.
The QTL component's lead coding and original EM shift classification are retained,
so a gain involving an unshifted companion CS can be distinguished from a gain
at the CS that triggered screening.
H1 means expression only; H2 means GWAS only. The conditional H4/(H3+H4) ratio is
reported only when H3+H4 >=0.1. Inspect the unconditional probabilities first.
Each case also has `snp_h4_posteriors.tsv.gz`, containing the full shared-SNP
distribution conditional on H4, and saved GWAS BFs/fits and QC metadata.

`paired_changes.csv` compares, for each fixed GWAS signal and prior, the largest
H4 pair in each expression method and retains the corresponding H0–H4 vectors,
their differences, component IDs, top-hypothesis transition, and H4-lead distance.
**A maximum across components is not a posterior probability of any shared
signal**, and component 1 in EM is not automatically matched to component 1 in
vanilla. All underlying pairs remain available. Unequal numbers of CSs can affect
this exploratory summary; assess that in calibration simulations.

A gain requires EM H4 >=0.8, vanilla H4 <0.8 and delta H4 >=0.1; losses use the
symmetrical rule. Thresholds are configurable. Neither direction is selected
away. If any component is unassessed, or vanilla has no CS, the method comparison
is unassessed rather than assigning zero H4 to the missing comparator.
`change_counts.csv` counts these labels separately by trait, mode and prior.
`audit.csv` accounts for every requested trait × case × mode, including failures,
missing data and pending LD. `completed.txt` means the traversal finished, not
that every requested analysis was assessable. Check the audit before interpretation.

## Recommended stress tests

Low observed eQTL–GWAS overlap can reflect weak association, imprecise fine-mapping,
the wrong tissue/cell state/age, or true biological differences. It is not itself
a calibrated negative control. See [Hukku et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33308443/).

1. **Known-truth simulations using GTEx genotypes.** Generate independent
   expression and phenotype samples with a shared causal SNP versus distinct
   causal SNPs spanning low to high LD. Vary additive/recessive/dominant effects,
   MAF, minor-homozygote count, N, effect strength, multiple signals and missing
   causal variants. Include additive truth and complete/one-trait nulls. Repeat
   the entire P/shift selection step; compare power at matched false-colocalization
   rates, posterior calibration and location coverage, not just numbers of gains.
2. **Independent donor splits of the same GTEx tissue.** Fit expression in one
   split and test colocalization with the same gene in another. Learn EM priors
   and choose candidates using training donors only. Repeated splits estimate
   stability, not independent replication counts; reduced N particularly hurts
   recessive signals. Independent external eQTL cohorts are stronger confirmation.
3. **GWAS–GWAS replication controls.** Use non-overlapping cohorts measuring the
   same phenotype to test the harmonization/LD/coloc workflow independently of
   expression model choice. Sample overlap and correlated traits are not independent
   replication. Cross-tissue QTL sharing is useful but is not guaranteed truth.
4. **Perturb the assumptions.** Repeat the prior grid (implemented), widen/narrow
   loci, remove the putative causal SNP, thin coverage, vary LD reference ancestry
   and sample size, and compare single- versus multi-signal analyses. Stable H3→H4
   changes at independent loci with adequate GWAS/QTL evidence are more persuasive
   than H0→H4 alone. Use negative controls simulated with realistic LD instead of
   permuting individual SNP Z scores, which destroys the LD structure.

Report results by independent genomic locus as well as gene–tissue–trait pair.
The >10-kb screen deliberately enriches for EM changes and misses other classes
of vanilla advantage; gains/losses here cannot estimate overall method performance.
Add an unselected or matched unshifted eQTL cohort for that claim. Growth traits
against adult GTEx also have an age/context mismatch. Colocalization supports a
shared variant, not expression-mediated causality or proof that a GWAS effect is
recessive. The same evaluation design can later be applied to slider BFs.

Method references: [coloc BF interface](https://chr1swallace.github.io/coloc/reference/coloc.bf_bf.html),
[SuSiE colocalization](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html),
[prior sensitivity](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html).

## Validation

```bash
python3 -m unittest discover -s script/analysis/tests -p test_download_coloc_gwas.py
Rscript --vanilla script/analysis/tests/test_em_gwas_coloc.R
Rscript --vanilla script/analysis/tests/test_em_coloc_screen.R
```

Tests use pinned public header fixtures and synthetic data: all 32 parsers,
reverse-strand liftover, ZIP/PGC headers, immutable download checksums, exact H0–H4
enumeration, gain/loss symmetry, null evidence, SNP order, causal-SNP loss, absent
CSs, effect-scale invariance, actual SuSiE RSS, and frozen-job disk I/O. Real RCC
colocalizations require the screening bundles, full downloads and a suitable LD
reference; local tests do not establish biological findings.

On this Windows installation Python could not validate UCSC's certificate chain,
while Windows curl verified it successfully. If the RCC Python has the same CA
issue, download the UCSC `hg19ToHg38.over.chain.gz` and `hg18ToHg38.over.chain.gz`
using a client with a working trusted CA store, then set `chain_grch37` and
`chain_grch36` in the download settings (or the corresponding Python CLI flags).
Local chain checksums are recorded too. Do not disable TLS verification.
