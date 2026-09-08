# GTEx simulation benchmark

Compare **SuSiE** (additive predictors) and **ordinary SuSiE-mix** (additive,
recessive and dominant predictors with uniform predictor priors). Both use
**N = 200 and L = 10 in every scenario**. The weighted mixture is excluded.

| Generating effects | Number of distinct causal SNPs |
| --- | --- |
| Additive only | 1, 2, 3, 4, 5 |
| Recessive only | 1, 2, 3, 4, 5 |
| Dominant only | 1, 2, 3, 4, 5 |
| Additive + recessive | 2, 3, 4, 5 |
| Additive + dominant | 2, 3, 4, 5 |
| Recessive + dominant | 2, 3, 4, 5 |
| Additive + recessive + dominant | 3, 4, 5 |

Cross each row/count with PVE = 0.05, 0.10, 0.20: **90 scenarios**, each with
at least **200 replications**, giving **18,000 datasets and 36,000 fits**.
Each mixed scenario contains at least one SNP of every named effect type.
Each causal SNP has one generating coding. There are no interaction effects.

## Genotypes and phenotype generation

For every replication independently:

1. Draw a gene uniformly from the project's protein-coding GTEx gene list,
   restricted to autosomes and genes in the workhorse's annotation helper.
   As in the workhorse, use the first annotation for duplicate gene symbols.
2. Extract the TSS +/- 500 kb region from the same GTEx PLINK files. Preserve
   `--snps-only --max-alleles 2 --rm-dup exclude-all --maf 0 --recode A`.
3. Across all genotype donors, remove SNPs with any missing calls or no variation,
   then call **the existing `qc_filter_geno()`**: orient to the minor allele,
   retain MAF strictly greater than 0.05, and exclude HWE chi-square P < 1e-8.
4. Draw 200 distinct genotype donors without replacement. Retain the full-donor
   allele orientation. Use **the existing `recode_snp_matrix()`**, then apply the
   workhorse's sample-specific filters separately to each predictor: column sum
   >= 5 and positive SD. For recessive predictors, this requires >= 5 minor-allele
   homozygotes. No additional MAF/HWE or LD pruning is applied to the 200 donors.
5. Draw coding assignments uniformly conditional on all requested types appearing.
   Draw causal SNPs uniformly from the retained predictors of each assigned type,
   without reusing a biological SNP. Thus causal recessive effects are conditional
   on surviving the workhorse's recessive filter. Rejected loci/samples are logged.
6. Give each causal predictor equal effect magnitude **after scaling by its sample
   SD**, with independent random positive/negative signs. Rescale the total genetic
   value, accounting for LD, so `var(g) = PVE`. Generate independent normal errors
   with variance `1 - PVE`, then `y = g + error`.
7. Fit both methods to the **same y and donor sample**, using the full retained
   additive or mixed matrix. Both use L = 10, `standardize = FALSE`, EM prior
   estimation, 95% CSs, and `min_abs_corr = 0`, as in the workhorse. Both are allowed
   1,000 iterations (more than the workhorse's implicit default) with tolerance
   0.001. Neither is told the true causal count, effect types, or residual variance.

This is a genotype-based simulation, without tissue-specific donor selection,
expression normalization or sex residualization: the simulated outcome is already
a quantitative phenotype with a specified generating model. Sex chromosomes are
excluded because they require different genotype/ploidy handling.

PVE is `var(g) / (var(g) + sigma_error^2)` conditional on the sampled genotypes.
The *realized* `var(g) / var(y)` and `cor(g, y)^2` fluctuate at N = 200; both, along
with genetic/error covariance, are saved. Errors are not projected off genotypes
to force exact realized PVE, which would alter their independence.

## Metrics

The primary evaluation unit is the **biological SNP**. For SuSiE-mix, sum coding
probabilities within each single-effect component, then combine components:

`PIP(SNP j) = 1 - product_l(1 - sum_{coding c of j} alpha[l,c])`.

Only components with estimated prior variance > 1e-9 contribute, consistent with
the [SuSiE PIP accessor](https://stephenslab.github.io/susieR/reference/susie_get_methods.html).
Do not add or maximize the final coding-specific PIPs. Predictor PIPs are also saved.

- **Power:** fraction of true causal SNPs with SNP PIP >= 0.95. Threshold configurable.
- **CS power:** fraction of causal SNPs appearing in the union of reported CSs.
  This can be high for large diffuse CSs; interpret alongside size and purity.
- **CS size:** number of distinct biological SNPs per reported CS; also report
  number of coding-specific predictors. Mixed CSs are mapped from the fitted
  predictor CSs, not rebuilt at SNP level.
- **Coverage:** fraction of reported CSs containing at least one true causal SNP.
  Also report correct SNP-and-coding coverage. A CS containing two causal SNPs is
  counted once. This is empirical set coverage, not effect-size interval coverage
  or the nominal posterior mass. Coverage/size/purity are undefined for a fit with
  no CS; CS power is zero. Both pooled-CS and equal-locus summaries are saved.
- **Purity:** minimum absolute correlation among CS members. Report native fitted
  coding purity and additive-genotype LD purity over the unique SNPs, giving both
  methods a common biological-SNP measure. Like SuSiE's default, sets larger than
  100 members use a reproducible random subset; approximation indicators are saved.
  Increase `n_purity` for more accurate minima. Singleton purity is 1. No purity
  threshold is imposed, matching the workhorse.
- **PIP calibration:** bin *all* SNP PIPs into 20 bins; compare mean PIP with the
  empirical causal fraction, recording counts and Brier scores. Zero and one are
  included. Small bins can be noisy. Monte Carlo uncertainty uses replicate-level
  clusters, accounting for within-locus dependence.
- **Power vs FDR:** sweep SNP PIP thresholds from 0 to 1 (plus 0.999). For each
  replication, FDP = false discoveries / max(1, discoveries). FDR is the mean FDP
  over replications, including zero-discovery replicates; power is mean causal
  recall. Pooled false discoveries / discoveries is reported separately. Plot
  observed FDR, not `1 - threshold`; empirical curves need not be monotonic.

The main paired summaries include a replication only when **both methods converge**.
Completion requires the requested number of paired successes for *every* scenario,
without failures. Failures/nonconvergence are saved, reported and cause nonzero exit;
they are never replaced with easier loci. Sparse or difficult scenarios are not
silently dropped. Monte Carlo SEs and paired SuSiE-mix minus SuSiE differences are
reported. Figure error bars are mean +/- 1.96 Monte Carlo SE, clipped to valid bounds.

## Running

Required R packages: `susieR`, `data.table`, `matrixStats`, `R.utils`. Production also
requires PLINK2, the GTEx PLINK dataset and GTF used by the workhorse. Paths default
to `/project2/mstephens/gtex`. Set `GTEX_DATADIR` and `PLINK2` or supply an R config
file. Run commands below from the project root; R script paths themselves also work
from other directories.

```sh
# Write the 90-scenario plan without accessing genotypes or fitting models.
Rscript --vanilla script/simulation/run_simulation.R --dry-run

# Serial production run; completed replication checkpoints are skipped on restart.
Rscript --vanilla script/simulation/run_simulation.R

# Preferred cluster launch: 200 array jobs, maximum 20 simultaneous jobs.
# Run sbatch FROM THE PROJECT ROOT. The 'simulation results' folder must exist.
sbatch script/simulation/submit_simulation.sh

# First check seven architectures in a single RCC pilot job (14 fits).
sbatch --array=1 --time=01:00:00 script/simulation/submit_simulation.sh \
  --pilot --replications 1 --n-shards 1 --scenario 1,16,31,46,58,70,82

# After all array jobs finish, aggregate all 90 scenarios and produce PNG figures.
Rscript --vanilla script/simulation/summarize_simulation.R

# Small local validation on the existing single-locus extract (NOT a benchmark).
Rscript --vanilla script/simulation/run_simulation.R --pilot --replications 1 \
  --scenario 1,16,31,46,58,70,82 --raw-file plink2.raw
Rscript --vanilla script/simulation/summarize_simulation.R \
  --output "simulation results/pilot" --allow-incomplete

# Scientific correctness tests, including real L=10 model fits.
Rscript --vanilla script/simulation/tests/test_simulation.R
```

`--pilot` permits fewer than 200 replications and writes to `simulation results/pilot`.
`--raw-file` is allowed only in pilot mode because a single locus cannot establish
performance across random GTEx genes. All pilot and partial figures are labelled.
For a cluster pilot across random genes, omit `--raw-file`.

Other arguments: `--output PATH`, `--config FILE`, `--replications INTEGER`,
`--scenario 1,2,...`, `--replicate 1,2,...`, `--shard-id I --n-shards M`, `--seed INTEGER`,
`--save-fits`, and `--retry-failed`. Scenario IDs are in the generated plan.
Shards partition trial IDs, so all array workers cover disjoint trials. Keep the
same config and output directory across shards. The defaults produce 90 trials
per shard; resource needs depend on locus size. Adapt partition/memory/time to the
cluster if necessary. Fit objects are optional because they can be large.

A custom config file modifies the existing `config` list, for example:

```r
config$genotype_prefix <- "/path/to/GTEx_genotypes"
config$gtf_file <- "/path/to/annotation.gtf.gz"
config$plink_exec <- "/path/to/plink2"
config$replications <- 300L
```

Seeds depend on scenario and replication, not scheduling order. Successful
checkpoints are never rerun. `--retry-failed` uses the *same seed and data draw*;
changing the model/data/code requires a new output directory. Source hashes,
package versions, input paths/sizes/timestamps (plus hashes for small inputs),
session information, simulation truth, sampled donors, and PVE diagnostics are saved.
Large genotype files are not fully hashed. A per-trial lock prevents overlapping
workers from writing the same result. After a killed job, remove its stale `.lock`
directory only after verifying that no worker is running that trial.

## Outputs

Under `simulation results/production` (or `pilot`):

- `plans/`: the complete scenario grid, replication target, N and L.
- `manifest_shard_*.rds`: configuration, code/data provenance and R environment.
- `replicates/scenario_*/rep_*.rds`: gene, donor IDs, filtering counts, truth,
  phenotype, genetic value, PVE checks, method warnings/status, CSs and all PIPs.
- `summary/completion.csv`, `STATUS.txt`, `replicate_status.csv`, `method_status.csv`:
  completion/failure audit; inspect these before interpreting results.
- `summary/metrics_summary.csv`, `paired_differences.csv`, `cs_summary.csv`,
  `power_fdr.csv`, `pip_calibration.csv`: main performance tables.
- `summary/replicate_metrics.csv`, `credible_sets.csv`, `phenotype_variance.csv`,
  `rejected_loci.csv`: detailed audit tables.
- `summary/figures/`: PNG plots by PVE and architecture, with separate K panels/files
  for calibration and power-FDR curves. Aggregation writes partial outputs for
  inspection but exits with an error unless complete or `--allow-incomplete` is set.

Generated runs are ignored by Git. Production GTEx inputs are not bundled with
this project, and local validation does not constitute the 18,000-dataset benchmark.
