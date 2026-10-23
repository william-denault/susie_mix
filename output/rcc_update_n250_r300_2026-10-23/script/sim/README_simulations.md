# Fresh n=250 simulation benchmark

Prepared 23 October 2026. This replaces the active n=500, 400-replicate,
five-method workflow. Its results and source snapshots are preserved in
`archive/simulations_before_n250_2026-10-23/`.

## Design

- n = **250** individuals per dataset; **300** replicates per condition.
- Exactly **SuSiE**, **SuSiE-mix**, and **SuSiE-slide-prior** are fitted.
- PVE = **10%, 20%, 30%, 40%**. There is no 5% condition.
- **1–5 distinct causal SNPs**, with fitted upper bound **L = 10**.
- All **225 allocations** among at most three of five generating effect types:
  additive, recessive, dominant, partial recessive and partial dominant.
- 5 pure families, 10 pair families and 10 triple families; the allocation
  counts for K=1..5 are 5, 15, 35, 65 and 105.
- **900 conditions, 270,000 datasets, 810,000 fits**, in three 300-job batches.
- Seeds 1,000,001–1,000,300 are shared across conditions.

The generating/QC rules remain unchanged: minor-allele orientation, MAF >= .05,
HWE threshold 1e-8, coding-specific variability/count filtering, distinct causal
SNPs, equal standardized effect magnitudes with random signs, and joint scaling
to var(g)=PVE. Independent Gaussian noise has variance 1-PVE. Partial effects
use x + delta*(x==1), with delta=-.5 or +.5, and require five observations in
each genotype class at causal SNPs. Exact-coding truth is unavailable for
partial effects in SuSiE-mix; biological-SNP evaluation remains valid.

All methods use standardize=TRUE, estimate_prior_method="optim", max_iter=1000,
coverage=.95 and min_abs_corr=.5. SuSiE-mix uses its three original coding
blocks and default uniform predictor priors. SuSiE-slide-prior uses min_obs=5,
the fixed grid seq(-1,1,length.out=17) and uniform fixed probabilities 1/17.
There is no continuous-slide fit, initialized-additive fit or EM weight update.

## Genotype source

The production runner explicitly selects random GTEx gene regions, even if
an old temp_plink pool is present. Each seed selects an autosomal protein-coding
gene from data/genes_protein_coding.txt, then extracts TSS +/-500 kb using the
existing GTEx v8 genotype, annotation and PLINK paths. Each extraction is
private and cleaned up. The chosen gene and region are recorded.

Use SUSIE_MIX_GTEX_DIR and SUSIE_MIX_PLINK only if the RCC paths differ.
SUSIE_MIX_PLINK_MEMORY_MB defaults to 8000. Leave SUSIE_MIX_SIM_GENE unset for
random genes. SUSIE_MIX_SIM_INPUT=raw is an explicit fixture/alternative-study
override; it is not the requested benchmark and preflight rejects it.

## Prepare and submit on RCC

See ../../RCC_SIMULATION_UPDATE.md for the complete upload procedure.
After old simulation jobs have finished or been cancelled, from the project root:

```sh
module load R/4.2.0
Rscript --vanilla script/sim/write_jobs.R
Rscript --vanilla script/sim/preflight.R
sbatch job/launch_simulation_slide
```

The R library needs susieR, susieRSlidePrior, data.table and matrixStats.
The existing finite-prior package must expose min_obs, delta_grid and delta_prior.
Preflight checks versions and input metadata without extracting genotypes or
fitting models. Generating jobs also does not submit or fit anything.

The launcher submits one array at a time and chains the next with afterok.
There are at most 300 array tasks plus continuation/submission jobs for this
run. Other running jobs still count against your account quota. A failed task
stops the chain. Resume with `sbatch job/launch_simulation_slide BATCH_NUMBER`.
Completed compatible checkpoints are skipped; unfinished ones resume by seed.
For the alternative allocation use:
`SBATCH_ACCOUNT=pi-wdenault sbatch --export=ALL job/launch_simulation_slide`.

## Results and reproducibility

New saves are isolated in **simulation results/n250_r300_v1/chunks/** with
schema **4**, n and replicate count in filenames. Old checkpoints cannot resume.
Checkpoints record the design, three methods, package versions, R version,
RNG kind, fitting-source checksums, input signatures and gene pool. A mismatch
requires a separate output directory. PLINK workspace size alone can change
on resume. Keep source files stable while the batch chain is running.

Compact records contain PIPs, biological-SNP CS maps, convergence, ELBO,
residual variance, truth, fixed-prior delta diagnostics and generating variance.
They do not contain full genotypes, phenotypes or fitted models. Extraction
errors stop before a replicate is marked done; simulation/fitting errors are
saved and count toward 300 attempted replicates, as in the earlier benchmark.
Plots exclude these error records consistently across methods and audit them.

## Plot and validate

`source("generate_sim_figures.R")` produces figures in the new results folder.
See README_plots.md for the point-panel style and metric definitions.

Local checks use synthetic genotypes, not the protected RCC data:

```sh
Rscript --vanilla script/sim/tests/test_slide_simulation.R
Rscript --vanilla script/sim/tests/test_slide_genotypes.R
Rscript --vanilla script/sim/tests/test_metric_intervals.R
Rscript --vanilla script/sim/tests/test_plot_calibration.R
bash script/sim/tests/test_slide_launchers.sh
```

The main check independently enumerates the design, verifies every generated
driver, fits all 25 scenario families plus K=5/high-PVE controls, reconstructs
generating effects, checks save/resume/incompatible checkpoints, and exercises
plots and paired metrics. The genotype test mocks PLINK but runs real model
fits. The launcher test mocks Slurm and submits no real jobs.
