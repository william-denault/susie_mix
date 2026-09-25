# Additive simulation: initialization from SuSiE-slide

This focused experiment uses **two true additive causal SNPs**, n = 500,
total PVE = 0.05, and **fitted L = 10** for all methods, matching the existing
simulation settings. Four array tasks run 100 replicates each (400 total).
Seeds are 1000001-1000400. Genotype selection, QC, causal sampling, equal
effect magnitudes with random signs, and LD-adjusted genetic variance match
the pure-additive generator in `sim_workhorse.R`.

Each replicate fits exactly the same additive genotype matrix and phenotype:

1. `susieR::susie(X, y, L = 10)` from its default initialization.
2. `susieSlide::susie(X, y, L = 10)` with estimated genotype sliders.
3. `susieR::susie(X, y, L = 10, model_init = slide_fit)`.

The third fit is an additive refit; it does not optimize slider deltas.
On older susieR installations the runner uses the equivalent `s_init`
argument and records that choice. All fits use `estimate_prior_method =
"optim"`, 95% credible sets, minimum absolute correlation 0.5, tolerance
0.001 and up to 1000 iterations. Namespace-qualified calls select the
intended package. Package versions and settings are saved in checkpoints.

## Run on RCC

Sync the new R script, `simulation_metric_helpers.R`, the existing
`script/scan_tissue_attempt/workhorse_utils.R`, and the job file to the project.
Use a stable directory of PLINK `.raw` files. From the project root:

```sh
sbatch job/run_additive_slide_init
```

The defaults are `/project2/mstephens/wdenault/susie_mix` and its `temp_plink`
subdirectory. Set `SUSIE_MIX_PROJECT_DIR` / `SUSIE_MIX_GENOTYPE_DIR` to override
them. The R/4.2.0 module must have susieR, susieSlide, data.table and matrixStats
installed. The runner checks the packages and initialization API before fitting.

Results go to `simulation results/additive_slide_init_v1/chunks`. Each replicate
is checkpointed. Resubmitting skips successful replicates and retries failures;
incompatible settings, genotype file signatures or packages are rejected.
Do not run two copies of the same chunk concurrently. Errors and nonconvergence
are recorded separately. No simulations are submitted when sourcing the R file.

## Summarize locally after copying the results

```r
project_dir <- "C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix"
Sys.setenv(SUSIE_MIX_PROJECT_DIR = project_dir)
source(file.path(project_dir, "script/sim/sim_additive_slide_init.R"))
summarize_additive_initialization()
```

Alternatively: `Rscript --vanilla script/sim/sim_additive_slide_init.R summarize`.
Outputs are `metric_summary.csv`, `method_comparison.csv`,
`replicate_metrics.csv`, `optimization_comparison.csv` and `errors.csv`.
Coverage/power/purity use the requested denominator-based normal intervals;
CS size uses Gaussian intervals. Successful nonconverged fits remain included
and their counts are reported. Partial checkpoints summarize completed records.

Compare **SuSiE versus SuSiE-init-slide ELBOs** in `optimization_comparison.csv`:
both optimize the same additive model. A positive gain indicates a better
additive objective from slide initialization. The slider has a different model,
so its raw ELBO is not an additive-optimization comparison. Power improvement
after initialization would support an optimization explanation, but the
experiment alone does not establish why the original performance gap occurred.

Checkpoints retain PIPs, credible sets, objective histories, convergence,
iteration counts, variance estimates, timing and slider deltas at causal SNPs.
Use `save_fits = TRUE` in `run_additive_initialization()` if full fitted models
are needed. For a small local run after sourcing:

```r
run_additive_initialization(chunk = 1L, reps_per_chunk = 2L,
  genotype_dir = "path/to/stable/genotypes",
  output_dir = file.path(project_dir, "simulation results/additive_init_pilot"))
```

`fit_L` is adjustable in that function; the true causal count remains two.
