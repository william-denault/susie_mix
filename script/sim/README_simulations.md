# Five-effect simulations with SuSiE-slide

The simulation compares `susieR::susie()` on additive genotypes (SuSiE),
`susieR::susie()` on the existing additive/recessive/dominant predictor blocks
(SuSiE-mix), and `susieSlide::susie()` on the original genotypes (SuSiE-slide).
All methods fit the same phenotype and use fitted L = 10, 95% credible sets,
minimum CS absolute correlation 0.5, `estimate_prior_method = "optim"`, and
`max_iter = 1000`. Namespace-qualified calls prevent the packages' two `susie`
functions from masking one another.

## Design

| Generating effect | Count argument | Delta | Relative genotype effects (0, 1, 2) |
|---|---|---:|---|
| Additive | `L_add` | 0 | 0, 1, 2 |
| Recessive | `L_rec` | -1 | 0, 0, 2 |
| Dominant | `L_dom` | 1 | 0, 2, 2 |
| Partially recessive | `L_prec` | -0.5 | 0, 0.5, 2 |
| Partially dominant | `L_pdom` | 0.5 | 0, 1.5, 2 |

Partial effects use `x + delta * (x == 1)`. The original recessive and dominant
generating columns remain their 0/1 indicator codings; the factor of two in the
table has no effect after each causal column is standardized.

Every allocation of **1–5 distinct causal SNPs** among **at most three of the
five effect types** is included: five pure families, ten pairs, and ten triples.
There are 5, 15, 35, 65, and 105 allocations at K = 1, 2, 3, 4, and 5,
respectively, for **225 configurations**. All 55 original configurations remain.

Each configuration has n = 500, PVE = **0.05, 0.10, 0.20, 0.30, 0.40**, and
400 replicates: **1,125 jobs, 450,000 datasets, and 1,350,000 fits**.
PVE is total genetic PVE. Causal predictors have equal individual contribution
variances and random effect signs. The combined genetic signal is rescaled to
`var(g) = pve`, accounting for LD; independent Gaussian errors have variance
`1 - pve`. Finite-sample realized `var(g) / var(y)` need not equal the target.
The same seed is reused across conditions and PVE values, as in the original design.

Genotype filtering and the original scenarios' causal eligibility rules are
preserved. When partial effects are present, candidate causal SNPs additionally
require at least five observations in **each** genotype class after donor
subsampling, so their intermediate shape is identifiable. SuSiE-slide uses
`min_obs = 5`: at other candidate SNPs with a rare/absent genotype class it
forces delta to zero. This status is saved rather than dropping those SNPs.
The slider uses its native fixed additive-SD scaling; SuSiE-mix uses its native
per-column standardization. Generating predictors are standardized separately.

## Prepare and run

Install `susieR`, the implemented **susieSlide** package, `data.table`, and
`matrixStats` in the R library used by the cluster's `R/4.2.0` module.
The local slider source is the root package in the `susieR` repository on the
`susie_slide` branch; its package name is `susieSlide`, and its fitting entry
point is `susieSlide::susie()`. The runner checks required packages before any
replicates are started. Local installation does not install it on the cluster.

### Check your submission quota first

`QOSMaxSubmitJobPerUserLimit` means the requested tasks would exceed the QoS
limit on your running **plus pending** jobs. Each array task counts separately.
An array `0-200` contains 201 tasks; appending `%20` limits concurrency but does
not reduce the number submitted. See [Slurm resource limits](https://slurm.schedmd.com/resource_limits.html)
and [array limits](https://slurm.schedmd.com/job_array.html).

On RCC, inspect your actual quota and existing tasks:

```sh
rcchelp qos
squeue -u "$USER" --array -h -o '%i %q %T'
```

If `rcchelp` is unavailable on Midway2, try `accounts qos`. RCC documents
`rcchelp qos` in its [FAQ](https://docs.rcc.uchicago.edu/slurm/faq/).
Choose a batch size no larger than the available slots for the applicable QoS.
Even a single task can be rejected while the quota is full; wait for existing
jobs to finish. Failed submissions in this state do not start simulations.

### Submit the slide simulations

The dedicated launcher is **`job/run_simulation_slide`**. `job/run_simulation`
is currently the legacy launcher and must not be used for these manifest jobs.
Copy the new launcher to RCC along with the updated simulation scripts.

The default plan now uses **50 tasks per batch** (22 batches of 50, then 25).
This is a configurable batch size, not an assertion about RCC's actual quota.
From the project root, generate the manifest and, if 50 slots are available,
submit the first batch:

```sh
Rscript script/sim/write_jobs.R
sbatch --array=0-49 job/run_simulation_slide 0
```

**Wait until that batch finishes**, then submit the second batch:

```sh
sbatch --array=0-49 job/run_simulation_slide 50
```

Continue through the offsets printed by the generator: 100, 150, ... 1050.
**Wait for each batch to finish** before submitting the next. The final batch is:

```sh
sbatch --array=0-24 job/run_simulation_slide 1100
```

The complete generated plan covers manifest jobs **1–1125** exactly once.
Each default array has at most 50 tasks and uses indices no higher than 49. The final
argument is the manifest offset: `job_id = offset + SLURM_ARRAY_TASK_ID + 1`.
Submit one batch at a time so these simulations do not queue more than 50
tasks at once. If other jobs use your quota, use a smaller batch size below.
Submitting `sbatch job/run_simulation_slide` defaults to the first 50 tasks only.
To resume an interrupted batch, repeat its command; completed checkpoints are
detected and skipped by the R runner.

The generator writes `script/sim/jobs_slide/conditions.csv`, `manifest.csv`,
and `submission_batches.csv`, and prints the batch submission commands.
The checked-in files are already generated for the design above. The launcher
calls `script/sim/run_job.R` with the mapped manifest job and uses separate log
files per array task. Set `SUSIE_MIX_PROJECT_DIR` if its project path differs from
`/project2/mstephens/wdenault/susie_mix`. Genotypes default to `temp_plink/*.raw`.
The old `script/sim/jobs/sim_job_*.R` files are legacy jobs and are no longer
used by this launcher.

In R, the generator can also be used explicitly:

```r
source("script/sim/write_jobs.R")
write_simulation_jobs()  # Call after sourcing; sourcing alone does not generate jobs.
# If needed, print and save a plan with smaller submission batches:
# write_simulation_jobs(array_batch_size = 20L)
```

For one job, without submitting the array:

```sh
Rscript script/sim/run_job.R 1 /path/to/susie_mix /path/to/genotypes
```

## Checkpoints and output

New saves go to `simulation results/slide_v1/chunks`, leaving the original
`simulation results/chunks` intact. Each filename encodes all five causal
counts, sample size, fitted L, PVE, seed base, replicate count, and chunk.
Each checkpoint includes the design, genotype directory, and package versions.
Resumption requires an exact match and valid three-method records. Changing
the design or packages requires a separate results directory; do not resume
old two-method saves as if the slider had already been fitted.

Successful records contain SNP-level PIPs and CSs for all three methods,
convergence and discovery metrics, true effect labels and deltas, and slider
CS deltas, component-by-causal-SNP deltas, and genotype-count fallback flags.
Phenotypes, genotypes, and full fitted objects are not saved. `genetic_variance`
records the achieved `var(g)` for a scaling check.

`true_pos` always identifies biological causal SNPs. `true_pos_mix` is **NA for
partial effects**, because no exact partial predictor exists in the fitted
SuSiE-mix dictionary. Partial effects must not be relabelled as their nearest
endpoint. SNP-level metrics remain well-defined. The older
`plot_coding_mixture.R` and `plot_lead_pip_coding.R` are still analyses of the
original three-coding checkpoints; their exact-coding metrics are not extended
to these partial-effect scenarios.

Simulation errors are saved with their seed and count toward the requested
replicate total, as before. Nonconvergence is recorded separately. Plotting
audits expose failures and incomplete cells; error records are not interpreted
as null fits. See [README_plots.md](README_plots.md) for the three-method figures.

## Local verification

With the packages available in `.libPaths()`, run from the project root:

```sh
Rscript script/sim/tests/test_slide_simulation.R
```

The test independently enumerates the grid; fits all 25 scenario families,
an unequal five-SNP triple, and a high-PVE paired repeat; checks the generating
effects and variance scaling; exercises interrupted and completed checkpoint
resumption; and validates plots and all three paired method contrasts.
Fixtures and figure previews are written only to `tmp/slide_simulation_validation`.
These are small validation runs, not the full cluster simulation.
