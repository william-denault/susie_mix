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

From R at the project root, source the writer. Sourcing **generates the files**;
it does not fit models or submit jobs:

```r
source("script/sim/write_jobs.R")
```

Alternatively, run `Rscript --vanilla script/sim/write_jobs.R` in a terminal.
The writer creates:

- `script/sim/jobs_slide/sim_job_1.R` through `sim_job_1125.R`. Each is an
  executable R script containing its scenario counts, PVE, seeds, replicate
  count and output filename. It uses `run_job.R` for fitting and checkpointing.
- `job/run_simulation_slide_batch_1` through `run_simulation_slide_batch_4`.
- `conditions.csv`, `manifest.csv` and `submission_batches.csv` in `jobs_slide`
  for inspecting the full design and batch mapping. Generated R scripts carry
  their own settings and do not need the manifest at execution time.

The default batch size is **300 tasks**. On RCC, from the project root, start
the whole sequence with one command:

```sh
sbatch job/launch_simulation_slide
```

The global launcher submits batch 1 and queues one small continuation job. After
all tasks in batch 1 finish successfully, that continuation submits batch 2,
then repeats for batches 3 and 4. It does **not** queue all 1,125 array tasks.
At most one simulation batch, one pending continuation and the briefly running
submission job count toward your quota (302 jobs with a 300-task batch).
No terminal needs to remain open.

This follows the EM launcher's `afterok` dependency pattern. A failed, cancelled
or timed-out task stops the chain. Once the cause is fixed and the previous
jobs have stopped, resume from that batch, for example:

```sh
sbatch job/launch_simulation_slide 2
```

Completed compatible checkpoints are skipped. If only the continuation failed
to submit, the log gives the next batch number to launch after the active array
finishes. Active-job checks and a submission lock prevent duplicate launches
through this global launcher. It records job IDs under
`simulation results/slide_v1/launcher/` and prints a `scancel` command for stopping
future batches while leaving the current array running.

For manual submission instead, the same four launchers are available:

```sh
sbatch job/run_simulation_slide_batch_1  # jobs 1–300
# After batch 1 finishes:
sbatch job/run_simulation_slide_batch_2  # jobs 301–600
# After batch 2 finishes:
sbatch job/run_simulation_slide_batch_3  # jobs 601–900
# After batch 3 finishes:
sbatch job/run_simulation_slide_batch_4  # jobs 901–1125
```

When submitting manually, submit one batch at a time. Other running and pending jobs also consume your
QoS quota; 300 is a configurable batch size, not a guarantee of available slots.
To generate smaller batches, set the option before sourcing:

```r
options(susie.sim.array_batch_size = 50L)
source("script/sim/write_jobs.R")
```

The writer prints the corresponding commands and removes obsolete generated
batch launchers and numbered R scripts if the design or batch count shrinks.
Regenerate only when no jobs from the previous layout are pending or running.
For custom generation without the automatic default run:

```r
options(susie.sim.generate_jobs = FALSE)
source("script/sim/write_jobs.R")
write_simulation_jobs(array_batch_size = 300L)
```

Each batch uses local array indices starting at zero and a built-in offset,
with `job_id = offset + SLURM_ARRAY_TASK_ID + 1`. No array index exceeds 299
with the default plan. `job/run_simulation_slide` is also generated and defaults
to the first batch; it accepts an optional offset for manual submissions.
The legacy `job/run_simulation` and `script/sim/jobs` are separate and are not
used by the new batches.

The launchers default to `/project2/mstephens/wdenault/susie_mix`; set
`SUSIE_MIX_PROJECT_DIR` if the project is elsewhere. Generated R scripts find
the project from their own location, so they can be generated locally and synced
to RCC. Genotypes default to `temp_plink/*.raw`; set `SUSIE_MIX_GENOTYPE_DIR` to
use another genotype directory. With the global launcher, task logs are
`job/logs/susie_slide_sim_ARRAYID_TASKID.out` and `.err`, and continuation logs
are `job/logs/susie_slide_launch_JOBID.out` and `.err`. The initial global launch
logs are written where `sbatch` was invoked, as are task logs for manual batches.
Shell scripts are generated with Unix line endings, including on Windows.

To resume an interrupted batch, repeat its submission command. Completed
compatible checkpoints are skipped. For one job, without submitting an array:

```sh
Rscript --vanilla script/sim/jobs_slide/sim_job_1.R
```

`source("script/sim/jobs_slide/sim_job_1.R")` also runs that individual job.

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
bash script/sim/tests/test_slide_launchers.sh
```

The R test checks every generated job's settings and independently enumerates
the grid; fits all 25 scenario families,
an unequal five-SNP triple, and a high-PVE paired repeat; checks the generating
effects and variance scaling; exercises interrupted and completed checkpoint
resumption through a generated R script; and validates plots and all three
paired method contrasts. The shell test checks batch boundaries, the four-stage
continuation chain, duplicate prevention and failure propagation using fake
Slurm commands and Rscript, without submitting anything to Slurm.
Fixtures and figure previews are written only to `tmp/slide_simulation_validation`.
These are small validation runs, not the full cluster simulation.
