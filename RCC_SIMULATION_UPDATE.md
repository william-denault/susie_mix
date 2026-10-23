# RCC update: n=250, 300 replicates, three methods

The prepared update bundle is `rcc_update_n250_r300_2026-10-23.zip` at the project
root. It contains only the updated simulation/plotting sources, 900 generated
jobs, launchers, required shared helpers, gene list and documentation. It
contains no old results and no test/synthetic results.

The local backup is `archive/simulations_before_n250_2026-10-23/`: 1,061 result
files (35.90 GiB), frozen current scripts and publication panels, and a separate
September source snapshot. This is a local archive; it has not been uploaded
to RCC. Check that any RCC-only results you want are backed up before removing
them. No RCC jobs or files were changed by preparing this update.

**Pushing the code to GitHub does not back up the raw simulation archive:**
those large result directories are ignored by Git. Preserve a verified copy of
the full local archive on independent storage before deleting the last RCC copy.
See `SIMULATION_ROLLBACK.md` for code rollback and data recovery instructions.

## Apply the update

1. Finish or cancel old simulation jobs before replacing their shared scripts.
2. Upload the zip to `/project2/mstephens/wdenault/susie_mix/` and extract it
   there with `unzip -o rcc_update_n250_r300_2026-10-23.zip`.
3. From that project root, run:

```sh
module load R/4.2.0
unset SUSIE_MIX_SIM_INPUT SUSIE_MIX_SIM_GENE
Rscript --vanilla script/sim/write_jobs.R
Rscript --vanilla script/sim/preflight.R
sbatch job/launch_simulation_slide
```

Regenerating the jobs removes obsolete numbered jobs above 900 and stale batch
launchers. Preflight checks the three-method setup and real GTEx inputs without
extracting data or fitting models. The launcher automatically runs batches
1–300, 301–600 and 601–900 in sequence. There is no batch 4.

Required packages: susieR, susieRSlidePrior, data.table and matrixStats. The
continuous susieSlide package is not needed for this benchmark. Keep the
finite-prior grid uniform and fixed; no coding-prior EM is run.

If needed for your allocation, replace the final command with:

```sh
SBATCH_ACCOUNT=pi-wdenault sbatch --export=ALL job/launch_simulation_slide
```

The new output is `simulation results/n250_r300_v1/`. It cannot accidentally
resume the old n=500 runs. Removing the old RCC results is not required for
correctness; it is your storage decision. Keep the GTEx source data, gene list,
PLINK executable and annotation file. Cached temp_plink files are not used by
the new default genotype sampling.

## Generate figures

After downloading the new chunks into the same project-relative location:

```r
source("generate_sim_figures.R")
```

Outputs go to `simulation results/n250_r300_v1/figures/`. Four PVE columns,
K=1–5, 25 scenario families, three methods; September-style point panels.
The full design is 225 allocations x 4 PVEs x 300 datasets = 270,000 datasets
and 810,000 fits. Fitted L stays 10, distinct from the true causal count K.
