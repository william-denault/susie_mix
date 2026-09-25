# Test SuSiE initialization using the original simulations

The focused experiment calls the existing `sim_mix()` generator in
`script/sim/sim_workhorse.R` with `L_add = 2`, all other causal counts zero,
and PVE 0.025 or 0.05. It uses the original genotype file pool when available;
otherwise it extracts real GTEx genotypes for a gene using the scanning
workhorse's procedure. Both paths then use the original sim_mix() QC, donor
sampling, causal selection, random signs, LD-adjusted effect scaling and
Gaussian noise.

The new optional `return_data = TRUE` argument to `sim_mix()` returns its X/y
immediately before fitting. Its default is FALSE, so the original simulation
jobs retain their existing behavior and output.

The focused runner fits these three methods to each generated dataset:

1. SuSiE with its default initialization.
2. SuSiE-slide.
3. Additive SuSiE on the same original X/y, with `model_init = slide_fit`
   (or the equivalent `s_init` on older susieR versions).

Both baseline calls use the same arguments as `sim_mix()`: fitted L = 10,
standardization, prior method "optim", coverage 0.95, minimum absolute
correlation 0.5 and max_iter = 1000. Package defaults, including convergence
tolerance, are used as in the original workhorse. Slide uses min_obs = 5.
The warm-start call uses the same additive fitting arguments plus the
initialization object. SuSiE-mix is not fitted in this focused comparison.

## Run in RCC RStudio

Upload/extract the current ZIP into the project root on RCC, including the
updated `sim_workhorse.R`. Then run in the R Console:

```r
project_dir <- "/project2/mstephens/wdenault/susie_mix"
Sys.setenv(SUSIE_MIX_PROJECT_DIR = project_dir)
source(file.path(project_dir, "script/sim/sim_additive_slide_init.R"))

results <- run_additive_initialization_experiment()
```

This uses n = 500, two true additive causal SNPs, fitted L = 10, both PVE
values, and 400 replicates per PVE (four chunks of 100). Seeds
1000001-1000400 match the original 400-replicate jobs. Reusing the same seed,
genotype files and R/package versions reproduces the original generator's
dataset at the matching PVE. Across PVEs, genotypes, causal SNPs, signs and
standardized noise are paired.

The genotype directory is resolved exactly as in the original generated jobs:
`SUSIE_MIX_GENOTYPE_DIR`, if set, otherwise the project's `temp_plink` folder.
If your original jobs use an override, use that same value in RStudio.
If that folder contains .raw files, they are used as before. If it is empty or
missing, the runner automatically extracts real cis-genotypes from GTEx. For
each seed it samples an autosomal protein-coding gene from
`data/genes_protein_coding.txt`, obtains its annotation using the scanning
workhorse's parser, and exports SNPs within 500 kb of its TSS with PLINK.
The same seed picks the same gene at both PVEs and on reruns. The selected gene
and coordinates are saved with the results.

The fallback uses the scanning workhorse's existing paths:

- `/project2/mstephens/gtex/plink2`
- `/project2/mstephens/gtex/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.{bed,bim,fam}`
- `/project2/mstephens/gtex/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz`

PLINK uses the workhorse's SNP/duplicate filtering and additive export flags,
with one thread and a 2000 MB workspace. The resulting .raw file is passed to
the original sim_mix() generator in a private temporary folder. Its temporary
export and logs are removed after use, including on failure. They are not added
to the original genotype pool. The genotype data remain real; only the phenotype
is simulated. Since the fallback samples genes afresh, its regions need not be
identical to those in an older cached-file run.

No new arguments are needed for the fallback. If desired, specify one gene
instead of randomly selecting genes (this applies when no cached files exist):

```r
Sys.setenv(SUSIE_MIX_SIM_GENE = "GTF2H2")
results <- run_additive_initialization_experiment()
```

Unset `SUSIE_MIX_SIM_GENE` to restore random gene selection. If the data paths
differ from the scanning workhorse defaults, set `SUSIE_MIX_GTEX_DIR` or
`SUSIE_MIX_PLINK` before running. A change of input mode or selected gene requires
a separate output directory when checkpoints already exist.

Sourcing or pasting the main script only loads functions. The explicit
experiment call runs sequentially inside the existing RStudio session.
It does not submit a Slurm job. Packages are checked before fitting.

To start with the first 100 replicates at each PVE:

```r
results <- run_additive_initialization_experiment(chunks = 1L)
```

Later, repeat the full command to skip successes and retry failures, or use
`chunks = 2:4`. Keep reps_per_chunk = 100 when resuming because it determines
the seed numbering. For a two-replicate pilot, use a separate output folder:

```r
pilot <- run_additive_initialization_experiment(chunks = 1L, reps_per_chunk = 2L,
  output_dir = file.path(project_dir, "simulation results/additive_init_pilot"))
```

## Results

Checkpoints and summaries go to `simulation results/additive_slide_init_v3/`.
Each PVE has its own `pve_0.025/` or `pve_0.05/` subfolder. This separates the
original-generator experiment from the previous runner's outputs.

Each replicate is checkpointed. Successful records are skipped on rerun;
failures are retried. Changes in settings, genotype file signatures or package
versions require a new output directory. Do not run the same PVE/chunk twice
concurrently. All three methods use the same X and y within a replicate.

```r
results <- summarize_additive_experiment()
results$summary
```

The per-PVE and combined CSV files are `replicate_metrics.csv`,
`metric_summary.csv`, `method_comparison.csv`, `optimization_comparison.csv`,
and `errors.csv`. Coverage, power and purity use the requested denominator-based
normal intervals; CS size uses Gaussian intervals. Successful nonconverged fits
remain included and their convergence counts are reported.

Compare SuSiE versus SuSiE-init-slide in `optimization_comparison.csv`:
a higher additive ELBO after initialization indicates a better solution for
the additive objective. The saved results also include PIPs, credible sets,
ELBO histories, convergence, timings and causal slider estimates.
Set `save_fits = TRUE` to retain complete models as well.

## Optional Slurm execution

```sh
sbatch job/run_additive_slide_init
```

Tasks 1-4 run chunks 1-4 at 2.5% PVE; tasks 5-8 run them at 5% PVE.

```sh
Rscript --vanilla script/sim/sim_additive_slide_init.R 1 100 0.025
Rscript --vanilla script/sim/sim_additive_slide_init.R 1 100 0.05
Rscript --vanilla script/sim/sim_additive_slide_init.R summarize
```

## Validation

The focused tests compare X, y, causal SNPs and effects against the original
generator and compare both baseline fits against the original `sim_mix()`
route. They also verify the direct warm-start call, pairing across PVEs,
checkpoint resume/retry, summaries and Console/source loading. GTEx extraction
tests use fixture annotations and a mocked PLINK export to check selection,
arguments, temporary-file cleanup, reproducibility and unchanged phenotype
generation. Running against the actual GTEx files requires the RCC environment.
