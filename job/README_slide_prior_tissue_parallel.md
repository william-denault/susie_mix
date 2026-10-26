# Optional tissue jobs for long slider-prior genes

This patch adds an opt-in execution path to the existing shared EM controller.
It does not submit another controller, change the 20-update target, cancel jobs,
expand the chromosome scope, or change the slider model, convergence tolerance,
prior weights, component-count rule, or simulation sources.

## Install, inspect, then select genes

Upload `output/slide_prior_tissue_parallel_rcc_2026-10-26.zip` into the RCC project
and run:

```bash
cd /project2/mstephens/wdenault/susie_mix &&
unzip -o slide_prior_tissue_parallel_rcc_2026-10-26.zip &&
python3 updates/slide_prior_tissue_parallel/install.py "$PWD" &&
module load R/4.2.0 &&
Rscript --vanilla script/scan_tissue_attempt/run_slide_prior_parallel.R preflight "$PWD" latest
```

The installer validates every payload and installed source before any replacement.
It accepts the October 23 shared bundle or the optional October 25 target patch,
and backs up replaced files. Each replacement is atomic, with helpers and the
worker installed before the CLI. Repeating the same installation is safe.
The target-control helper is included as a dependency, but installation never
writes `run_control.rds` or changes the target to 40.

Show measured whole-gene runtimes for a completed iteration (for example 2):

```bash
Rscript --vanilla script/scan_tissue_attempt/run_slide_prior_parallel.R split-report "$PWD" 2
```

An unfinished receipt is not proof of a timeout: inspect its Slurm log. Select
specific eligible, previously fitted genes using their exact manifest names:

```bash
Rscript --vanilla script/scan_tissue_attempt/run_slide_prior_parallel.R split-enable "$PWD" latest GENE_NAME
```

Replace `GENE_NAME` with a real name from the report; several names may be supplied.
The next existing controller adopts the selection. An already submitted whole-gene
task finishes normally. If it succeeds, the selected gene uses tissue jobs starting
with its next EM iteration. If it fails, tissue preparation becomes eligible when
the whole-gene task leaves the queue. Previous attempt receipts are preserved;
each new prepare/tissue/assembly unit has its own three-attempt budget. Selection
persists for subsequent iterations and can be extended by repeating split-enable.
Do not submit an additional EM controller or edit the frozen scope files.

## Execution and saved files

All units are ordinary `slide_parallel_fit` array tasks counted by the existing
shared scheduler. The 300-job total budget, up to 298 combined workers, and
five-minute refill eligibility remain unchanged. Prepared inputs are made once
per gene per iteration, avoiding repeated loading of the full expression dataset
and repeated PLINK extraction in every tissue worker.

```text
results_slide_prior_em/iteration_NNN/
  tissue_results/GENE/
    manifest.rds                 frozen inputs, sources, tissues and checksums
    tissues.csv                  mapping of numbered units to tissue names
    inputs/tissue_001.rds         prepared data and this tissue's warm start
    results/tissue_001.rds        validated, converged tissue result
    completed/tissue_001.rds      compact receipt with checksum and file stamp
    attempts/prepare.rds
    attempts/tissue_001.rds
    attempts/assemble.rds
  results/GENE.rds                usual combined gene result, written last
```

Preparation uses the unchanged QC and normalization operations and exactly the
preceding iteration's fitted tissues, samples and SNPs. Independent tissue workers
then call the same package fitter with the current iteration's frozen tissue prior.
Their successful outputs survive failures of other tissues. Input changes,
unconverged fits, foreign tissues, and exhausted retries fail closed.

Tissue jobs share a gene lock and each hold a separate exclusive unit lock.
Whole-gene, preparation and assembly jobs hold the exclusive gene lock. The planner
also checks ownership, permits different tissues of the same gene concurrently,
and prevents duplicate units or overlap with a whole-gene job.

Assembly runs only after all expected tissue tasks leave the queue and their
completion receipts validate. It rechecks the tissue results and atomically writes
the original gene-level list format with the original metadata. The existing gene
validator then collects the EM counts. No partial gene contributes to a prior
update. The whole iteration still waits for every required gene, so one genuinely
slow tissue can delay its M-step, but cannot discard completed neighboring tissues.

Each unit retains the existing 48 GB / 23-hour allocation. Splitting does not solve
a single tissue that itself exceeds this allocation; that requires a separately
reviewed resource change or within-fit checkpointing. After three failures the
controller stops and retains completed work. No gene or tissue is automatically
excluded. Prepared inputs and tissue results use additional disk space and are
retained for diagnosis; this patch does not delete them.

The regional posterior and EM count calculations are unchanged. CS purity may
use random subsampling for large credible sets in the installed package. Separate
tissue processes reset to the saved preparation RNG state, so those sampled purity
annotations (and threshold-based CS reporting) can differ from a sequential gene
run. PIPs, posterior counts and EM weights are independent of that post-fit sample.
Retries of the same tissue use the same saved RNG state.

## The 243 old files are a different issue

The saved pre-scope audit records 243 PLINK extraction failures, all with zero fits:
220 X genes, 11 Y genes, and 12 MT genes. They are not successful slow fits.
Some logs report that no variants survive the chromosome/SNP filters; others have
a different PLINK exit status. Their presence as RDS files did not mean they had
entered the EM. These genes are among 851 genes excluded by the frozen X/Y/MT policy.

This patch intentionally refuses split-enable for those excluded genes. Their
recovery first needs usable genotype inputs and chromosome-appropriate modeling.
The input files on RCC are not available in the synced local workspace. A first,
read-only check on RCC is:

```bash
awk '{n[$1]++} END {for (c in n) print c, n[c]}' \
  /project2/mstephens/gtex/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.bim | sort -V
```

If only autosomes are present, a genotype source including the relevant chromosomes
is needed. Chromosome X requires sex and PAR-aware handling; the old PLINK log
reports all 866 samples as ambiguous sex. The current downstream MAF/HWE calculation
assumes diploid autosomal genotypes. Haploid variants without heterozygotes cannot
inform the present heterozygote-slider prior: min_obs forces them additive and their
counts are excluded. More workers cannot create that information.

An expanded cohort must have its own versioned EM history. One possible initialization
is the current learned prior, with current autosomal fits reused as warm starts;
all eligible gene/tissue pairs must then be fitted under those same frozen priors
before pooling. This avoids replaying every old update, but is a new EM run, not a
retroactive correction of existing history. Alternatively, unchanged autosomal
iteration-000 fits could be reused with new genes fitted under the original uniform
prior. Neither expansion is activated by this scheduling patch.

## Validation

The real-fitter test compares the refactored worker against the deployed October 23
worker, then compares whole-gene and separate-tissue outputs, EM counts and priors
on synthetic GTEx inputs. It also tests a failed tissue, nonconvergence, exhausted
retries, preserved completed tissues, active whole-gene ownership, concurrent tissue
jobs, the assembly barrier, frozen cohort and changed-input detection. Shell tests
cover worker dispatch and locks; existing controller tests cover the shared budget.
These local tests do not substitute for a pilot on a selected real RCC gene.
