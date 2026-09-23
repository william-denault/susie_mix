# Simulation and workhorse audit — 23 September 2026

Neither run is complete. This audit uses the locally synchronized logs and checkpoint inventory, with logs through 22 September. It does not query the live cluster scheduler.

## Simulation

Expected design: 1,125 jobs, 400 replicates each, 450,000 datasets, three fitted methods per successful dataset.

| Batch | Array ID | Expected jobs | Checkpoints present / jobs reaching 400 records | Replicate errors |
|---|---:|---:|---:|---:|
| 1 | 49041298 | 300 | 300 | 22 |
| 2 | 49045474 | 300 | 300 | 272 |
| 3 | 49047600 | 300 | 97 | 32,194 |
| 4 | No array submission in local logs | 225 | 0 | Not run in this snapshot |
| Total | | 1,125 | 697 | 32,488 |

The 697 checkpoints correspond to jobs 1–696 and 706. Jobs 697–705 failed at startup with `No genotype .raw files found`; jobs 707–1125 have no local logs or checkpoints. Thus 428 expected checkpoints are absent. Batch 4's continuation was queued as 49047601, but the nine startup failures in batch 3 prevent its `afterok` dependency from succeeding. Its live scheduler state is not available here.

There are 278,800 saved replicate slots: 246,312 without a recorded error and 32,488 errors. This is 54.7% of the full design with no recorded error, not a certification of input integrity or convergence. Of the 697 checkpoint jobs, 550 have no recorded replicate errors and 147 have errors. Batch 3 alone has only 6,606 successful returns among 38,800 attempts (83.0% failed).

Counts were reconciled from all 706 simulation log pairs and checkpoint filenames. Ten jobs resumed earlier checkpoints, leaving 1,775 prior records outside the new progress logs. Those ten checkpoints were inspected, revealing one additional earlier error in job 305 beyond the 32,487 errors in the new logs.

### Input directory is unsuitable for simulations

The saved input path is `/project2/mstephens/wdenault/susie_mix/temp_plink`. The simulation lists and randomly samples `.raw` files from this directory on every replicate (`script/sim/sim_workhorse.R:45`). The workhorse writes gene extracts to the same directory and removes them on exit (`script/scan_tissue_attempt/workhorse.R:61,119`). This is strong evidence of interference between a simulation input pool and temporary workhorse outputs.

The new logs contain:

- 32,288 errors saying no `.raw` files were found.
- Five missing/unreadable-file errors, plus one additional such error in the resumed portion of job 305.
- 154 errors for insufficient donors.
- 32 errors saying `object 'target_gene' not found`.
- Six zero-genetic-variance errors, one insufficient-causal-SNP error, and one connection error.
- 236 truncated-read warnings (`Discarded single-line footer`). These warnings can occur on records that otherwise finish successfully.

The `target_gene` messages originate from QC helper error paths that reference a variable not passed to the helper (`workhorse_utils.R:216,289`). They obscure the underlying no-SNPs-after-QC condition.

The successful records in 18 selected checkpoints were inspected: all 5,608 non-error records contained the three method rows and the expected slider PIP, CS, causal-delta and fallback fields; all were marked converged. This was a targeted inspection, not a full convergence/structural audit of all saved records. Saved package versions in the first inspected checkpoint were susieR 0.14.2, susieSlide 0.2.0, data.table 1.18.6.1 and matrixStats 1.5.0.

Across those selected conditions, the same seed drew 4–15 different genotype files. The mutable input pool therefore also breaks the intended reuse of the same input locus across conditions. Within each successful replicate the three methods still share the generated dataset.

### Completion and retry semantics

`script/sim/run_job.R:72` catches replicate errors and saves them as records. At line 64, a checkpoint with 400 records is considered complete even if some are errors. Therefore shell success and `Finished simulation job` do not certify 400 successful datasets. Simply submitting again will skip these checkpoints rather than repair the failed records.

Use a fixed, validated genotype corpus outside the workhorse scratch directory before further simulation. Preserve the existing checkpoints. Recovery needs either explicit failed-replicate retry logic or a new output directory; because truncated reads and a changing locus pool also affect nominal successes, a clean rerun from the fixed corpus is preferable for final comparisons.

## New real-data workhorse (array 49040700)

All 298 task log pairs are present:

- Chunks 1–99 fail immediately because `run_chunk_001.R` through `run_chunk_099.R` are absent on the cluster.
- Chunks 100–298 report completion, attempting 12,330 distinct genes. This leaves 6,138 of the 18,468 genes unattempted by this run.
- Every one of the 199 chunks that ran shows susieSlide loading. The current workhorse fits and stores both `fit_slide` and `fit_slide_perm` (`workhorse.R:777,786,863,864`). No missing-susieSlide-package errors or explicit tissue-level errors were found.
- Completion still includes gene-level failures: 461 PLINK exports failed because all variants were excluded for the requested chromosome. Matching the output gene/extraction blocks gives 433 X genes, 16 Y genes and 12 mitochondrial genes. These require a decision about chromosome coverage and the underlying genotype input.
- Chunk 172 contains nonconvergence warnings at 100 iterations. Logs do not identify the affected saved fit reliably enough to certify convergence by method.

The filename failure is explained by `script/scan_tissue_attempt/write_job.R:27`, which writes unpadded filenames such as `run_chunk_1.R`, while `job/test1:18–20` requests three-digit filenames. These agree only from chunk 100 onward. Make the writer and launcher consistent before retrying chunks 1–99.

The local generated inputs also mix generations: both padded and unpadded drivers exist, while `data/temp_index` has 185 old gene lists and its first list contains 100 genes, whereas the new design uses 298 chunks of up to 62 genes. Regenerate/synchronize a consistent set of drivers and lists; do not use the local old lists blindly to rerun the new chunk IDs.

The local `results` directory is empty. Consequently, the gene RDS files, their error records and slider fits cannot be verified here. Old outputs, if present on the cluster, do not establish that the newly failed chunks contain slider results.

## Downstream gap

`script/analysis/generate_summary_results.R` still summarizes additive, mixed and weighted-mixed fits. It does not consume `fit_slide` or `fit_slide_perm`. The simulation plotting pipeline supports three methods, but the real-data summaries/plots still need slider integration.

## Recovery order

1. Prepare a stable simulation genotype corpus separate from workhorse temporary outputs.
2. Reconcile the workhorse filename convention and regenerate matching 62-gene lists/drivers; retry the missing chunks.
3. Resolve or explicitly exclude the X/Y/MT extraction failures and inspect chunk 172's nonconverged fits in the saved RDS files.
4. Rerun simulations with appropriate checkpoint handling; resubmission alone will retain saved errors. Confirm cluster status before restarting any array.
5. Audit the actual gene results and extend real-data summaries to include both slider fits.

No production scripts, existing results, or cluster jobs were changed by this audit. The helper `tmp/audit_slide_saved.R` can inspect selected checkpoint IDs supplied as arguments; without arguments it reads all locally present checkpoints.
