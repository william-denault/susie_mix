# Resume 005, expand at 006, continue through 040 in the same folder

Project: `/project2/mstephens/wdenault/susie_mix`.
Results remain in `results_slide_prior_em/iteration_005`, `iteration_006`, ...,
`iteration_040`. Target 040 means the folder numbered 040; 000 is the initial pass.

This bundle implements the approved deletion policy: keep the current and
previous full checkpoints, and preserve every iteration's priors, counts,
summaries and history. It does **not** keep an archive of deleted full results.
Do not manually remove iteration folders.

## Upload and start once

Upload `slider_inplace_rcc_2026-10-29.zip` to `/home/wdenault`, then:

```bash
cd "$HOME" &&
unzip slider_inplace_rcc_2026-10-29.zip &&
cd "$HOME/slider_inplace_2026-10-29" &&
sha256sum -c SHA256SUMS &&
bash start_005_to_040.sh
```

This is a combined bundle; no separate partial-recovery or rolling bundle
needs to be launched as well. If the earlier partial-recovery bundle has
already been launched, its receipt is reused and it is not launched twice.
Keep the extracted bundle and HOME state folders while the run is active.

1. If bulk compaction is still unfinished and no partial handoff exists,
   submit a held finalizer, cancel only the recorded compactor and its pending
   resume (original receipt: 49638385 and 49638386), then release the finalizer.
   No simulation jobs are cancelled. An already running old resume causes a
   protective stop rather than cancellation.
2. Finalize the interrupted compaction transactions and metadata, install the
   saved-call storage fix, and resume 005. Valid completed genes are reused;
   saved files without cached summaries receive parallel validation. No new
   bulk compaction pass is started.
3. A small `slide_next_gate` waits for the original recovery and slider jobs to
   exit. It verifies that 005 is complete before submitting setup. If 005 stops
   incomplete, the continuation stops with a diagnostic instead of skipping it.
4. Setup checks the 680 expression IDs (679 renamed symbols) against RCC
   expression rows, annotations and genotype windows. It prepares 006 and a
   private runtime, without replacing the original EM/simulation scripts.
5. The controller retires full results through 004, then fits 006. It continues
   through 040 with up to 298 workers, one gene per fit job and a refill check
   eligible every five minutes. All other user jobs count toward 300 slots.
   Slurm priority can delay when an eligible controller actually starts.

The initial 006 prior is exactly the ordinary update from completed 005.
Existing genes use 005 warm starts; added genes start cold in 006. Their counts
first enter the prior update for **007**. From then on the expanded gene/tissue
cohort is held fixed. Input availability does not guarantee an eligible fit
for every gene: the same sample/genotype QC and convergence checks still apply.
Expected manifest: 19,148 entries, including 851 excluded non-autosomal genes;
18,297 autosomal expression IDs are eligible for the scan.

## Retention sequence

| Before fitting | Full checkpoints kept | Full results retired |
| --- | --- | --- |
| 006, after validated 005 | 005 and 006 | 000â€“004 |
| 007, after validated 006 | 006 and 007 | 005 |
| 008, after validated 007 | 007 and 008 | 006 |
| 040, after validated 039 | 039 and 040 | through 038 |

Cleanup never fits a model or decompresses old full posterior files. It checks
the completed previous checkpoint and frozen next priors, records a durable
file plan, deletes only approved old `results/*.rds`, then repairs retained
validation stamps. An interruption resumes from the recorded plan. New workers
cannot start until cleanup finishes. Retired checkpoint certificates are
archived; their priors, counts, summaries, audits and folders remain.

Global prior history stays at `results_slide_prior_em/prior_history.csv`.
The original history through 005 is also preserved byte-for-byte under
`results_slide_prior_em/inplace_expansion_2026-10-29/prior_history.before_expansion.csv`.
Deleting old full results means they cannot later be used for detailed
posterior comparisons; 039 and 040 remain available when the run finishes.
The saved-call patch removes only the embedded `fit$call`, including from warm
starts before fitting. Posterior arrays remain, so storage is smaller but not zero.

## Monitor

```bash
squeue -u "$USER" -o "%.18i %.28j %.12T %R"
cat "$HOME/slider_inplace_launch_2026-10-29/last_submission.json"
tail -n 30 "$HOME/slider_partial_finalize_2026-10-29/logs/"*.out
tail -n 30 "$HOME/slider_inplace_launch_2026-10-29/logs/"*.out
tail -n 30 /project2/mstephens/wdenault/susie_mix/slider_inplace_runtime_2026-10-29/job/slide_inplace_control_*.out
du -h --max-depth=1 /project2/mstephens/wdenault/susie_mix/results_slide_prior_em | sort -h
```

Some log paths do not exist until the corresponding job starts. Worker logs
remain inside each iteration's `recovery/parallel/batches` directory. Cleanup
receipts are under `results_slide_prior_em/rolling_cleanup_inplace/before_NNN`.
Setup reports are under `results_slide_prior_em/inplace_expansion_2026-10-29`.
The small gate occupies one job slot until 005 is finished.

If the old shared controller still owns submitted simulations after 005, the
gate waits for that controller to finish as well; it does not interrupt its
simulation chain. Once started, the new controller uses available slots
alongside unrelated simulation jobs.

For a stopped continuation, inspect its errors and resolve the cause first.
When no continuation controller/workers remain, rerun:

```bash
cd "$HOME/slider_inplace_2026-10-29" && bash start_005_to_040.sh --resume
```

An ambiguous `sbatch` submission leaves a pending receipt and deliberately
blocks retries. Inspect Slurm and the recorded command before resolving it;
do not blindly delete receipts. A running continuation is never duplicated.
The original recovery target remains 005; only the private continuation targets 040.

## Validation

Local synthetic R integration used the real validators, prior update and
warm-start readers through 005â†’006â†’007â†’008, including interrupted deletion and
all 680 candidate IDs. Tests verified unchanged original 005 fits/priors/counts,
unchanged EM numerical helpers, repaired certificates, complete history, and
added counts first updating 007. Python tests cover deletion scope, changed
files, crash recovery, queue guards, finalized recovery state and the 300-job
submission budget. Shell tests mock Slurm; actual RCC fits and scheduling have
not been run locally. RCC setup repeats the input and package checks before 006.
