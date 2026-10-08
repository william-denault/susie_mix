# Recover slider-prior EM with bounded batches

This scheduler changes execution and recovery, not the statistical model. It
keeps the frozen grid, priors, fitting settings, package version, input QC,
warm starts and `spe_pool()` M-step. PLINK uses one thread to match the one-CPU
allocation. New tasks handle at most **4 genes by default**, each in a **fresh
Rscript process**, and save each gene before starting the next. Set
`SUSIE_SLIDE_GENES_PER_TASK=5` for a five-gene cap; values from 1 to 15 are
supported. The coordinator exports this setting to its continuations.

The original audit manifest retains its 1,847 groups of at most 10 genes.
Planning subdivides each group's unresolved genes into smaller execution tasks
before applying the available job capacity. Existing manifests, batches, fits
and retry counts are preserved. Changing this execution cap needs no new audit
and does not alter already-submitted tasks. Each new array has at most 298 tasks.

The first recovery target is the existing `iteration_000`. Do not delete or
replace `results_slide_prior_em`, its gene results, or its original manifest.
The old completion markers are retained as historical records and are not
trusted by recovery. No mixed-model files are modified.

## Upload and audit first

Extract `output/slide_prior_recovery_rcc.zip` into the RCC project, preserving
the paths. The archive contains code, documentation, tests and the 36-gene
pilot list at `output/slide_prior_pilot_genes.txt`; it contains no input datasets
or fitted results. It includes updated guards for the historical slider launcher
and preparation script; upload these too.

First verify the old slider array and continuation have stopped. Allow the
mixed-model run to reach your chosen stopping point before the full recovery:
its launcher independently requests up to 298 tasks. The recovery coordinator
counts other queued/running jobs, but cannot reserve slots against a different
launcher that submits new jobs after that count.

On RCC:

```bash
cd /project2/mstephens/wdenault/susie_mix/job
sbatch recover_slide_prior_em audit 0
```

This is an audit only. It reads the actual RCC gene files, including genes
inside the five old chunks marked complete, and writes:

```text
results_slide_prior_em/iteration_000/recovery/
  state.rds                  # hashes of the frozen original files and plan
  manifest.csv               # task, gene, original_chunk; <=10 genes per task
  audit.csv                  # status, detail, result checksum, number of fits
```

The local copy without `results/` cannot establish which genes are reusable;
run this audit on RCC. Progress is logged every 100 genes. A full audit can
take substantial time because it reads all saved posteriors. Auditing does
not run the fitter or submit workers.

Statuses:

| Status | Meaning |
|---|---|
| `valid` | Converged fits with the correct iteration, priors and posterior structure |
| `qc_excluded` | An existing result has no eligible tissues and no recorded errors, under the existing input QC |
| `retry` | Missing, unreadable, unconverged or invalid result, or an unresolved gene/tissue error |
| `reviewed_exclusion` | An exact saved error record explicitly reviewed for exclusion, or that frozen exclusion inherited in a later iteration |
| `scope_excluded` | Gene on X, Y or MT, excluded by the explicitly selected chromosome scope |

Old `.done` markers do not override any of these checks. Missing files are
retried, and successful genes in interrupted chunks are reused. Partially
successful genes with unresolved tissue errors are retried as whole genes.

## Run a pilot, then resume

### Exclude X/Y/MT for this analysis

The user selected this scope on October 7, 2026. After uploading the updated
bundle, activate it once, while no slider workers/continuations are active:

```bash
cd /project2/mstephens/wdenault/susie_mix/job
sbatch recover_slide_prior_em exclude-sex-mt 0
```

This reads the same GTEx gene annotation as the worker and freezes the first
chromosome assignment for every original gene, normalizing `chrX`, `chrY`,
`chrM`/`MT` and numeric aliases. An unmapped gene stops activation for inspection.
It records the mapping and annotation checksum in
`results_slide_prior_em/chromosome_scope.rds`, writes
`results_slide_prior_em/chromosome_exclusions.csv`, and updates the existing
recovery audit without rereading all fitted posteriors. The pre-scope audit is
preserved as `recovery/audit_before_chromosome_scope.csv`.

X/Y/MT genes are skipped even when no result exists, a retry budget was exhausted,
or an old fit succeeded. Their existing files remain intact. They contribute
no counts to the prior update and stay excluded in subsequent iterations.
Autosomal errors still require recovery; the fitter and `spe_pool()` calculation
are unchanged. Original manifests and frozen fitting settings are retained.
Old batches/certificates cannot be used across this scope change. Activation
can be safely repeated after interruption; changing scope after learned-prior
updates have begun is refused. No chromosome-inventory command is needed.

Wait for this small setup job to finish before the pilot. The updated
`output/slide_prior_pilot_genes.txt` contains 36 candidates after removing
the diagnostic X-chromosome case; the scope policy will also filter any other
excluded candidates identified by the annotation.
The previously synced pilot list also works: its X-chromosome gene is filtered
out automatically after scope activation.

### Pilot

For the supplied 36-gene diagnostic pilot, verify that the bundled list exists
on RCC, then submit from the project's `job` directory. With a four-gene cap,
these candidates occupy at most 12 tasks across their original audit groups:

```bash
test -s ../output/slide_prior_pilot_genes.txt && \
SUSIE_SLIDE_PILOT=1 \
SUSIE_SLIDE_PILOT_GENES=/project2/mstephens/wdenault/susie_mix/output/slide_prior_pilot_genes.txt \
SUSIE_SLIDE_GENES_PER_TASK=4 \
SUSIE_SLIDE_MAX_ITER=5000 \
sbatch recover_slide_prior_em run 0 12
```

If the list is missing, extract the updated bundle or sync that file before
submission. A controller that cannot read this list exits before planning or
launching any gene workers. The 5000 setting raises the inner solver's iteration
budget; it does not request 5000 EM updates.

Alternatively, after the audit job ends, start a **single five-task pilot**
(at most 20 genes at the default cap):

```bash
SUSIE_SLIDE_PILOT=1 sbatch recover_slide_prior_em run 0 5
```

To target known slow or failed genes, create a plain-text list with one gene
per line, selected from `audit.csv`, and provide its absolute RCC path:

```bash
SUSIE_SLIDE_PILOT=1 \
SUSIE_SLIDE_PILOT_GENES=/absolute/path/pilot_genes.txt \
sbatch recover_slide_prior_em run 0 5
```

Only unresolved listed genes are scheduled; the final command argument caps
the number of execution tasks, not genes. Include prior memory-failure and slow cases, not only
easy genes. Already validated genes are reused without refitting. The pilot
does not schedule a continuation or change any priors. Inspect its per-gene
logs, Slurm elapsed times and peak memory before scaling up. Fresh processes
limit cross-gene memory retention; they do not guarantee that a single large
gene fits in 40 GB or 23 hours.

Then, after the pilot is no longer active:

```bash
SUSIE_SLIDE_PILOT=0 SUSIE_SLIDE_GENES_PER_TASK=4 SUSIE_SLIDE_MAX_ITER=5000 \
sbatch recover_slide_prior_em resume 0
```

This retains the pilot's increased inner-solver budget. Use a five-gene cap
instead by changing only `SUSIE_SLIDE_GENES_PER_TASK=5`. Do not edit
`recovery/manifest.csv` or rerun the audit with a different manifest size.
The October 8 pilot snapshot had 35 validated genes (945 tissue fits) and one
unfinished gene, ACSL1. The slow task spent 16.3 hours on five preceding genes;
its first four alone took 14.3 hours. A smaller cap reduces accumulated runtime
but cannot guarantee that every future task will fit within 23 hours.

The default batch has at most 298 worker tasks. The coordinator subtracts
other submitted jobs and reserves two slots for the current and next control
jobs, keeping its submission within a 300-job budget for your user at that snapshot. It
counts **pending and running array elements**, not just parent array IDs.
Lower RCC QOS limits still apply; it does not bypass them. If no slot is
available it stops with a manual resume command. Submission failures preserve
all files and are reported; an ambiguous response requires queue inspection
before retrying. Duplicate launches are blocked by the shared preparation
lock, recorded job IDs and recovery job names.

Each batch queues one inspection coordinator with `afterany`, so a worker
failure does not discard the opportunity to recover. Batches run sequentially;
the next batch waits for the slowest task in the current batch. This first
implementation favors simple, auditable recovery over a rolling worker pool.

## Errors and bounded retries

An actual fitting attempt is recorded before it starts. Default maximum:
three attempts per gene across the recovery, including interrupted attempts.
Unattempted genes are prioritized ahead of retries. Successful checkpoints
are skipped. A worker error does not prevent the task from attempting its
other genes, but the task exits unsuccessfully if any gene failed.

The coordinator stops for inspection when the retry limit is reached. An
automatic continuation also stops when two genes report memory/system
allocation failures, when one gene reports such a failure twice, or if no worker in
the previous batch reached fitting; this prevents endless resubmission after
module/package/launch failures. A manual resume can retry a launch after its
cause has been fixed; manual resume also retries resource failures within the
remaining budget. Gene retry budgets are not silently reset. An exhausted
gene requires diagnosis and an explicit repair/review, rather than repeated
blind submission.

Errors are never automatically converted to biological exclusions. If an
initialization error is legitimately explained by input eligibility, review
the saved record and create `recovery/reviewed_exclusions.csv` with columns
`gene,result_md5,reason`. Use the checksum from `audit.csv` and document the
actual reason. This authorizes only that exact saved error record; a changed
record invalidates the approval. A partial gene retains its successful tissue
fits and its explicitly reviewed tissue-error exclusions. Do not approve
timeouts, allocation failures, missing packages or other infrastructure
failures as input exclusions. Ordinary zero-eligible-tissue results require
no manual approval. Keep this file empty/absent when there are no reviewed
exclusions. Re-run `audit` after adding reviews, then `resume`.

`SUSIE_SLIDE_MAX_ITER` remains available for a diagnosed inner numerical
iteration limit, as in the original workflow. Raising it does not cure a
wall-clock limit or memory failure and does not reset the retry counter.

## Completion and the next EM update

When all genes appear resolved, the coordinator rereads **every full result**
and writes `recovery/COMPLETE.rds` only if validation succeeds and at least one
tissue fit exists. The completion certificate is tied to the frozen files and
audit. A prior update revalidates all results again before calling the original
pooling/M-step. Original chunk markers are neither rewritten nor reused for
this decision. Later iterations must preserve the validated initialization
cohort, including reviewed exclusions.

**The recovery stops after finishing the requested iteration.** It does not
automatically launch 20 updates. Once initialization is complete and reviewed,
prepare and fit exactly one new update with:

```bash
sbatch recover_slide_prior_em new
```

This invokes the existing prior-estimation function, then audits and schedules
the new iteration using the same configurable execution cap. The original
`em_susie_slide_prior` launcher refuses to run once recovery is adopted, to
prevent an accidental return to 61-62-gene jobs. Use the new launcher for
subsequent updates and resumes.

## Inspect or stop a recovery

```text
recovery/audit.csv                      # latest audit/receipt status
recovery/attempts/<gene>.rds             # attempt count, timestamps, errors/warnings
recovery/blocked.csv                    # genes needing inspection, when applicable
recovery/batches/batch_*/manifest.csv    # exact genes in each submitted batch
recovery/batches/batch_*/worker_*.out    # array-task progress
recovery/batches/batch_*/gene_logs/*     # separate R stdout/stderr for each gene
job/slide_recover_control_*.out/.err     # audit/coordinator logs
```

The root `last_array_job_id.txt` and `last_continuation_job_id.txt` track this
workflow too. Read and verify their IDs in `squeue`, then cancel the queued
coordinator first and its array second. Cancelled work keeps its saved gene
results. A pilot has no new continuation; its recorded continuation ID may
belong to an older, completed job. Do not launch competing workers manually.

## Local verification

```bash
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_recovery.R
bash script/scan_tissue_attempt/tests/test_slide_prior_recovery_launcher.sh
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_em.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_em_worker.R
bash script/scan_tissue_attempt/tests/test_slide_prior_em_launcher.sh
```

The recovery tests use synthetic saved fits and mocked Slurm commands. The
real-worker test uses synthetic GTEx-format inputs and the installed discrete
slider package. Only an RCC pilot can verify real workload memory and timing.

The scheduler follows [Slurm's documented array accounting and dependency
semantics](https://slurm.schedmd.com/job_array.html): array tasks count against
job limits individually, and `afterany` on an array waits for all its tasks to
terminate. Rebuild the upload archive from the repository root with:

```bash
python script/scan_tissue_attempt/build_slide_prior_recovery_bundle.py
```
