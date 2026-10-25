# Resume the tissue-wide slider-prior EM with parallel count collection

This is the continuation for the existing `results_slide_prior_em` run. It
adopts the saved `iteration_000/results/*.rds` and does not restart successful
gene fits. The old 23-hour serial validation step is replaced by parallel
validation tasks with per-gene checkpoints. Later iterations fit one gene per
Slurm task, retaining the existing full posterior warm starts.

## What was already implemented, and what changes

The previous scan already implemented the outer, tissue-specific EM. The new
driver changes execution, validation and collection, not its statistical update.
The old rolling launcher defaulted to four genes per task and refilled below
200 workers; its final audit and subsequent pooling were serial and repeated
when preparing the next iteration.

The separate n250/r300 simulations use a fixed uniform 17-point slider prior.
They benchmark the regional fitter, without learning shared weights across the
simulation replicates. The shared-scheduling update changes launchers and adds
a guard against regenerating the frozen simulation jobs. It contains no
numbered simulation scripts, simulation outputs, R package binaries or shared
`workhorse_utils.R`/annotation changes. No package reinstall is required.

The scientific update is Section 6, equation (36), of `slider_model-3.pdf`:

1. Fit every eligible gene/tissue with the iteration's frozen prior weights.
2. Sum the joint SNP/slider posterior mass over free SNPs and positive-variance
   components. Exclude forced-additive SNPs, explicit null candidates and
   exactly zero-variance components. Do not filter by credible sets or
   renormalize each retained component's free mass.
3. Normalize the 17 counts separately for each tissue. If a tissue has zero
   informative mass, retain its previous weights.
4. Refit the same frozen cohort with those weights. Final reported fits are
   under the last requested iteration's weights.

The same `spe_pool_summaries()` M-step is used by both the serial reference and
the parallel driver. The original fitting settings, X/Y/MT exclusions,
reviewed exclusions and fitting retry allowances remain in force.

## Enable sharing for the controller already submitted as 49452023

Upload **`slide_prior_shared_rcc_2026-10-23.zip`** to the project root on RCC.
This is a different archive from the first, wait-for-simulations bundle.
Job 49452023 was submitted successfully; Ctrl-C after the returned job ID did
not cancel it. Its five-minute continuation may now have a newer ID.
**Do not submit a second EM controller and do not cancel the simulations.**

```bash
cd /project2/mstephens/wdenault/susie_mix
unzip -o slide_prior_shared_rcc_2026-10-23.zip
sha256sum -c job/slide_prior_shared_SHA256SUMS
module load R/4.2.0
Rscript --vanilla script/scan_tissue_attempt/run_slide_prior_parallel.R preflight "$PWD" latest
python3 script/scan_tissue_attempt/slide_shared_scheduler.py enable "$PWD"
```

Run each command only if the preceding one succeeds (use `&&` when pasting a
single chain). Python 3.6 or newer is required, with no additional packages.
Activation finds the current controller and simulation submission records;
no job IDs need to be hard-coded. If a launcher is in its short RUNNING step,
activation stops without changing jobs; rerun the last command five minutes
later. Existing user/admin holds are respected.

The durable sharing flag overrides the existing controller's inherited
`SUSIE_SLIDE_WAIT_FOR_SIMULATIONS=1`. Slurm saves submitted scripts, so the
already queued old controller first submits one updated continuation. Sharing
normally begins within one or two checks, plus Slurm scheduling delay.

`latest` resumes the latest prepared iteration, initially 000. Target 20 remains
inherited from the submitted controller. The target is an absolute
iteration index: initialization 000, then 20 learned-prior updates, ending with
fits under iteration 020. It is a maximum requested stage, not an automatic
claim of statistical convergence. The larger numerical fitting budget
does not change the frozen convergence tolerance or prior model.

The node exclusions are specific to the failures already seen in this run;
they propagate to workers and controllers. They can be omitted on a future
launch if those nodes have been repaired.

## Scheduling alongside the submitted simulations

One controller now fills **up to 298 combined pending/running worker tasks**
within a 300-job budget, reserving the current and next controllers. All other
user jobs reduce available capacity. When both workloads have enough work,
new allocations move toward approximately 149 EM + 149 simulation workers.
Existing jobs are never preempted to achieve that split. If one workload cannot
use its allocation, the other can use it. A single freed slot is replenished.
Slurm controls actual start times; 298 submitted tasks do not guarantee 298
running tasks. Five minutes is eligibility time, with planning and queue delay
added. Other independent/manual submissions can consume slots after a snapshot.

Activation reads the actual last submitted simulation batch. If batch 1 is
still current, conditions 1-300 remain in that existing array and conditions
301-900 enter the shared queue. If batch 2 or 3 has already been submitted,
those arrays are also left intact and only the remaining conditions are adopted.
Future simulations execute the unchanged numbered R scripts, with the same
300-replicate checkpoints, inputs and seed settings. They remain simulations
with uniform slider weights; they do not consume the learned EM priors.

During migration, only the tiny pending EM/simulation launcher jobs are held
briefly, then released. Simulation worker jobs are never held or cancelled.
Original future batch launchers are backed up under
`results_slide_prior_em/shared_scheduler/original_launchers/`, then replaced
with single lightweight receipt jobs. This handles the old simulation launcher
already saved inside Slurm: it cannot accidentally submit another 300 workers.
New launchers see the sharing flag and exit without submitting duplicate work.
Two extra slots are reserved until those legacy launcher/receipt jobs drain.
If the existing queue already exceeds 300, no worker additions are made until
space becomes available; keeping the controller chain alive can temporarily
add one small continuation job to that inherited queue.

An interrupted migration is recorded in `shared_scheduler/migration.json`.
Rerun `enable` to finish it and release only the holds it owns. It does not
release unrelated holds or submit a second EM chain. Do not regenerate
`script/sim/jobs_slide` or replace simulation source files during this run;
the scheduler checks their frozen hashes before submitting further work.

Validation tasks handle up to 25 genes, one isolated R process at a time, and
checkpoint each gene immediately. Fitting tasks handle exactly one gene.
Workers retain the original conservative 48 GB / one CPU allocation; the new
summary-only coordinator requests 8 GB. No memory reduction is assumed for
the large existing posterior files.

Simulation workers retain 20 GB / one CPU / 23 hours and resume the existing
replicate checkpoints after interruptions, with at most three attempts per
adopted condition. The first bundle's standalone mode still waits for simulations
until sharing is explicitly enabled. Do not manually switch that old mode to
`WAIT_FOR_SIMULATIONS=0`; use the handover above.

## Progress and restart

### Extend an already running chain to 40 updates

The October 25 target-control patch adds `set-target` and a read-only metadata
`audit` mode to the R CLI. After installing that patch atomically, run from the
RCC project root:

```bash
module load R/4.2.0
Rscript --vanilla script/scan_tissue_attempt/run_slide_prior_parallel.R audit "$PWD" latest
Rscript --vanilla script/scan_tissue_attempt/run_slide_prior_parallel.R set-target "$PWD" latest 40
```

Chain these with `&&` to extend only if the audit succeeds. The target is saved
in `results_slide_prior_em/run_control.rds`. Every planning check reads this
file, overriding the old inherited `SUSIE_SLIDE_TARGET_ITERATION=20`. No worker,
controller or simulation needs cancelling or resubmitting. This is **40 total
learned-prior updates**, ending with fits under iteration 040, plus initial 000.
It is an iteration cap; the controller still stops for unresolved failures.
`set-target` only extends targets. It changes no priors, fits, fitting tolerance
or posterior calculations. Controller stderr prints `Run target iteration: 40`.

The metadata audit compares every prepared iteration's tissue weights with
its normalized pooled counts, checks history/snapshot agreement and the frozen
successful gene/tissue cohort, and samples compact result summaries. It does
not reopen every large posterior, and it does not establish biological validity
or convergence of the outer EM. Stale summaries can occur while workers write
new results; their count is reported separately. Monitor the root
`prior_history.csv`, especially per-tissue `prior_max_change` and `w_01..w_17`.

```bash
squeue --user="$USER" --array --format='%.20i %.28j %.10T %.12M %R'
cat results_slide_prior_em/parallel_continuation_job_id.txt
tail -n 30 job/slide_parallel_control_JOBID.out job/slide_parallel_control_JOBID.err
```

Controller messages report the iteration, resolved genes, active tasks and
available slots. R progress messages in `.err` are not themselves errors.
Validation can temporarily use one job for several genes; the one-gene task
rule applies to fitting. When validation detects an unresolved fit, the normal
bounded fitting retry policy applies. It never counts a failure as zero data.

Per iteration, new files are placed in:

```text
iteration_NNN/recovery/parallel/summaries/GENE.rds
iteration_NNN/recovery/parallel/index.rds
iteration_NNN/recovery/parallel/batches/
iteration_NNN/recovery/parallel/pooled.rds
iteration_NNN/recovery/parallel/COMPLETE.rds
iteration_NNN/slider_counts.csv
iteration_NNN/proposed_priors.csv
iteration_NNN/priors.csv
iteration_NNN/prior_history.csv
```

`proposed_priors.csv` contains the weights implied by that iteration's fits.
The weights actually used for those fits remain in its existing `priors.csv`.
The next prepared iteration records the update in `prior_history.csv` and its
own `priors.csv`. Every newly created iteration also gets a `prior_history.csv`
snapshot through that iteration, one row per tissue per iteration, containing
`w_01` through `w_17`, informative mass, maximum weight change and provenance.
For example, `iteration_001` contains weights learned from iteration 000 and
fits made using those weights; its history contains stages 000 and 001.
The root history continues accumulating later updates. A proposal at the final
target does not launch another fit.
Compare `prior_max_change` in the proposed priors to assess weight stability.

The completion certificate belongs to `recovery/parallel`, not the old serial
`recovery/COMPLETE.rds`. Once adopted, use `resume_slide_prior_em`; the old
`recover_slide_prior_em` launcher refuses to switch the chain back to serial
verification. The preparation lock is shared with the legacy launchers.

To stop future refills, cancel the specific controller ID printed in the last
controller log (also saved in `parallel_continuation_job_id.txt`). Let existing
workers finish. This stops refills for both workloads. Resume a genuinely
stopped chain only after confirming there is no live `slide_parallel_control`:

```bash
SUSIE_SLIDE_TARGET_ITERATION=20 SUSIE_SLIDE_WORKERS=298 \
SUSIE_SLIDE_POLL_MINUTES=5 SUSIE_SLIDE_MAX_ITER=5000 \
SUSIE_SLIDE_EXCLUDE_NODES=midway2-0003,midway2-0011 \
sbatch --exclude=midway2-0003,midway2-0011 --chdir="$PWD/job" job/resume_slide_prior_em latest
```

Completed summaries and fits are reused. Do not remove the gene results, frozen metadata or attempt
receipts. Validation interruptions have their own default three-attempt budget
and do not consume/reset fitting retries. Exhausted retries or workers that
never entered R stop for log inspection.

If submission was interrupted before its returned job ID could be recorded,
`parallel_pending_submission.txt` blocks further launches. Inspect Slurm and
the recorded batch before repairing that registration; this prevents duplicate
jobs after an ambiguous submission.

## Validation and integrity

Workers retain the existing full-posterior validity checks and save counts in
the same pass. Each summary contains the result checksum, file size/mtime,
frozen prior/settings/manifest identities, exclusion policy and tissue cohort.
The coordinator checks source and summary file size/mtime and only reads new
small summaries, maintaining a compact index. The M-step checks all required
genes and the fixed gene/tissue cohort, and does not re-open full posteriors.

This relies on completed result files being immutable under the gene locks.
Normal replacements change size/mtime and trigger recollection. External edits
that deliberately preserve both size and timestamp are outside that contract;
they require a new distributed validation, not trust in the old summaries.

Aggregate ELBO reporting is off by default (`SUSIE_SLIDE_RECORD_ELBO=0`). Saved
finite fit objectives and convergence checks remain part of fit validity; they
do not require a separate pass through results. To record the optional objective
history, enable `SUSIE_SLIDE_RECORD_ELBO=1` consistently from iteration 000.

## Local validation shipped with the bundle

The fixture tests cover checkpointed adoption, interrupted validation, exact
serial/parallel count and prior equivalence, forced/null/zero-variance exclusions,
zero-information tissues, unchanged initialization files, stale/corrupt results,
frozen cohorts, duplicate jobs and fitting retry preservation. A real-fitter
smoke test runs two updates with warm starts. Mock Slurm tests exercise slot
refills, job budgets and submission failures. The shared-scheduler fixtures
cover migration with existing queued jobs, inherited wait settings, balanced
allocation, single-slot refill, interrupted handover, checkpoints and retries.

```bash
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_parallel.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_parallel_smoke.R
bash script/scan_tissue_attempt/tests/test_slide_prior_parallel_launcher.sh
python3 script/scan_tissue_attempt/tests/test_slide_shared_scheduler.py
```

Tests use temporary fixtures, never the production fits. Full RCC performance
and filesystem throughput have not been benchmarked locally.
