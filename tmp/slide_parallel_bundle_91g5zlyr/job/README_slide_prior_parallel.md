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
simulation replicates. This bundle contains no simulation scripts, simulation
outputs, R package binaries or shared `workhorse_utils.R`/annotation changes.
It does not require reinstalling the package used by the running simulations.

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

## Upload and resume

Upload `slide_prior_parallel_rcc_2026-10-23.zip` to the project root on RCC.
The recorded old coordinator 49435779 timed out. If another *slider EM*
controller or worker has since been submitted, let it finish or stop that
specific chain before replacing the EM code. Leave the simulation jobs running.
The new launcher refuses to compete with a live legacy slider chain.

```bash
cd /project2/mstephens/wdenault/susie_mix
unzip -o slide_prior_parallel_rcc_2026-10-23.zip
sha256sum -c job/slide_prior_parallel_SHA256SUMS
module load R/4.2.0
Rscript --vanilla script/scan_tissue_attempt/run_slide_prior_parallel.R preflight "$PWD" latest

cd job
SUSIE_SLIDE_TARGET_ITERATION=20 \
SUSIE_SLIDE_WORKERS=298 \
SUSIE_SLIDE_POLL_MINUTES=5 \
SUSIE_SLIDE_WAIT_FOR_SIMULATIONS=1 \
SUSIE_SLIDE_MAX_ITER=5000 \
SUSIE_SLIDE_EXCLUDE_NODES=midway2-0003,midway2-0011 \
sbatch --exclude=midway2-0003,midway2-0011 resume_slide_prior_em latest
```

Run the submission only after preflight succeeds. `latest` resumes the latest
prepared iteration, initially 000. It also recovers a stopped chain after 001
or a later iteration has already been prepared. The target 20 is an absolute
iteration index: initialization 000, then 20 learned-prior updates, ending with
fits under iteration 020. It is a maximum requested stage, not an automatic
claim of statistical convergence. Set a smaller target before the first
submission if fewer updates are desired. The larger numerical fitting budget
does not change the frozen convergence tolerance or prior model.

The node exclusions are specific to the failures already seen in this run;
they propagate to workers and controllers. They can be omitted on a future
launch if those nodes have been repaired.

## Scheduling alongside the submitted simulations

The simulation chain currently submits a full 300-task array at each batch
boundary. By default the EM controller waits while any `susie_slide_sim` or
`susie_slide_launch` jobs are active, including pending dependency launchers.
It wakes about every five minutes and starts EM work once that chain is gone.
This prevents EM from taking the slots needed by the next simulation batch.
Nothing cancels, resubmits or edits simulations.

After that, the driver targets **298 pending plus running EM worker tasks**.
It counts all other user jobs and reserves two controller slots within its
conservative 300-job submission budget. It tops up every available slot, even
when just one worker has finished. Slurm still controls actual start times,
so 298 submitted tasks does not guarantee 298 running tasks. Five minutes is
the next controller's eligibility time; scheduler delay and planning add time.

Validation tasks handle up to 25 genes, one isolated R process at a time, and
checkpoint each gene immediately. Fitting tasks handle exactly one gene.
Workers retain the original conservative 48 GB / one CPU allocation; the new
summary-only coordinator requests 8 GB. No memory reduction is assumed for
the large existing posterior files.

`SUSIE_SLIDE_WAIT_FOR_SIMULATIONS=0` permits sharing available slots, but the
current simulation launcher does not reserve its future 300-task batch. Keep
the default for the already-submitted simulation chain.

## Progress and restart

```bash
squeue --user="$USER" --array --format='%.20i %.28j %.10T %.12M %R'
cat ../results_slide_prior_em/parallel_continuation_job_id.txt
tail -n 30 slide_parallel_control_JOBID.out slide_parallel_control_JOBID.err
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
```

`proposed_priors.csv` contains the weights implied by that iteration's fits.
The weights actually used for those fits remain in its existing `priors.csv`.
The next prepared iteration records the update in `prior_history.csv` and its
own `priors.csv`. A proposal at the final target does not launch another fit.
Compare `prior_max_change` in the proposed priors to assess weight stability.

The completion certificate belongs to `recovery/parallel`, not the old serial
`recovery/COMPLETE.rds`. Once adopted, use `resume_slide_prior_em`; the old
`recover_slide_prior_em` launcher refuses to switch the chain back to serial
verification. The preparation lock is shared with the legacy launchers.

To stop future refills, cancel the specific controller ID printed in the last
controller log (also saved in `parallel_continuation_job_id.txt`). Let existing
workers finish. Resume using the same submission above; completed summaries
and fits are reused. Do not remove the gene results, frozen metadata or attempt
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
refills, simulation priority, job budgets and submission failures.

```bash
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_parallel.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_parallel_smoke.R
bash script/scan_tissue_attempt/tests/test_slide_prior_parallel_launcher.sh
```

Tests use temporary fixtures, never the production fits. Full RCC performance
and filesystem throughput have not been benchmarked locally.
