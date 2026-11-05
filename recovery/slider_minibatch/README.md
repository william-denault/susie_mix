# Separate slider minibatch EM experiment — 2026-11-04

This bundle creates a separate branch of the expanded GTEx slider analysis. It
does not install anything in, cancel jobs from, or modify the current full EM.
It reads that run's latest retained, certified completed checkpoint when its
first controller runs. An incomplete newer iteration is not used.

## Run on RCC

Upload `slider_minibatch_rcc_2026-11-04.zip` to `/home/wdenault`, then:

```bash
cd /home/wdenault &&
unzip slider_minibatch_rcc_2026-11-04.zip &&
cd /home/wdenault/slider_minibatch_2026-11-04 &&
sha256sum -c SHA256SUMS &&
bash start.sh
```

The source is `/project2/mstephens/wdenault/susie_mix`. All new code, results,
logs and state are stored in its sibling:

```
/project2/mstephens/wdenault/susie_mix_slider_minibatch_2026-11-04
```

The launcher can wait up to ten minutes for a safe free submission slot. The
first controller has a five-minute delayed start. Existing jobs keep running.
If no safe slot becomes available, inspect the queue and run `bash start.sh
--resume` later; the already-created configuration is retained.

For status:

```bash
cd /home/wdenault/slider_minibatch_2026-11-04
bash status.sh
squeue --user=wdenault -o "%.18i %.28j %.12T %.12l %R"
```

`slide_mb_c_...` are new branch controllers; `slide_mb_f_...` are its gene
workers. The original `slide_inplace_...` jobs are the full EM. Branch controller
logs live under the new folder's `logs/`; worker logs under `tasks/task_*/`.

After an interruption, inspect `ERROR.json`, logs, and the queue, then resume:

```bash
cd /home/wdenault/slider_minibatch_2026-11-04
bash start.sh --resume
```

Do not restart while a branch controller is already queued: the launcher will
refuse duplicates. Valid completed gene results are adopted without fitting
again. A failed gene gets at most three submitted attempts per batch; repeated
failures stop the branch rather than silently removing the gene. An ambiguous
`sbatch` response leaves `pending_submission.json` and stops automatic
resubmission: inspect that command, the queue and accounting before resolving it.

## Exact statistical schedule

The required expanded source cohort is 18,263 retained autosomal gene IDs from
the 19,114-ID manifest, including the 646 added IDs. The 34 deferred mappings
remain excluded. Source iteration 007 or later must have a valid completion
certificate, cached per-gene summaries and retained results. The source's own
validator and frozen code hashes are checked; the scientific fitting function,
model settings, count calculation and convergence validator must match this
bundle. The existing renamed-gene map must match exactly.

For tissue t and slider grid value k, let C[t,k] be the source checkpoint's
pooled expected counts and N=18,263. Initialization uses:

```
S[t,k] = C[t,k] / N
prior[t,k] = S[t,k] / sum_k S[t,k]
```

Thus the starting prior is the ordinary next M-step from the completed source
fits. The source checkpoint is not refitted just to initialize the branch.

There are six complete passes. In each pass every gene appears exactly once
in nine randomly shuffled batches, balanced by gene count: two batches of 2,030
and seven of 2,029. The saved shuffle is deterministic (seed 20261104 plus pass)
and changes between passes. It does not select batches by runtime or credible
set count. These batches are not stratified by tissue or chromosome.

Each batch uses one frozen set of tissue priors. Every gene runs the unchanged
fitter, warm-started from that gene's latest validated posterior. Once every
gene has finished and been validated, its expected counts are pooled into B:

```
S_new[t,k] = 0.7 * S_old[t,k] + 0.3 * B[t,k] / batch_gene_count
prior_new[t,k] = S_new[t,k] / sum_k S_new[t,k]
```

This averages sufficient statistics per sampled gene. It does not average
individually normalized gene priors or change the weighting according to which
jobs finish first. Zero-count tissue rows retain the preceding prior. The same
active-component/free-SNP rules apply as in the full EM: no credible-set filter
and no added pseudocounts. Fixed eligibility/exclusion metadata is preserved.

After 54 minibatch updates, five ordinary full-data E/M steps run on all 18,263
genes. Each full step uses one frozen prior and the established full-data M-step
and Q-increase check. The first full step warm-starts from the completed sixth
pass; subsequent full steps warm-start from the preceding full step. There are
59 committed updates in total. The schedule stops after five full steps even
if convergence has not yet been reached; inspect prior changes to assess it.

The final fitted results use `full_em/iteration_005/priors.csv`. The M-step
computed from those results is `full_em/iteration_005/prior_after.csv` (also
`latest_priors.csv`). These are intentionally distinct: no extra sixth full
E-step is silently added to refit under the last proposed prior.

## Storage and prior history

The branch pins source results with hard links, not symlinks or a bulk copy.
The original cleanup can unlink its own filenames without breaking branch
warm starts. Source and branch must be sibling directories on the same
filesystem; a failed hard link stops setup. Pinned files are never rewritten.
Their recorded checksums are verified by the gene workers before first use.
If the original cleanup is refreshing its certificate or removing a checkpoint
during capture, the already-queued successor retries capture from a newer
certified checkpoint. Scientific-code or cohort mismatches stop setup.

Hard links initially consume negligible additional result storage. They keep
the source checkpoint's underlying bytes allocated until the branch releases
its pins, even if the original run has already retired that checkpoint. Pins
are released at the start of pass 2, after all genes have branch warm starts.

Large results are kept for only the current and previous pass/full-step groups.
For example, pass 3 retires pass 1; full step 2 retires pass 6. At completion,
full steps 004 and 005 retain results. Cleanup removes only `.rds` files in the
new branch's old `results` directories. Prior files, per-gene summaries, pooled
counts, validation receipts, audit tables, commits and history remain. The
original EM's cleanup policy is unchanged.

Result objects use the existing storage fix: remove the captured R `$call`,
which could otherwise contain nested previous fits and input matrices, while
retaining the fitted posterior needed for warm starts. They are saved with the
existing RDS writer; there is no serial bulk recompression pass.

Useful files in the new branch:

- `seed/READY.json`: chosen source iteration and gene/tissue-fit totals.
- `seed/snapshot_*/source_prior_history.csv`: preserved original prior history.
- `prior_history.csv`: starting prior plus every committed new tissue update;
  columns identify stage, phase, pass, batch, rho and maximum prior change.
- `latest_priors.csv`: latest committed proposed tissue priors.
- `stages.csv` and `plan/`: complete frozen 59-stage schedule and assignments.
- `passes/pass_001/batch_01/`, etc.: minibatch outputs and histories.
- `full_em/iteration_001/` through `iteration_005/`: full-step outputs.
- Each stage's `priors.csv`, `prior_after.csv`, `counts.csv`,
  `running_counts_per_gene.rds`, `audit.csv`, and `history.csv`.
- `DONE.json`: all six passes and five full steps committed successfully.

## Queue sharing and recovery

The new branch has a default maximum of 148 submitted/running gene workers,
one gene per array task, with 48 GB and 23 hours per worker. Controllers request
16 GB and 23 hours. The limit is a cap, not a guaranteed equal allocation; the
existing full EM may already occupy most slots. All user jobs, including pending
array elements, count toward the shared 300-job budget. Cluster resources and
priority still determine how many can run at once.
One extra slot is reserved for each external controller's not-yet-submitted
successor. With the original controller present, the new branch therefore
fills at most 299 total occupied slots; the remaining slot lets the original
queue its next controller without stopping.

The current full-EM controller uses a nonblocking prepare lock; taking its lock
from the new branch could make it exit. This bundle leaves that lock untouched.
It submits only while external controllers are absent or pending specifically
for `BeginTime`, at least 90 seconds before their earliest start. Each sbatch
RPC has a ten-second timeout. Pending jobs with `Priority`/`Resources` do not
provide a safe submission window. Existing controllers already account for
unrelated user jobs when calculating their own available slots. Concurrent
manual submissions can still consume slots; avoid launching another independent
bulk submission at the same instant.

A successor controller is queued before potentially lengthy capture or pooling
work, with `afterany` and a five-minute delay. A timeout can resume at the last
durable boundary. A recorded scientific/validation error stops the chain and
preserves the original error. No prior is updated from a partially completed
batch. Finalization writes a commit payload before advancing the state pointer,
so recovery does not apply rho twice. Refill is checked approximately every five
minutes, subject to queue availability and safe submission windows.

Use `bash start.sh --workers 100` to choose a lower cap at first launch. The cap
is frozen in the experiment configuration; use the same option when resuming
that configuration. Default launches/resumes need no option. An independent
branch uses its own target and ignores the original EM's target-iteration env
variable. No convergence-speed improvement or equal wall-clock sharing is
guaranteed by the chosen stochastic schedule.

## Validation performed before packaging

- Eleven scheduler tests: all 18,263 IDs covered once per pass, deterministic
  reshuffling, 59-stage schedule, job-array expansion, safe submission windows,
  global capacity reservation including both controllers' successors, immutable completion
  receipts, ambiguous submission handling and Python 3.6-compatible syntax.
- Real hard-link/retirement tests: original unlink leaves the pinned file usable;
  cleanup refuses outside-branch targets and preserves the preceding group and
  all history.
- Complete synthetic 59-stage execution with exact count-update comparisons,
  renamed-ID dispatch, warm starts, null/forced exclusions, zero-mass fallback,
  failed-convergence/prior/cohort checks, interrupted result adoption, and
  exactly-once commit recovery. Source fixtures remained unchanged.
- Source capture against a previously certified expanded-cohort fixture using
  its actual private runtime and completion validator.
- Tiny real `susieRSlidePrior` fits through the new worker, including warm starts
  and both minibatch and full-data M-steps.

No RCC jobs or production GTEx fits were launched during local validation.
The production checkpoint/cohort/package checks run on RCC during setup.
The included R tests use local authoring fixture paths; they are development
tests, not commands that need to run on the RCC login node.
