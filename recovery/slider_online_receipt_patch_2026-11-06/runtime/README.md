# Slider recovery: finish 014, online 015-020, full EM 021-025

This bundle resumes the same project and results folder:
`/project2/mstephens/wdenault/susie_mix/results_slide_prior_em`.
It preserves existing 014 results and all running/pending gene workers.
It requires iteration 014 to be prepared and 015 and later folders to be absent.
The old failed `slider_online_014_024` installation is left unchanged.

## Launch on RCC

Upload `slider_online_015_rcc_2026-11-05.zip` to `/home/wdenault`, then run:

```bash
cd /home/wdenault &&
unzip slider_online_015_rcc_2026-11-05.zip &&
cd /home/wdenault/slider_online_015_2026-11-05 &&
sha256sum -c SHA256SUMS &&
bash start_detached.sh
```

`start_detached.sh` checks the bundle, starts a nohup launcher, and returns a
local process ID. This is not a Slurm job ID or proof that the handoff succeeded.
It automatically supplies `--resume` if this exact workflow already exists.
The terminal may be closed; the launcher can wait up to 24 hours for a pending
legacy controller. Only one launcher can own the workflow lock at a time.

```bash
tail -n 40 /home/wdenault/slider_online_015_2026-11-05/handoff_launch.log
```

Successful launch prints `Handoff secured` and `Submitted controller JOBID`.
The persistent receipt is:
`PROJECT/workflows/slider_online_015_025/HANDOFF.json`.
Inspect `last_controller.json` for the current Slurm controller ID. A genuine
error or an ambiguous submission stops; repeated launches do not erase errors.
Do not use the old 014 bundle or old full-EM launcher to resume this workflow.

## Why the earlier handoff failed, and what changed

The old code requested `squeue %R`, which prints pending reasons with
parentheses, e.g. `(BeginTime)`, but compared that result to `BeginTime`.
Consequently the safe-window check rejected pending controllers, and the
overnight retries kept encountering the same error. The original full EM
continued into 014, which the online launcher correctly refused to overwrite.
The original tests used unparenthesized reasons and missed this defect.

This recovery requests raw `%r` reasons and also normalizes parentheses.
A regression test uses the real Slurm output shape. The handoff additionally
uses an explicit verified hold, so it no longer depends on an estimated start
time or a narrow 90-second submission window. Original cleanup before 014
must already be certified before any controller is held.

The new launcher verifies the pending controller's job ID, user, name,
project command and recorded continuation ID. It places that one controller
on hold with `scontrol hold`, confirms it is still PENDING with Priority=0,
and then cancels only that held controller using exact ID/name/user/state
filters. A controller which starts during the initial hold race is released
and allowed to finish; it is never cancelled. Every hold/cancel operation is
journaled in `handoff_pending.json` before the operation. If interrupted after
a hold, resume this bundle; the held job and recorded state are inspected again.
Running and pending gene workers are not held or cancelled.

One new controller now handles both finishing 014 and the online run. While
finishing 014 it invokes the unchanged private `run_slide_prior_parallel.R`
planner with TARGET_ITERATION=14 and reuses its exact gene worker script.
Saved genes and existing active tasks are recognized by that planner. No
ordinary full iteration 015 is prepared. Once all tasks resolve, the original
completion validator certifies 014 before online seed capture is permitted.
The new controller records this transition in `SOURCE_014_COMPLETE.json`.

The original private runtime is verified against `frozen.json`, not edited.
No changes are made to the statistical fitter, count rules or validation.

## Schedule, history and storage

- Checkpoint: completed, certified full EM 014, including its learned counts.
- Online passes: 015-020, nine shuffled balanced batches per pass, rho=0.3.
- Each pass uses all 18,263 active gene IDs exactly once: two batches of 2,030
  and seven of 2,029, fixed seed 20261104 plus pass.
- Full EM: 021-025, five steps, then stop. Convergence is assessed from history.
- 298-worker maximum and a 300-job all-user budget, including controllers and
  unrelated queued jobs. A successor is queued before long metadata work.
- Refill checks are eligible every ONE minute; Slurm may delay their start.
  This changes scheduling, not when prior updates happen.
- Frozen batch prior; every batch gene must resolve before its update.
- Scientific update: S_new = 0.7*S_old + 0.3*batch_counts/batch_gene_count;
  normalize S within each tissue to get its 17 slider probabilities.
- S starts at certified 014 counts divided by 18,263; the initial 015 prior
  is the ordinary full-data update from 014. Full steps use the original
  pooled-count M-step. No per-gene normalized-prior averaging is introduced.

Results stay in `iteration_XXX/results/GENE.rds`. Online batch metadata stays
under `iteration_XXX/batches/batch_01` through `batch_09`; full steps use
`batches/full_step`. Each online-owned iteration has `ONLINE_OWNER.json`.

Root and iteration-specific `prior_history.csv` preserve the original columns
and values and add `iteration`, `batch`, phase, input/after_update role, rho,
stage and all slider grid values. Online input uses batch 0; updates use 1-9.
All full-EM rows, past and future, use batch 0. The long history contains one
row per tissue/grid value/update. The original CSV is preserved byte-for-byte
in the workflow before publication. A commit is applied exactly once after
interruptions. A final full-step M-step estimates the next prior; final 025
fits themselves used that step's input prior.

Rolling cleanup keeps full results for the current and previous scans.
Before 015, retained full results are 014 and 015; before 016, they are 015
and 016. Older full results and no-longer-needed hard-linked seed fits are
retired only after validated successors exist. Priors, histories, counts,
audits and provenance remain. Final retained full scans are 024 and 025.

## Status and reproducibility

```bash
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
squeue --user=wdenault -o "%.18i %.32j %.12T %R"
```

The installed immutable code, configuration, cohort, batch plan, inputs,
submission journals, logs, seed snapshot and prior updates are all under
`PROJECT/workflows/slider_online_015_025`. The original full-EM runtime remains
an explicit dependency while 014 finishes. Raw GTEx data and installed R
packages remain external dependencies; their identities/settings are checked.
See `REPRODUCIBILITY.md` for the existing analysis entry points. No broad
repository reorganization is performed during the active run.

## Validation

Python tests cover missing/old start estimates, held-controller verification,
a controller starting during hold, exact job identity, interrupted handoff,
source-014 worker delegation and reuse, the shared 300-job budget, one-minute
continuations, ambiguous submissions, retention and Linux-readable checksums.
R tests cover all 54 online updates and five full steps, real small slider
warm-start fits, unchanged scientific helper checks, adopted saved outputs,
exactly-once commits and full/batch history numbering. These local tests do
not claim a production RCC execution or a full GTEx fit.
