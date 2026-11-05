# Repair the online 015 seed-retention failure

The 015 recovery controller mistakenly released the hard-linked 014 warm
starts at the beginning of 015. The threshold still referred to the previous
014-start schedule. Workers then failed before fitting, and the three-attempt
limit stopped the run. Source 014 was certified (18,263 genes / 490,509 tissue
fits), and its original results remain in iteration_014/results.

This repair is specifically for the reported stage-zero failure. It refuses
active slider jobs, an ambiguous submission, changed code/state, any completed
015 results/receipts, or an already-committed online update.

## Run on RCC

Upload slider_online_seed_patch_rcc_2026-11-05.zip to /home/wdenault:

```bash
cd /home/wdenault &&
unzip slider_online_seed_patch_rcc_2026-11-05.zip &&
cd /home/wdenault/slider_online_seed_patch_2026-11-05 &&
sha256sum -c SHA256SUMS &&
bash start_repair.sh
```

The last line prints the Slurm repair-job ID. The terminal can be closed once
the job has been submitted. After it starts:

```bash
tail -n 40 /home/wdenault/slider_online_seed_patch_2026-11-05/logs/repair_*.out
tail -n 40 /home/wdenault/slider_online_seed_patch_2026-11-05/logs/repair_*.err
```

Expected messages include `Repair verified: 18263 seed links restored` and
`Submitted controller JOBID`. A queued repair job has no output file yet.
The repair resumes the existing workflow automatically. Do not relaunch the
old HOME bundle: its manifest intentionally no longer matches the patched
installed workflow.

```bash
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
```

## What is preserved and changed

The repair checks all source file sizes/timestamps against the frozen capture
manifest, then restores missing hard links from iteration_014/results. It
does not copy/recompress the RDS data or refit 014. Each worker still checks
the saved MD5 before using its restored warm start. Hard links require only
directory metadata; both names refer to the same existing data.

State, sufficient statistics, priors, the batch plan, input.rds and history
are hashed before/after and must remain unchanged. A base-R preflight checks
stage 0, the source cohort, warm-start mappings, initial counts/weights, and
015 batch 1's frozen input. It does not fit any model.

The source seed is now retained throughout all nine batches of 015. It can
be released at 016 only after every stage of that first pass is committed.
The two-full-iteration storage policy otherwise stays unchanged. Missing
warm starts also produce an explicit error during preparation, before a
new worker array is submitted, and in the worker itself.

Failed task submissions and logs are moved, not deleted, into:

```
PROJECT/workflows/slider_online_015_025/repairs/seed_retention_2026-11-05/
    before/                      # original config, frozen hashes and ERROR
    code.before/                 # original installed code
    failed_tasks/                # all failed 015 task records and logs
    ERROR.before_repair.json
    PLAN.json
    COMPLETE.json
```

Moving failed tasks out of the active task directory resets only their
exhausted retry counters. The specific recorded error is archived after all
checks pass. No prior update is undone or replayed. The code manifest and the
configuration's code fingerprint are updated together, with backups; the
new configuration hash is recorded in frozen_inputs.json. Other frozen
inputs are not rewritten. Interrupted installation is resumable by submitting
this repair again; an already-complete repair does not reset tasks again.

The schedule remains online 015-020, rho=0.3, nine balanced batches per pass,
then full 021-025, one-minute refill eligibility, 298-worker maximum and
300 shared queued/running slots. Scientific fitting and count functions
are unchanged. Existing input checksums and batch assignments stay fixed.

## Local validation

Twenty runtime Python tests and seven repair tests passed, including all
nine first-pass batches retaining seed links, refusal of early release at
016, later release, full repair with histories unchanged, changed-source and
active-job rejection, and resuming an interruption during code replacement.
R checks validated actual saved stage-zero objects and rejected a committed
stage. A real small slider fit succeeded after invoking the Python cleanup
in the production order; both online/full updates and history checks passed.
Missing-input error messages were also verified. These are local tests,
not a claim that the repair has already executed on RCC.
