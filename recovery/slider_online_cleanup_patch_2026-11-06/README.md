# Safe cleanup upgrade for the running slider online workflow

Install only on the queue-patched workflow confirmed on RCC on 2026-11-06:
minibatch.py SHA256 ae0cb0b4b506c7701e77c70e05a9657d0ef1b5412076e4b6e914933f38bf3a6b.
The installer verifies the entire baseline manifest and frozen configuration.

The run was healthy: seven updates committed, batch 8 in progress. Files named
.em_write_* in active iteration 015 may still be in use. Do not delete them.
If a killed worker leaves such a file after a successful retry, the old cleanup
would reject it when retiring that iteration. The new cleanup handles only the
reserved alphanumeric temporary-file format, and only in obsolete results.

## Launch once, while the current gene jobs keep running

Upload slider_online_cleanup_patch_rcc_2026-11-06.zip to /home/wdenault:

```bash
cd /home/wdenault &&
unzip slider_online_cleanup_patch_rcc_2026-11-06.zip &&
cd /home/wdenault/slider_online_cleanup_patch_2026-11-06 &&
sha256sum -c SHA256SUMS &&
bash start_upgrade.sh
```

Two free submission slots are needed. If the launcher reports no room, the EM
is unchanged; retry when jobs finish. Once it prints both job IDs, Slurm manages
the upgrade and you can close the terminal. Do not submit duplicate upgrade jobs.

## How the upgrade preserves active work

1. Submit the upgrade job held, then its afterany recovery job; release the
   upgrade only after both submission receipts are saved.
2. Verify and temporarily hold exactly one pending online controller. Existing
   running/queued gene jobs and a currently running controller finish normally.
3. Wait until no branch fitting jobs or running controllers remain. New gene
   submissions pause during this drain. No gene job is cancelled.
4. Under the workflow locks, preserve analysis metadata, archive the old code
   and configuration, install the new code and corresponding fingerprints, and
   verify that all saved priors, counts and state are unchanged.
5. Release the same pending controller. It reuses saved result receipts, commits
   any completed batch once, and continues the original schedule.

The recovery job restores continuation after an interrupted upgrade: it releases
the held controller if installation never began, or finishes the recorded code
transaction before releasing it. It is normally a short no-op after success.
Unknown code, an existing ERROR.json, ambiguous submissions, unexpected jobs or
changed analysis metadata cause a refusal rather than an unsafe overwrite.

The gate allows 18 hours for worker draining within a 23-hour allocation. If it
cannot drain in that time, it restores the old continuation without installing.
Inspect its logs rather than resubmitting blindly. Slurm/node/filesystem failures
cannot be ruled out, and an integrity error still requires inspection.

## Monitor

```bash
tail -n 40 /home/wdenault/slider_online_cleanup_patch_2026-11-06/logs/gate_*.out
tail -n 40 /home/wdenault/slider_online_cleanup_patch_2026-11-06/logs/gate_*.err
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
```

During draining, the pending controller intentionally shows JobHeldUser.
Success prints "Cleanup hardening installed" and "Released existing controller".
Backup/transaction records are under
PROJECT/workflows/slider_online_015_025/repairs/cleanup_transition_2026-11-06/.

Rho remains 0.3; two batches of 2,030 and seven of 2,029 genes per pass;
online 015-020, then full EM 021-025. Current and previous full results remain;
every prior history/count file remains. No R fitting/update code changes.

The local audit combined the actual Python cleanup with all 59 R update stages,
interrupted cleanup and commit recovery, and 99 real slider fits on small
synthetic data. Tests also injected abandoned R temporary files at every scan
transition. They do not constitute a rerun of GTEx or a guarantee of RCC uptime.
