# Scheduler query recovery for slider online EM

This bundle repairs the observed squeue timeout from controller 49863377.
It requires the receipt-patched code, five committed updates, and six fully
validated batches (12,176 saved genes) in iteration 015. It refuses a different
error, active slider jobs, changed outputs, or an unresolved submission journal.

Upload slider_online_queue_patch_rcc_2026-11-06.zip to /home/wdenault, then:

```bash
cd /home/wdenault &&
unzip slider_online_queue_patch_rcc_2026-11-06.zip &&
cd /home/wdenault/slider_online_queue_patch_2026-11-06 &&
sha256sum -c SHA256SUMS &&
bash start_repair.sh
```

The last line prints a Slurm repair job ID. The submitted job survives closing
the terminal. Do not submit another launcher while this repair is running.

## What changes

Only runtime/minibatch.py changes executable behavior. Each squeue call has
a 60-second response deadline instead of 15 seconds. Known temporary queue
failures retry after 10, 20, 40, then 60 seconds between attempts. A controller
keeps waiting within its existing 23-hour allocation; it does not write a
permanent ERROR.json for these failures. Its already-queued afterany successor
can continue if the controller reaches its wall time. An outage spanning the
entire allocation before any successor can be submitted can still require a
manual resume; cluster outages cannot be eliminated by this script.

Retries repeat only the read-only squeue operation, not earlier submissions,
prior updates, cleanup or result writes. Failed/partial output never means an
empty queue. Invalid command responses and analysis-integrity errors still stop.
Interactive status and launch commands retry for approximately five minutes,
then report the unavailable scheduler without changing analysis state.

sbatch behavior is unchanged: its existing 120-second timeout and pending
submission journal remain. An ambiguous submission is never blindly resubmitted.
The earlier warm-start retention and literal T/TBXT receipt fixes are included.

## Preservation and resumption

The repair verifies code and frozen configuration, hashes saved metadata, and
checks every result against its saved validation receipt. It dry-runs the real
batch-6 commit in memory, validating all 2,029 summaries. It does not reread the
large fitted objects or refit genes. The dry run prints "Committed stage 6" from
the calculation function, but redirects all seven writes into memory; on-disk
state remains at stage 5 until the normal controller actually commits it.

After validation, the repair archives the original code/config/error, installs
the new runtime with matching checksums, and submits a controller. It commits
batch 6 once, then proceeds to batch 7. Rho remains 0.3; online passes remain
015-020 with nine batches each, followed by full EM 021-025. Prior history,
gene membership, scientific fitting/count code, and rolling cleanup are unchanged.

Backup and dry-run previews:
PROJECT/workflows/slider_online_015_025/repairs/scheduler_queue_retry_2026-11-06/

## Monitoring

```bash
tail -n 40 /home/wdenault/slider_online_queue_patch_2026-11-06/logs/repair_*.out
tail -n 40 /home/wdenault/slider_online_queue_patch_2026-11-06/logs/repair_*.err
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
```

Expected repair output includes "12176 validated results preserved; five prior
updates unchanged" and a submitted controller ID. Status should subsequently
advance to 6/59 committed stages and batch 7. RETRY messages in controller logs
now indicate a temporary scheduler delay, followed by a recovery message when
the queue becomes available.

## Tests

Python fault-injection tests cover timeouts, repeated communication failures,
permanent errors, fresh queue reads, no repeated commit or duplicate submission,
existing-fit reuse, interrupted installation, and preservation guards.
The R fixture exercises the production batch-6 dry run and all 54 online updates
plus five full steps, with independent count/prior checks and literal T/F genes.

The deployed runtime is PROJECT/workflows/slider_online_015_025/code/; include
that folder when copying scripts for reproducibility.
