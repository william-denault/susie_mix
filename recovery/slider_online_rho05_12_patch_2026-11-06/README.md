# Live schedule upgrade: rho 0.5, twelve batches

Upgrade the existing `workflows/slider_online_015_025` run. All already-started
iterations keep their original plan and priors. From the next unstarted online
iteration through 020, use rho 0.5 and twelve shuffled, equal-size batches:
eleven batches of 1,522 genes and one of 1,521. Each of the 18,263 genes appears
once per pass. The balancing is by gene count, not credible-set count.

With 016 already running and 017 absent at installation, the schedule is:

| Iterations | Batches | rho |
|---|---:|---:|
| 015–016 | 9 | 0.3 |
| 017–020 | 12 | 0.5 |
| 021–025 | one full-data step | 1 |

This is 71 update stages, including the stages already committed. If 017 has
started before the gate secures the workflow, it too is preserved and the
change begins at 018. The audit records the actual boundary. If all online
iterations have started, the upgrade refuses to change the plan and releases
the existing controller.

## Upload and launch

Upload the ZIP to `/home/wdenault`, then run this block once:

```bash
cd /home/wdenault &&
unzip slider_online_rho05_12_patch_rcc_2026-11-06.zip &&
cd /home/wdenault/slider_online_rho05_12_patch_2026-11-06 &&
sha256sum -c SHA256SUMS &&
bash start_upgrade.sh
```

The launcher submits a gate and an afterany recovery job to Slurm. After it
prints both job IDs, the terminal can be closed. Two free submission slots
are required; if unavailable, it makes no change and asks you to retry later.

The expected installed baseline is the cleanup patch, with these hashes:

```
665d688f21cd2397bb5dd54e966320c55a86b054990ba4b57f37ad9685c7b8bd  code/minibatch.py
24f1a74aff8fcb4c669d1ccc7c25c7cb2aeacca0ca00992b14d1cc1e9b72e1be  code/engine.R
```

The gate verifies all baseline code and frozen inputs. It holds only the exact
pending controller, lets submitted gene jobs and the running controller finish,
then installs under the workflow locks. No jobs are cancelled. Completed fits,
running sufficient statistics, current batch inputs and prior histories are
preserved. No R fitting or M-step is performed by the installer.

Changes to code, configuration, future plan files and checksums are journaled.
The recovery job completes an interrupted installation before releasing the
same controller. An existing ERROR.json or ambiguous submission is never erased.
Do not manually edit checksums, configuration or stage files.

## Check progress

```bash
squeue --user=wdenault -o "%.18i %.32j %.12T %R"
tail -n 30 /home/wdenault/slider_online_rho05_12_patch_2026-11-06/logs/gate_*.out
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
cat /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/repairs/rho05_12_batches_2026-11-06/INSTALLED.json
```

Before the gate runs, its log may not exist. During the drain, the controller
will show JobHeldUser and gene jobs continue. After installation, status lists
the batch count/rho for every iteration and the correct total stage count.
The same results folder, one-minute refill eligibility, 298-worker cap,
300-user-job budget, scheduler retries, and current/previous result retention
remain in effect. No guarantee is made about convergence speed or Slurm uptime.

## Histories and plots

Histories retain every existing prior value and add `batches_in_iteration`
(9 or 12 for online passes; 0 for full EM). `rho`, `iteration`, `batch` and
`update_method` reflect each actual update. Full EM keeps public batch zero.
The new column is added when the released controller next publishes history;
the installer does not rewrite histories.

Optional plotting scripts in `analysis/` recognize both 9- and 12-batch ends.
They are not automatically installed over project scripts. The statistical
fitter, likelihood, eligible-count rules and warm-start validation are unchanged.
The update is `S_new = (1-rho)*S_old + rho*batch_counts/batch_gene_count`, then
normalization within tissue. Five ordinary pooled-count full EM updates follow.

## Reproducibility and local validation

`repairs/rho05_12_batches_2026-11-06/` retains the old code, original config,
original stages/checksums, checkpoint metadata, new schedule and transaction.
Old plan files remain; new future plans use `plan/rho05_12/`. The modified
runtime remains self-contained under the existing project's workflow `code/`.

Tests cover the complete 18,263-gene partition and preserved plan prefix,
held-controller draining, recovery after code/config/stage-write interruption,
dynamic completion through 025, and the existing scheduler timeout handling.
R transition fixtures exercise 71 updates, per-gene warm-start priors,
exact sufficient-statistic updates, zero-information fallback, old-result
cleanup, temporary files, and interrupted commits. A separate small synthetic
model test performs 198 actual slider fits across 015–025. These are local
tests, not a full GTEx execution or a live RCC scheduler test.
