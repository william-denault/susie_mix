# Slider schedule upgrade: 18 batches, rho 0.5

This bundle upgrades the already-installed rho05_12 patch for
`/project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025`.
It preserves every started iteration, its batches, completed fits, running
sufficient statistics, warm starts and history. Only unstarted online passes
are repartitioned. It does not cancel any job or refit a completed gene.

With 016 running and 017 absent when the gate secures the workflow:

| Iterations | Updates per pass | rho |
|---|---:|---:|
| 015–016 | 9 | 0.3 |
| 017–020 | 18 | 0.5 |
| 021–025 | one full-data update | 1 |

Each new online pass contains all 18,263 genes exactly once: eleven batches of
1,015 and seven of 1,014 genes. Batches are shuffled with the recorded seed and
balanced by gene count. This is approximately 1,000 genes per batch, not a strict
1,000-gene maximum. The expected total plan has 95 updates, including completed
updates. If 017 starts before installation, its existing 12-batch plan is kept
and 18 batches begin at 018 instead. The actual boundary is recorded.

## Launch on RCC

Upload the ZIP to `/home/wdenault` and run once:

```bash
cd /home/wdenault &&
unzip slider_online_rho05_18_patch_rcc_2026-11-06.zip &&
cd /home/wdenault/slider_online_rho05_18_patch_2026-11-06 &&
sha256sum -c SHA256SUMS &&
bash start_upgrade.sh
```

Two free submission slots are required. After the launcher prints both gate
and recovery job IDs, Slurm manages the upgrade and the terminal can be closed.
The gate waits for the previous upgrade jobs, verifies installed code/config,
holds only the exact pending online controller, and waits for submitted gene
jobs and the current controller to finish. It then installs under workflow
locks and releases the same held controller. Existing jobs are never cancelled.
The separate afterany recovery job can finish an interrupted installation.

## Check installation and progress

```bash
squeue --user=wdenault -o "%.18i %.32j %.12T %R"
tail -n 30 /home/wdenault/slider_online_rho05_18_patch_2026-11-06/logs/gate_*.out
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
cat /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/repairs/rho05_18_batches_2026-11-06/INSTALLED.json
```

The gate log should report `Schedule upgrade installed: rho=0.5, 18 batches`
and `Released existing controller`. While it drains workers, the pending
controller may show JobHeldUser. Before installation, status shows the old
schedule; afterwards it lists the actual settings for each iteration.

The result folder, one-minute refill eligibility, 298-worker maximum, all-user
300-job budget and existing scheduler retry handling are unchanged. Cleanup
still keeps current and previous full results, preserving all histories,
counts and audits. Full EM ends at 025. Scheduling and convergence speed are
not guaranteed by reducing the batch size.

## Scientific behavior and provenance

The executable runtime is identical to the installed twelve-batch version;
only its README is updated. Model fitting, eligibility/count rules, warm starts,
validation, commits and cleanup are unchanged. New batches use
`S_new = 0.5*S_old + 0.5*batch_counts/batch_gene_count`, followed by within-tissue
normalization. There is no extra batch-size/total-gene-count multiplier.

History automatically records batch 1–18, rho 0.5 and
`batches_in_iteration=18`. Existing nine- or twelve-batch history remains labeled
correctly. Full-data EM uses public batch 0. The plot scripts from the preceding
bundle already use this column to identify completed online passes.

The previous plan/code/configuration and the new transaction are retained under
`repairs/rho05_18_batches_2026-11-06`. All old gene-list files remain; new plans
use `plan/rho05_18/`. The earlier twelve-batch patch's audit remains intact.
An existing error or ambiguous submission is not cleared by this upgrade.

## Local validation

Tests cover the 18,263-gene partition, preservation of started nine/twelve-batch
passes, all five full steps, previous-upgrade coordination, uninterrupted and
interrupted code/config/plan installation, and continuation beyond the old
71-stage stopping point. Synthetic R fixtures exercise the complete 95-stage
plan, per-gene warm starts, exact sufficient-statistic updates, literal T/F
identifiers, zero-information fallback, history and cleanup through 025.
A separate small real-model test performs 198 slider fits across 015–025.
These are local tests, not a live RCC or full GTEx execution.
