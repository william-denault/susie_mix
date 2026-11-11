# Extend the slider workflow through iteration 038

This bundle is an upgrade for the existing workflow:
`/project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025`.
The historical directory name stays unchanged so saved paths remain valid.
Results continue in `results_slide_prior_em/iteration_XXX`.

| Iteration | Method | Batches / updates | rho |
|---|---|---:|---:|
| 023 | Finish the already-started full EM | 1 | 1 |
| 024–025 | Online EM | 61 per scan | 0.8 |
| 026–028 | Online EM | 61 per scan | 0.5 |
| 029–033 | Online EM | 61 per scan | 0.3 |
| 034–038 | Full EM | 1 per scan | 1 |

Every online scan covers all 18,263 genes once: 24 batches of 300 genes and
37 batches of 299. Balanced here means equal gene counts, using seeded random
assignment; it does not imply equal statistical information across tissues.
The original 93 stages through 023 stay unchanged. The complete extended
plan has 708 stages, including historical completed stages.

## Install on RCC

Upload `slider_online_extend_024_038_rcc_2026-11-09.zip` to `/home/wdenault`.
Run once:

```bash
cd /home/wdenault &&
unzip slider_online_extend_024_038_rcc_2026-11-09.zip &&
cd slider_online_extend_024_038_2026-11-09 &&
sha256sum -c SHA256SUMS &&
bash start_upgrade.sh
```

It submits a gate and an `afterany` recovery job. After both job IDs are
printed, Slurm owns the work and the terminal can be closed. Do not launch
an older bundle or delete ERROR.json/pending_submission.json by hand.
The gate refuses to replace 024 if any 024-or-later iteration or task exists.

The exact installed baseline is checked, including:

```
00b3727c73c5af92120ddeb90d0335bfba91e9f621af0f820facd8a600892cf5  minibatch.py
3f1998220ca4f83948eab3992999b307f7ab6a517f9b8003c23f4ff1193fc828  engine.R
```

## Recovery of the reported 023 stop

Controller 50013914 received `Unexpected message received` from sbatch, but
Slurm accepted task_001411 as array 50014103. RCC accounting reported all
130 workers completed with exit code 0:0. Before modifying the stopped state,
the gate rechecks the exact pending command, array identity, all 130 accounting
records, the frozen input checksum, and every worker's saved validation receipt.
It registers that accepted array and archives the matching error and pending
submission records. Unrelated errors are not cleared. Existing 023 outputs
are reused; uncompleted genes are scheduled normally.

The code/configuration/plan change and receipt repair are journaled together.
The recovery job can finish an interrupted installation. On a running workflow,
the gate holds only its pending controller and drains existing gene jobs before
swapping code. No worker job is cancelled. Current and previous fitted results,
running sufficient statistics and prior history are preserved by the upgrade.

Future worker submissions with a lost acknowledgement are reconciled using
their unique task log directory and Slurm accounting. An ambiguous submission
is never blindly repeated. If positive evidence remains unavailable for 20 hours,
the workflow keeps its pending journal and stops for diagnosis. No code can
guarantee that cluster outages or new input/model failures will never stop a run.

## Runtime ordering

The gate estimates each gene's duration from its latest successful worker
attempt and validation receipt timestamps. These measurements include R startup
and input/output, not just fitting. The newest usable observation is used;
genes with no usable timing use the observed median. The resulting CSV is frozen
at `plan/extend_024_038_2026-11-09/runtime_estimates.csv` inside the workflow.

Future batch membership stays fixed and randomized; within each batch the
longest predicted genes are submitted first. The same ordering also applies to
remaining unscheduled work in 023 and to the final full steps. Already-submitted
tasks are untouched. Slurm determines actual start times and available resources.
With batches close to the 298-worker cap, a slow final gene can still delay an
update; this ordering reduces that risk but does not guarantee shorter scans.

## History and storage

`results_slide_prior_em/iteration_summary.csv` lists each iteration's method,
rho, batch count, gene counts, updates completed, and planned/running/completed
status. Historical pre-workflow iterations are labelled `historical` and 014
`source_checkpoint`, rather than claiming to revalidate their retired results.
The file is refreshed by the controller, and copies through the relevant scan
are saved alongside each iteration's prior history. Counts describe the fixed
gene cohort, not a guarantee that every tissue has an eligible fit.

Existing `prior_history.csv` and `prior_history_long.csv` retain per-tissue
weights, iteration, batch, rho, phase, prior_role, update_method, and slider grid.
Full EM rows use batch=0; online rows include input batch=0 and every completed
batch update. Compare `prior_role="after_update"` rows when measuring changes.
At a scan's start the input prior equals the previous update by design.

Rolling cleanup continues to keep current and previous full gene result files,
preserving priors, counts, history, audit and checkpoint metadata for all scans.
Each gene's warm start is validated against its own previous batch prior.
The scientific fitter and count/update rules are unchanged. A full scan's final
M-step prior is recorded separately from the input prior used to fit that scan.

## Check progress

```bash
cat /home/wdenault/slider_online_extend_024_038_2026-11-09/launch.json
tail -n 40 /home/wdenault/slider_online_extend_024_038_2026-11-09/logs/gate_*.out
tail -n 40 /home/wdenault/slider_online_extend_024_038_2026-11-09/logs/gate_*.err
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
```

After installation, status shows the extended schedule and eventual final
iteration 038. In RStudio:

```r
root <- "/project2/mstephens/wdenault/susie_mix/results_slide_prior_em"
read.csv(file.path(root, "iteration_summary.csv"))
length(list.files(file.path(root, "iteration_023", "results"), pattern = "\\.rds$"))
```

Use `results_slide_prior_em` for slider; `results_em` is the separate SuSiE-mix EM.
The audit, original code/configuration and archived submission error are under
`workflows/slider_online_015_025/repairs/extend_024_038_2026-11-09`.

## Validation

See `VALIDATION.md` for the executed local checks. All tests use synthetic
fixtures or mocked Slurm; no local test submits jobs to RCC. The bundle includes
the full runtime, model helpers, gene maps, job stubs and tests. No GTEx data or
fitted results are packaged, and no project directories are reorganized.
