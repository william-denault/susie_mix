# Shorten the slider schedule through 030

This bundle updates the already installed 2026-11-09 extension inside
`/project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025`.
The historical workflow directory name stays unchanged. All results continue
in `results_slide_prior_em/iteration_XXX`.

| Iteration | Method | Updates per scan | rho |
|---|---|---:|---:|
| 024 | Finish current online scan unchanged | 61 | 0.8 |
| 025–026 | Online | 61 | 0.8 |
| 027 | Online | 61 | 0.5 |
| 028 | Online | 61 | 0.3 |
| 029–030 | Full EM | 1 | 1 |

The workflow stops after completing 030. Each online scan still covers all
18,263 genes once: 24 batches of 300 genes and 37 batches of 299 genes.
There are 400 stages in total, including previously completed updates.

## Run on RCC

Upload `slider_online_taper_025_030_rcc_2026-11-10.zip` to `/home/wdenault`, then:

```bash
cd /home/wdenault &&
unzip slider_online_taper_025_030_rcc_2026-11-10.zip &&
cd slider_online_taper_025_030_2026-11-10 &&
sha256sum -c SHA256SUMS &&
bash start_upgrade.sh
```

This can run while 024 is fitting. It holds a pending controller, lets all
already-submitted workers finish, installs the new frozen schedule, and releases
that controller. It cancels no gene jobs and preserves existing results.
An afterany recovery job completes an interrupted installation or releases the
held controller when no change was made. Once both job IDs are printed, Slurm
owns the continuation and the terminal can be closed.

025 already used rho=0.8, so its complete plan is also preserved. It is safe if
025 starts before installation. The first altered iteration is 026; the gate
refuses to proceed if any 026-or-later iteration or task already exists.
Unrelated ERROR.json or an unresolved submission after workers have drained
is left for diagnosis, rather than cleared by this schedule-only upgrade.

## What is preserved

- The first 215 stage rows, through the end of 025, are unchanged.
- The assignments and within-batch gene order for 026–028 are unchanged;
  only the requested rho values change. No new random partition is generated.
- The frozen runtime estimates and longest-first submission policy are retained.
- The model fitter, scheduler code, R engine and count calculations are
  byte-identical to the installed extension. Only the runtime README changes.
- Current fits, warm starts, learned priors, running sufficient statistics,
  validation receipts, and all existing prior histories are preserved.
- Rolling cleanup retains current and previous full result files and all prior,
  count and audit metadata. The worker cap stays 298 within the 300-job user
  budget; memory remains 48 GB per worker and 16 GB per controller.

029–030 use the complete cohort and ordinary full-data M-steps. As before, the
final M-step prior is recorded separately from the input prior used to fit the
last scan. Completing the requested schedule is not a convergence certificate.

The prior history continues to record iteration, batch, rho, method and input
versus after_update roles. Full EM has batch=0. The next controller refreshes
`results_slide_prior_em/iteration_summary.csv` to the new plan and reports
completed/running/planned status; no statistical output is rewritten during
installation. The previous longer plan is archived for provenance.

## Check installation and progress

```bash
cat /home/wdenault/slider_online_taper_025_030_2026-11-10/launch.json
tail -n 40 /home/wdenault/slider_online_taper_025_030_2026-11-10/logs/gate_*.out
tail -n 40 /home/wdenault/slider_online_taper_025_030_2026-11-10/logs/gate_*.err
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
```

Successful installation prints the new schedule and `400 total update stages`.
Status will show online 025/026 at 0.8, 027 at 0.5, 028 at 0.3, and full 029/030.
While the gate waits or validates, status can still show the previous schedule.
Use the installed code for future status/resume, not an older HOME bundle.

The gate checks the complete previous manifest, including these baseline hashes:

```
af5c1b43551ab2eb8b8490db22cb5e2b64fe93c25ffb746c454eee3d08e7b308  minibatch.py
bf40d692a07fe96befa78523c5efed3cc13b5d2bbdc79387dbc5471276903b8c  engine.R
```

Original code/configuration, previous stages.csv, and the transaction record are
kept in `workflows/slider_online_015_025/repairs/taper_025_030_2026-11-10/`.
Old unused plan files remain for provenance. The bundle contains no GTEx data
or fitted results. Local tests and their scope are listed in VALIDATION.md.
