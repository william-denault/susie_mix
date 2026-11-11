# Installed slider EM workflow: final iteration 030

The frozen schedule in ../stages.csv is authoritative. The historical directory
name slider_online_015_025 remains so saved checkpoint paths stay valid.

| Iteration | Method | Updates per scan | rho |
|---|---|---:|---:|
| 015-016 | Online, historical | 9 | 0.3 |
| 017-020 | Online, historical | 18 | 0.5 |
| 021-023 | Full, historical | 1 | 1 |
| 024-026 | Online | 61 | 0.8 |
| 027 | Online | 61 | 0.5 |
| 028 | Online | 61 | 0.3 |
| 029-030 | Full | 1 | 1 |

The plan has 400 stages including historical updates and stops after 030.
Each online scan covers all 18,263 genes once: 24 batches of 300 and 37 of 299.
The 2026-11-10 change preserves all stages through 025, keeps the existing
026-028 batches and runtime ordering, and changes only future rho values and
the final two full steps. No learned count, fit, or prior is changed by installing
this schedule. Previous longer plans are retained in the repairs directory.

Run bash status.sh for progress. Resume only an inspected, error-free workflow
using this installed start.sh --resume; do not launch older HOME copies.
Worker acknowledgement recovery, scheduler retries, the 298-worker limit,
one-minute refill eligibility and 300 all-user job budget are unchanged.
The fitter, R engine, scheduler and statistical helpers are identical to those
in the 2026-11-09 extension; memory allocations remain 48 GB/worker, 16 GB/controller.

The online update is S_new = (1-rho)*S_old + rho*batch_counts/batch_gene_count,
normalized separately per tissue. Full steps use pooled counts. Every gene
warm-starts from its previous fit and validates against the prior that produced
that fit. No credible-set filter or extra batch/total-gene multiplier is added.

../plan/extend_024_038_2026-11-09/runtime_estimates.csv retains the historical
runtime estimates. Pending jobs are submitted longest predicted runtime first
within each fixed batch; Slurm controls actual start times. Unknown runtimes
use the observed median. This schedule update does not resample batches.

results_slide_prior_em/iteration_summary.csv shows method, rho, batch count and
completion status. Detailed histories preserve all tissue weights and slider
grid values. Full-EM rows use batch=0. For online EM, batch=0 is the input prior;
positive batch numbers with prior_role=after_update are learned updates.
Completing 030 does not by itself establish convergence. Its final M-step prior
is recorded separately from the prior used for the 030 fitted results.

Rolling cleanup keeps current/previous gene result files and all prior histories,
counts and audit metadata. Runtime code and analysis paths remain inside the
project. No simulation or SuSiE-mix analysis is modified.

The current schedule record is ../repairs/taper_025_030_2026-11-10. Earlier
repairs retain the old source and plans. REPRODUCIBILITY.md maps the other
analyses; references there to the earlier extension are historical provenance.
To export this self-contained code, use python3 build_bundle.py NEW_FOLDER.
Resuming also requires the existing workflow state, data and package environment.
