# Installed slider EM workflow: extended through 038

The frozen schedule in ../stages.csv is authoritative. The workflow retains
its historical directory name slider_online_015_025 to preserve checkpoint paths.

| Iteration | Method | Updates per scan | rho |
|---|---|---:|---:|
| 015-016 | Online, historical | 9 | 0.3 |
| 017-020 | Online, historical | 18 | 0.5 |
| 021-023 | Full | 1 | 1 |
| 024-025 | Online | 61 | 0.8 |
| 026-028 | Online | 61 | 0.5 |
| 029-033 | Online | 61 | 0.3 |
| 034-038 | Full | 1 | 1 |

There are 708 stages including the preserved 93 through 023. Each new online
scan uses 24 batches of 300 and 37 batches of 299 genes. The cohort is unchanged.
All final five full EM steps are retained. Results remain in the project's
results_slide_prior_em/iteration_XXX directories.

Run bash status.sh to inspect progress. Resume only an inspected, error-free
workflow using this installed start.sh --resume; do not launch older HOME copies.
An unacknowledged worker submission is reconciled from positive Slurm accounting
and its unique task logs. It is never blindly resubmitted. The reported accepted
array 50014103 is registered by the extension upgrade, with its evidence archived.

The model and sufficient-statistic update are unchanged:
S_new = (1-rho)*S_old + rho*batch_counts/batch_gene_count.
Weights are normalized separately for each tissue. Full steps use pooled counts.
No additional batch/total-gene multiplier or credible-set filter is introduced.
Warm starts retain the exact previous result and the prior that produced it.

../plan/extend_024_038_2026-11-09/runtime_estimates.csv contains frozen historical
worker durations. New batch membership is reproducibly randomized first; pending
jobs are then submitted longest predicted runtime first within each batch.
Unknown durations use the observed median. Slurm controls actual start order.

results_slide_prior_em/iteration_summary.csv reports method, rho, batch count,
updates completed and completion status. Original prior histories keep their
weights and within-scan updates. Full rows have batch=0; online batch=0 is the
input prior, while positive batch numbers with prior_role=after_update are learned
updates. At 038 the saved fits use its input prior; the final M-step is recorded
separately. Rolling cleanup keeps current/previous result files and all prior
history, counts, audits and metadata. One-minute refill eligibility, at most 298
workers and the all-user 300-job budget remain.

The upgrade record is ../repairs/extend_024_038_2026-11-09. Previous code and
configuration are retained there. README_original_launch.md is historical only.
REPRODUCIBILITY.md maps the other analyses. This bundle does not move them.

The runtime is self-contained. To export a copy of its code without old HOME
bundles, run python3 build_bundle.py /a/new/destination_folder; it also writes a
ZIP beside that folder. Resuming this analysis additionally needs the existing
workflow state, frozen plans, data and package versions, not just the code ZIP.
