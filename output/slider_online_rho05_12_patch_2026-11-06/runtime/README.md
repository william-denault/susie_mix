# Installed slider online EM workflow

This runtime supports a mixed schedule recorded in config.json and stages.csv.
The rho05_12_batches_2026-11-06 upgrade preserves every started iteration and
changes only the remaining online passes to rho=0.5 with 12 equal-size shuffled
batches (11 x 1,522 genes; 1 x 1,521 genes). Full steps 021-025 follow unchanged.
For an upgrade during 016, the plan has 71 stages: 18 old online updates,
48 new online updates, and five full updates. Status uses the actual plan length.

Run status.sh in this directory to inspect the per-iteration settings and progress.
Existing live jobs are managed by Slurm. Resume a stopped, error-free run only
using this installed code/start.sh --resume; do not relaunch older HOME bundles.
An ERROR.json or ambiguous submission still requires diagnosis before resuming.

The model, count rules and warm starts are unchanged. Each batch freezes its
prior before submission; every gene is validated before the count update.
S_new = (1-rho)*S_old + rho*batch_counts/batch_gene_count, normalized per tissue.
There is no additional batch_size/total_genes multiplier. Full steps use pooled
counts. The 025 result fits use the input prior; its final M-step is separately
recorded in the history.

History records iteration, batch, rho, phase, prior_role, update_method and
batches_in_iteration. Full-data rows use batch=0. Completed online passes are
identified by batch==batches_in_iteration and prior_role=="after_update";
a partial pass must not be marked complete just because it is the latest row.

Rolling cleanup keeps current/previous full results and every prior, count,
history and audit. Warm starts validate each gene against its own previous
batch's prior. One-minute refill eligibility, 298 worker maximum and 300 all-user
job budget remain. Transient scheduler query failures are retried; uncertain
submissions are reconciled rather than blindly repeated.

The installed migration record and original code/configuration are under
../repairs/rho05_12_batches_2026-11-06/. Old plan files are retained; replacements
for unstarted passes live under ../plan/rho05_12/. No fitting occurs during the
upgrade, and no submitted gene job is cancelled.

README_original_launch.md describes the original rho=.3 / 9-batch launch and
its earlier repairs, for provenance only. It does not override the current
frozen schedule. REPRODUCIBILITY.md maps the other project workflows.
