Expanded slider EM: iterations 006 through 040, with two full checkpoints
=======================================================================

This is a NEW bundle. Do not overwrite or reinstall slider_storage_2026-10-28.
The retirement/compaction/resume jobs already submitted remain unchanged.
No job is cancelled by these scripts. SuSiE-mix is not modified.

What runs
---------
1. A small gate checks every five minutes until the existing slider jobs finish.
   It only reads their state. If they stop without a validated original 005,
   the gate stops with an error instead of starting from incomplete results.
2. Setup verifies the storage completion receipt, deployed sources, R package,
   and all 680 candidate expression IDs/annotation matches/genotype windows.
3. A SEPARATE branch is initialized at iteration_005, using original 005's
   exact frozen priors and its validated gene results. The result files are
   HARD LINKS: no extra full-checkpoint copy, and no original result rewrite.
4. Only the added expression IDs are fitted in the expanded 005 seed. Every
   new entry must validate (including explicit no-eligible-tissue outcomes).
   Only then are original + added counts pooled for the same tissue-wide
   M-step. Original 005's proposed next priors are not used as if its existing
   fits had been computed under them.
5. The branch fits every eligible gene at iteration_006, then 007 ... 040.
   Final target is the actual directory iteration_040, not 40 extra iterations
   and not the previous bundle's target 035. Numbering does not restart at 000.
   Original history 000..005 is retained; expanded history continues to 040.

The candidate map has 680 distinct expression IDs under 679 historical symbols.
Expected in-scope entries: 17,617 existing + 680 added = 18,297. Availability
preflight is mandatory; it does not guarantee informative fits in every tissue.
Chromosome X/Y/MT exclusions remain. Fitting settings, numerical EM helpers,
and convergence validation are unchanged. Missing/failed genes block the
M-step rather than silently shrinking the cohort.

Automatic retention
-------------------
Before submitting the first workers for iteration N, the controller:
  - requires a fully validated checkpoint N-1 and frozen priors for N;
  - confirms no slider worker is still active;
  - snapshots retained validation metadata and hashes all prior/count/history
    files, then records a durable, exact deletion plan;
  - unlinks only old full results/*.rds at N-2 and earlier;
  - repairs N-1's cached dependency stamps and completion certificate;
  - verifies that N-1 remains valid and all prior/count/history files are unchanged;
  - records a completion receipt, then allows N's jobs to be submitted.

Examples:
  while 006 runs: keep full 005 and 006;
  while 007 runs: keep full 006 and 007;
  while 010 runs: keep full 009 and 010;
  after 040: keep full 039 and 040.

All iteration folders, frozen priors, complete per-tissue prior histories,
counts, cohort summaries, audits, and archived completion certificates remain.
Large historical gene/tissue posteriors are permanently removed. This is
retention, not a backup. A second unstripped fit$call copy is unnecessary: the
retained compact results already contain the full scientific fit fields.
The previous storage patch continues to strip only the duplicated saved call.

Source handoff: immediately before expanded 006 starts, original full results
000..005 can retire. Original 005's bytes are still owned by the hard links in
the expanded 005. Original priors/counts/history remain in their original
folders. Before 007 starts, those expanded 005 full files retire as well.
Thus an extra full original checkpoint is not kept forever outside the rolling
policy. Other branches depending on the original results block that handoff.
Unrelated analyses, backups, and simulation results are never deleted.

Interrupted cleanup resumes from its stored plan. It refuses changed files,
redirected/symlink directories, unexpected files, or changed prior history.
No new fits start until metadata repair finishes. Do not remove cleanup state.
Any active job from the original slider run prevents setup/source retirement.

Upload and launch (can be done now or after waiting)
-------------------------------------------------
Upload slider_rolling_rcc_2026-10-28.zip to /home/wdenault, then:

  cd "$HOME" &&
  unzip slider_rolling_rcc_2026-10-28.zip &&
  cd "$HOME/slider_rolling_2026-10-28" &&
  sha256sum -c SHA256SUMS &&
  bash submit_expanded_040.sh \
    /project2/mstephens/wdenault/susie_mix \
    "$HOME/slider_storage_state_2026-10-28"

Submit once. Do not also submit the older expansion command. If Slurm refuses
the initial gate because all 300 user slots are already occupied, inspect
squeue and the pending_submission.json record before retrying; no current job
needs to be cancelled. The scheduler counts all the user's other jobs against
its 300-slot budget, uses at most 298 fit workers, and refills every five minutes.

Paths
-----
Original project:
  /project2/mstephens/wdenault/susie_mix
Expanded branch (new directory):
  /project2/mstephens/wdenault/susie_mix/slider_expanded_040_2026-10-28
Expanded results/history:
  BRANCH/results_slide_prior_em/iteration_005 ... iteration_040
  BRANCH/results_slide_prior_em/prior_history.csv
  BRANCH/results_slide_prior_em/original_prior_history.csv
  BRANCH/results_slide_prior_em/seed_provenance.rds
Cleanup receipts and restart state:
  BRANCH/results_slide_prior_em/rolling_cleanup/before_NNN/
Gate/setup logs and submission receipts:
  /home/wdenault/slider_expanded_040_state_2026-10-28/
Controller logs:
  BRANCH/job/slide_expanded_control_JOBID.out and .err

Checking progress
-----------------
  squeue -u "$USER" -o "%.18i %.26j %.12T %R"
  cat "$HOME/slider_expanded_040_state_2026-10-28/last_submission.json"
  ls -lt "$HOME/slider_expanded_040_state_2026-10-28/logs"

  BRANCH=/project2/mstephens/wdenault/susie_mix/slider_expanded_040_2026-10-28
  cat "$BRANCH/seed_report.csv"
  cat "$BRANCH/results_slide_prior_em/rolling_cleanup/before_007/COMPLETE.json"
  du -h --max-depth=1 "$BRANCH/results_slide_prior_em" | sort -h

A before_007 receipt means full 005 was retired, 006 remains a validated warm
start, and 007 can be fitted. The two retained full checkpoints can still be
large: retention bounds the number, not a fixed number of bytes. Their size
depends on SNPs/tissues and the added genes. Monitor RCC group quota as usual.

Resume after a reported failure
-------------------------------
Inspect the failed job's logs and squeue first. If a pending_submission.json
exists, compare its command with Slurm before clearing it: an ambiguous sbatch
response must not cause duplicate controllers. Never delete results by hand.
Once the issue is resolved and no gate/setup/expanded worker/controller remains:

  bash "$HOME/slider_rolling_2026-10-28/submit_expanded_040.sh" \
    /project2/mstephens/wdenault/susie_mix \
    "$HOME/slider_storage_state_2026-10-28" --resume

An initialized branch is resumed directly; it is not reseeded. Saved valid
current-iteration genes are reused. Interrupted cleanup completes before any
new work is scheduled. Hard fitting failures and exhausted retry budgets still
require diagnosis; the target does not promise convergence by iteration 040.

Validation and limits
---------------------
See TEST_RESULTS.txt. Local synthetic integration tests exercise actual metadata
validation and M-step helpers, actual filesystem hard links, source retirement,
interrupted deletion, repaired completion certificates and warm starts through
005 -> 006 -> 007 -> 008. Slurm/RCC execution and actual GTEx data are validated
by the guarded launch on RCC; they were not accessed from the local machine.
