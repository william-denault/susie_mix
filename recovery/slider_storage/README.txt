Slider storage recovery and renamed-gene expansion (2026-10-28)
=============================================================

This bundle has not been run on RCC. It is built from the deployed Oct-23
shared-scheduler bundle, not the optional tissue-splitting bundle.

What is preserved
-----------------
The storage patch removes only fit_slide_prior$call. It preserves ALL other
fields: joint posteriors, all grid moments, alpha, PIPs, credible sets, variances,
convergence diagnostics, tissue/sample/SNP metadata and warm-start arrays.
Each rewritten result is deserialized and compared exactly with the original
minus the call. Existing gene-result hashes, summary timestamps and completion
certificates are repaired. No model fitting occurs during maintenance.

Frozen priors and prior histories are checked byte-for-byte before/after:
  results_slide_prior_em/prior_history.csv
  results_slide_prior_em/iteration_005/prior_history.csv
  results_slide_prior_em/iteration_005/priors.csv
The history is in the iteration directory, not inside its results/ subfolder.
Histories/counts/settings/manifests for earlier iterations are also preserved.

The installed workhorse patch clears call metadata from warm-start inputs and
new outputs. Package version, grid, convergence tolerance, QC and M-step do not
change. Known deployed source hashes are checked before any maintenance.
An unexpected source version stops without overwriting it.

The optional retirement removes ONLY *.rds inside results/ for iterations
000--003. It permanently removes those historical full posteriors, so use it
only if they are no longer needed or are independently backed up. It keeps
all iteration directories, priors, histories, pooled counts, audit information
and summaries, then repairs the retained iteration 004's dependency stamps.
Iteration 004 is retained for unfinished iteration 005 warm starts.

1. Upload and unpack in HOME (project2 is quota-full)
---------------------------------------------------
Upload slider_storage_rcc_2026-10-28.zip to your RCC home directory.

  cd "$HOME"
  test ! -e slider_storage_2026-10-28 && \
    unzip slider_storage_rcc_2026-10-28.zip
  BUNDLE="$HOME/slider_storage_2026-10-28"
  PROJECT=/project2/mstephens/wdenault/susie_mix
  STATE="$HOME/slider_storage_state_2026-10-28"
  SCRATCH="/scratch/midway2/$USER/slider_storage_2026-10-28"
  cd "$BUNDLE" && sha256sum -c SHA256SUMS

2. Choose ONE recovery submission
--------------------------------
Recommended if you want to discard full results from 000--003 first (saves
their compaction I/O and immediately frees project quota):

  bash "$BUNDLE/submit_recovery.sh" "$PROJECT" --retire-through-003

Alternatively, retain every iteration's full posterior and compact them all:

  bash "$BUNDLE/submit_recovery.sh" "$PROJECT"

The script queues optional retirement -> compaction -> resume, using afterok
dependencies. Failure prevents later steps from running; newly queued dependent
maintenance jobs are cancelled if their dependency becomes impossible (Slurm
--kill-on-invalid-dep=yes, https://slurm.schedmd.com/sbatch.html). Logs and state go
under HOME; original bytes during a replacement are temporarily escrowed in
SCRATCH. Two single-thread R processes compact different genes in one Slurm
job (2 CPUs, 48 GB, 23 hours). This fits broadwl's 64 GB nodes; see
https://docs.rcc.uchicago.edu/slurm/partitions/ . If memory is exceeded, resume
with --workers 1; completed transactions are reused. Override resources when
resubmitting if needed.
Existing jobs are never cancelled. Maintenance refuses to run while any
slider controller/worker is still active and holds the preparation lock.

The scratch design works even when project quota is full: the original file
is copied to scratch and checksum-verified, the compact file is round-trip
verified, and only then is the original project file removed to release its
quota. A compact temporary file is copied back, verified and atomically
renamed. Escrow is removed only after a durable committed receipt. A failed
write or killed job can temporarily leave the project path absent; the exact
original and compact copies remain in scratch for resumption. Do not remove
STATE or SCRATCH while an operation is incomplete. This is not a substitute
for a long-term backup. Check scratch/home quota before starting if those
allocations have changed (rcchelp quota).

3. Check completion, then iteration 005 resumes automatically
------------------------------------------------------------
  cat "$STATE/queued_jobs.tsv"
  cat "$STATE/COMPLETE.json"
  cat "$STATE/validation_report.csv"
  cat "$STATE/logs/"resume_*.out

COMPLETE.json is written only after compaction, metadata repair and history
checks finish. validation_report.csv records the exact number of valid 005
genes retained. Existing valid 005 genes are NOT refitted. Missing/invalid or
nonconverged results must be fitted. Saved quota-failure attempt receipts are
archived before their retry counters are cleared; successful files are kept.
Finished batch records are reconciled only after the no-live-jobs check.

The resume deliberately targets iteration 005, so the original cohort finishes
005 and stops before 006. It uses 298 workers, a five-minute refill interval,
max_iter=5000 and the same excluded nodes. All user jobs count against the
300-job budget. The existing simulation scheduler remains in place.

If a maintenance job hits 23 hours or fails, fix the reported issue and resume
THE SAME state/scratch paths. Do not launch the EM until maintenance finishes.
Cancelled dependent jobs will not restart themselves; use
the following explicit commands instead of rerunning submit_recovery.sh:

  # Only if retirement was interrupted:
  sbatch --exclude=midway2-0003,midway2-0011 --cpus-per-task=1 --mem=16G \
    --output="$STATE/logs/retire_%j.out" --error="$STATE/logs/retire_%j.err" \
    "$BUNDLE/storage_job" "$BUNDLE" retire "$PROJECT" \
    "$HOME/slider_retire_000_003_2026-10-28" --through 3 --apply

  # After retirement finished, or if retaining all historical results:
  sbatch --exclude=midway2-0003,midway2-0011 \
    --output="$STATE/logs/compact_%j.out" --error="$STATE/logs/compact_%j.err" \
    "$BUNDLE/storage_job" "$BUNDLE" compact "$PROJECT" "$STATE" "$SCRATCH" --workers 2

  # Only once STATE/COMPLETE.json exists and no other controller is queued:
  python3 "$BUNDLE/storage.py" resume "$PROJECT" "$STATE" --target 5

4. Check the 680 renamed expression IDs on RCC
---------------------------------------------
Once original 005 is complete and its controller has exited, this single job
can run ALL expansion steps in order. It stops on any preflight/init failure:

  BRANCH="$PROJECT/slider_expanded_2026-10-28"
  sbatch --exclude=midway2-0003,midway2-0011 \
    --output="$STATE/logs/expand_%j.out" --error="$STATE/logs/expand_%j.err" \
    "$BUNDLE/expansion_job" "$BUNDLE" all "$PROJECT" "$STATE" "$BRANCH"

Alternatively, use the separate steps below to inspect each stage first.

This preflight is read-only with respect to the analysis. It can run after
storage recovery while the remaining original 005 genes finish.

  sbatch --exclude=midway2-0003,midway2-0011 \
    --output="$STATE/logs/preflight_%j.out" --error="$STATE/logs/preflight_%j.err" \
    "$BUNDLE/expansion_job" "$BUNDLE" preflight "$PROJECT"

  cat "$BUNDLE/checks/summary.txt"
  cat "$BUNDLE/checks/genotype_summary.json"

Required: 680 exact expression/annotation matches, 0 failed mappings, 680
windows with variants, 0 windows without variants. These are targets, not
claimed RCC observations. Input presence does not guarantee every gene has
eligible tissues. This is slider-only expansion; SuSiE-mix is not changed.

5. Once ORIGINAL iteration 005 is complete, initialize a separate branch
-----------------------------------------------------------------------
  BRANCH="$PROJECT/slider_expanded_2026-10-28"
  sbatch --exclude=midway2-0003,midway2-0011 \
    --output="$STATE/logs/expand_init_%j.out" --error="$STATE/logs/expand_init_%j.err" \
    "$BUNDLE/expansion_job" "$BUNDLE" init "$PROJECT" "$STATE" "$BRANCH"

Initialization refuses an incomplete/stale source checkpoint and any active
slider jobs. It validates the source completion certificate and reuses the
17,617 existing in-scope result files via read-only use of filesystem symlinks.
No large copies, no refitting of those genes, and no change to their priors.
Their source checksums and provenance are recorded. Keep original iteration
005/results while the expanded seed links depend on it; retirement tools
refuse deletion of a registered source checkpoint.

  cat "$BRANCH/seed_report.csv"
  cat "$BRANCH/INITIALIZED.json"

Expected: 17,617 reused entries + 680 added expression IDs = 18,297 in-scope
entries (including any genes with no eligible tissues). The manifest also
retains the 851 frozen chromosome exclusions. The 680 IDs represent 679 old
symbols: the two C2orf61 IDs receive distinct output keys. CTAGE5 is selected
by its exact expression ID, separately from the existing MIA2 expression row;
this annotation collision remains documented in renamed_gene_map.csv.

6. Start the expanded slider EM
------------------------------
  python3 "$BUNDLE/expansion.py" start "$BRANCH" --target 15

The expanded branch has its own iteration_000, seeded at original 005's
FROZEN fitting prior. It initially fits only the 680 added IDs. Only when all
entries are validated does it pool original + new counts and perform the
unchanged tissue-wide M-step. From expanded iteration_001 onward every
eligible gene is refitted/warm-started under the updated prior. It never
recycles old posteriors as if they had been fitted under a new prior.
Errors/nonconvergence block advancement; no failed genes are silently dropped.

Target 15 means 15 additional prior updates after this expanded seed, following
the original five updates. It is configurable at submission. Scheduling is
one gene per fit job, 298 workers, five-minute refills, total user budget 300.
Jobs/logs are named slide_expanded_* to distinguish them from the original run.

Expanded prior history:
  $BRANCH/results_slide_prior_em/prior_history.csv
Original history snapshot carried into the expanded branch:
  $BRANCH/results_slide_prior_em/original_prior_history.csv
Parent iteration and frozen-prior hash:
  $BRANCH/results_slide_prior_em/seed_provenance.rds

Historical retirement later
---------------------------
Do not use rm -r on iteration folders. The storage tool can retire old full
results, repair the remaining summary dependencies and preserve all history.
It guards the latest incomplete iteration's warm-start checkpoint. A complete
latest checkpoint can release its predecessor; an incomplete one cannot.
Example: after original 005 completes, a NEW state path can retire through 004:

  python3 "$BUNDLE/storage.py" retire "$PROJECT" \
    "$HOME/slider_retire_004_2026-10-28" --through 4

This prints a plan only. Add --apply via storage_job to execute after reviewing
the listed scope. The operation does not remove original 005, and it does not
touch any expanded branch. Automatic old-iteration deletion is not enabled.

Local validation
----------------
See TEST_RESULTS.txt. Tests use synthetic data/files, not the RCC data. The
bundle performs actual content, history and source-version checks on RCC.
