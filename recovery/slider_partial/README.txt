Stop bulk compaction safely and resume existing iteration 005
============================================================

Full compaction of iteration 004 is optional. Once enough project quota is
available, the unconverted 004 files remain valid warm starts. The required
handoff is to resolve any interrupted rewrite, update cached validation for
the files already changed, and install the saved-call storage fix for new fits.

This bundle does NOT compact the remaining thousands of files, change EM
priors/settings, add genes, or prepare iteration 006. It resumes the original
005 in the existing results_slide_prior_em folder, with target 005 only.
The later same-folder 006--040 expansion remains a separate task.

What the command does
---------------------
1. Verify the new bundle and original maintenance identity.
2. Read the original queued_jobs.tsv. For the submitted run these identify
   compaction 49638385 and pending resume 49638386. Validate their live names
   and states; refuse competing slider jobs or an already-running resume.
3. Submit a finalization job ON HOLD, dependent on those jobs ending.
4. Cancel only the recorded pending resume and compaction. Simulation and
   unrelated jobs are not cancelled. Release the replacement finalizer.
   It cannot start until the old jobs have terminated and released the lock.
5. Resolve the small set of in-flight file transactions:
     START: retain the checksum-verified original; discard incomplete scratch output;
     PREPARED/VACATED/REPLACED: finish the already verified replacement;
     committed: retain the existing result and its durable receipt.
   No new full-file compaction is initiated.
6. Check committed result stamps and the complete original inventory. Retain
   all untouched files exactly as they are. Their size/mtime must match the
   original maintenance snapshot; existing certified hashes/summaries are
   reused under the same validation policy as the existing parallel scheduler.
   Untouched multi-terabyte fits are NOT all read and checksummed again.
7. Repair result hashes, dependency stamps and completion certificates for
   changed files. All priors/history/settings/counts are checked unchanged.
   Saved 005 files without a valid cached summary are deferred to the existing
   parallel validators instead of being read serially by this coordinator.
8. Install the already tested workhorse change that removes only fit$call
   from warm starts and new outputs; all scientific fitting fields remain.
9. Submit the original controller to finish 005: 298 maximum workers,
   five-minute refills, total user-job budget 300, frozen settings unchanged.
   Valid saved genes are reused; missing/invalid/nonconverged genes need fitting.

This avoids the remaining bulk conversion; exact runtime still depends on
queue wait, metadata I/O, and any in-flight large files that need recovery.
Check rcchelp quota for actual group headroom: directory size is not free quota.

Upload and run once
-------------------
Upload slider_partial_rcc_2026-10-29.zip to /home/wdenault, then:

  cd "$HOME" &&
  unzip slider_partial_rcc_2026-10-29.zip &&
  cd "$HOME/slider_partial_2026-10-29" &&
  sha256sum -c SHA256SUMS &&
  bash stop_compaction_resume_005.sh

Defaults:
  PROJECT=/project2/mstephens/wdenault/susie_mix
  ORIGINAL_STATE=/home/wdenault/slider_storage_state_2026-10-28
  SCRATCH=/scratch/midway2/wdenault/slider_storage_2026-10-28
  FINAL_STATE=/home/wdenault/slider_partial_finalize_2026-10-29

If those original paths were customized, pass the actual paths:

  bash stop_compaction_resume_005.sh PROJECT ORIGINAL_STATE SCRATCH

Do not manually cancel jobs first, delete iteration folders, or erase the
original state/scratch. The handoff arranges its replacement before stopping
the old compaction. Do not launch the earlier separate-folder expansion bundle.

Progress
--------
  STATE="$HOME/slider_partial_finalize_2026-10-29"
  cat "$STATE/handoff.json"
  squeue -u "$USER" -o "%.18i %.26j %.12T %R"
  tail -n 30 "$STATE"/logs/finalize_*.out "$STATE"/logs/finalize_*.err

After finalization:
  cat "$STATE/COMPLETE.json"
  cat "$STATE/validation_report.csv"
  cat "$STATE/resume_job.json"

The new COMPLETE.json explicitly reports mode=partial_compaction_finalized,
metadata_finalized=true, and compaction_complete=false. It certifies a valid
mixed original/compacted checkpoint, not that every file was converted.
The ORIGINAL state is not falsely marked fully compacted; it gets a pointer
PARTIAL_FINALIZED.json. Use the NEW final state for future recovery decisions.

After handoff, original 004 is still needed until 005 completes. The untouched
004 files remain large, but new successful fits omit the duplicated saved call.
Keep all prior-history folders. Later rolling retention can retire old full
posteriors once their warm-start role ends.

If interrupted or stopped
-------------------------
All original receipts and scratch transactions are preserved and resumable.
If submission/cancellation fails, the replacement may remain held. Inspect
handoff.json, handoff.pending.json and squeue before doing anything else.
Do not blindly rerun the submitter or delete a submission marker: an accepted
job with an ambiguous response must not be duplicated.

If the finalizer itself fails after the old jobs are terminal, inspect its log.
After resolving the issue and confirming no other maintenance/slider job is
active, rerun just the finalization job with the SAME paths:

  BUNDLE="$HOME/slider_partial_2026-10-29"
  STATE="$HOME/slider_partial_finalize_2026-10-29"
  sbatch --exclude=midway2-0003,midway2-0011 \
    --output="$STATE/logs/finalize_%j.out" --error="$STATE/logs/finalize_%j.err" \
    "$BUNDLE/partial_job" "$BUNDLE" \
    /project2/mstephens/wdenault/susie_mix \
    "$HOME/slider_storage_state_2026-10-28" \
    "/scratch/midway2/$USER/slider_storage_2026-10-28" "$STATE"

A durable resume_job.json prevents the helper from submitting the controller
twice. If that controller later fails, diagnose it using its own logs and
the existing resume mechanism; rerunning the finalizer is not an EM restart.

Validation
----------
See TEST_RESULTS.txt. Synthetic tests exercise interrupted file states, a
checkpoint mixing converted and untouched files, historical retirement,
metadata repair, and parallel validation of a saved gene without refitting.
Actual RCC data and Slurm execution are guarded at runtime, not simulated
as a claim of having already executed the handoff remotely.

Slurm behavior used here:
https://slurm.schedmd.com/sbatch.html (--hold and afterany dependencies)
https://slurm.schedmd.com/scancel.html (--state=PENDING and job termination)
