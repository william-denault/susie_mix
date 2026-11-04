GTEx renamed-gene preflight and FULL project backup
================================================

This is the checking/backup stage, before expanding either EM cohort. It does
not install changes into the active script or job directories, update any
prior, cancel any job, or remove any result. No expanded EM is launched here.

The map has 680 available GTEx v8 expression IDs under 679 old names. C2orf61
has two IDs and receives two filenames. The available HIST1H3D ID is used; its
other reference ID has no expression row. CTAGE5 is matched to its expression
ID ENSG00000150527, separately from the old MIA2 row ENSG00000150526. These are
distinct GTEx v26 expression phenotypes even though the newer annotation
merges/renames their locus. Its collision is explicitly flagged in the map.

1. Unpack into a NEW directory under the project (no unzip -o).

   cd /project2/mstephens/wdenault/susie_mix
   test ! -e renamed_gene_preflight_2026-10-26 && \
     unzip renamed_gene_preflight_rcc_2026-10-26.zip
   cd renamed_gene_preflight_2026-10-26
   sha256sum -c SHA256SUMS
   BUNDLE="$PWD"
   PROJECT=/project2/mstephens/wdenault/susie_mix

2. Submit one preflight job on a compute node. The metadata check reads only
   expression row IDs/names, not model fits. The BIM check streams the genotype
   variant list once and uses little RAM. It may take several minutes because
   the source files are large. Submit checking jobs one at a time. They count
   as ordinary jobs in the existing shared scheduler's user-wide queue count.
   If Slurm reports a job-limit error, let a slot open and retry; cancel nothing.

   sbatch --exclude=midway2-0003,midway2-0011 \
     "$BUNDLE/preflight_job" "$BUNDLE" "$PROJECT"

   Expected checks/summary.txt:
     Candidate symbols: 679
     Candidate expression IDs: 680
     Exact expression and annotation matches: 680
     Failed/ambiguous mappings: 0

   Expected checks/genotype_summary.json:
     candidate_ids: 680
     windows_with_variants: 680
     windows_without_variants: 0

   These are expected targets, NOT claimed RCC results. Inspect the actual
   outputs. mapping_checks.csv and genotype_windows.csv list every candidate.
   A nonzero exit stops the check and retains the reports. No production gene
   list is edited. Variant presence is before SNP-only, donor-specific MAF,
   HWE and tissue QC; these checks cannot promise 680 successful fitted genes.

3. Back up the ENTIRE susie_mix project before any future change.
   This includes every iteration, simulations, logs, scripts, data within the
   project and hidden files. First check size, filesystem space and RCC quota:

   du -sh "$PROJECT"
   df -h /project2/mstephens/wdenault

   BACKUP_PARENT=/project2/mstephens/wdenault/susie_mix_backups
   mkdir -p "$BACKUP_PARENT/logs"
   BACKUP="$BACKUP_PARENT/full_$(date -u +%Y%m%dT%H%M%SZ)"
   printf '%s\n' "$BACKUP" > "$BUNDLE/full_backup_destination.txt"
   sbatch --exclude=midway2-0003,midway2-0011 \
     --output="$BACKUP_PARENT/logs/full_backup_%j.out" \
     --error="$BACKUP_PARENT/logs/full_backup_%j.err" \
     "$BUNDLE/backup_job" "$BUNDLE" "$PROJECT" "$BACKUP"

   The copied project is $BACKUP/susie_mix/. The script refuses an existing
   destination unless --resume is supplied for an interrupted/live copy.
   It uses rsync twice, then checks regular-file contents with SHA-256 and
   checks symbolic-link targets. It writes:

     $BACKUP/backup_metadata/backup_info.json
     $BACKUP/backup_metadata/verification.json
     $BACKUP/backup_metadata/file_checksums.jsonl

   If all current source entries match and no changes are observed during
   verification, it writes backup_metadata/BACKUP_VERIFIED.txt. Otherwise it
   writes LIVE_COPY_REPORT.txt and identifies missing, changed or extra files.
   The slider run is never paused. This is a live copy, not an atomic snapshot.
   Resume a partial/live copy using the SAME saved destination:

   BACKUP=$(cat "$BUNDLE/full_backup_destination.txt")
   sbatch --exclude=midway2-0003,midway2-0011 \
     --output="$(dirname "$BACKUP")/logs/full_backup_%j.out" \
     --error="$(dirname "$BACKUP")/logs/full_backup_%j.err" \
     "$BUNDLE/backup_job" "$BUNDLE" "$PROJECT" "$BACKUP" --resume

   Verified backups are protected against resuming/overwriting. No files are
   deleted. Files deleted from the live source can therefore remain in the
   backup and are reported as backup-only entries. Interrupted rsync is
   resumable; checksum records for unchanged files are also reused. Keep its
   logs outside the source project, as above, to avoid copying the backup's
   own changing log. A timed-out job is not a verified backup. Full copying and
   checksum verification may take many hours and need roughly the full
   project's size in additional quota. Filesystem free space is not user quota.

   Symbolic links are preserved as links, rather than following paths outside
   the project. External GTEx files or R libraries are not copied unless they
   are regular files/directories physically inside susie_mix. A backup on the
   same RCC filesystem protects against project edits/deletion, not filesystem
   failure. Preserve the backup directory when changing the source project.

4. Optional actual-fitter pilot, after the metadata check passes.
   Six tasks test AAED1, both C2orf61 IDs, CTAGE5, HIST1H3D and WARS. Each fits
   both models at a frozen existing prior (mix 018 and the latest completed
   slider checkpoint observed by preflight). No M-step/history change occurs.
   Source paths are in checks/pilot_checkpoints.csv. Outputs go under pilot/.
   To avoid submitting six extra jobs together, submit one task at a time:

   sbatch --array=1 --exclude=midway2-0003,midway2-0011 \
     "$BUNDLE/pilot_job" "$BUNDLE"

   Repeat with --array=2 through --array=6 as slots become available. Each
   pilot/<gene>/summary.csv reports PASS, NO_ELIGIBLE_TISSUES, or FAIL separately
   for both models. A PASS requires converged fits under the requested prior.
   Re-running a successful task validates and reuses its saved result.
   No plots are made. Graphics-device warnings from an interactive R session
   are unrelated to PLINK input errors.

Why these are pilots rather than the final expanded EM:
   Existing gene fits can be reused only at the same frozen prior. Backfilling
   the added IDs completes an expanded E-step at that prior. Subsequent prior
   updates require refitting all retained genes, using previous results as
   warm starts. A fresh results directory and prior history will keep that
   expanded analysis separate from the active slider run.

Return/sync checks/summary.txt, checks/genotype_summary.json, the CSV reports,
and the backup status/report before finalizing the expanded cohort.

Rsync option/exit-code reference: https://download.samba.org/pub/rsync/rsync.1
Exit 24 (vanished live source files) proceeds to the second pass and validation.
Other nonzero copy exits stop and retain the partial backup for inspection.
