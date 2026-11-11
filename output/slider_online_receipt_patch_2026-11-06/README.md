# Repair literal gene-name parsing at the batch-3 commit

The GTEx label `T` is the older symbol for TBXT (ENSG00000164458), already
recorded in renamed_gene_map.csv. It is a real gene, not a random batch label.
All 2,029 batch-3 workers completed validation, but the R commit reader inferred
logical TRUE from the one-row CSV gene field "T". Controller 49856781 therefore
reported `Missing/changed validated summary: T` before updating the prior.
The first two updates remain committed; 6,089 validated fits are saved.

This repair targets that exact stopped state, after the submission-timeout
patch. It changes only success-receipt reading in engine.R: gene names, hashes
and paths are explicitly character columns; stage and fit count are integers.
Existing CSV receipts and RDS fits are never rewritten. TBXT keeps its existing
T filename and Ensembl mapping. Batch membership, order, random seeds, fitting
settings, rho=0.3, pooled-count calculations and history conventions are unchanged.

## Run on RCC

Upload slider_online_receipt_patch_rcc_2026-11-06.zip to /home/wdenault:

```bash
cd /home/wdenault &&
unzip slider_online_receipt_patch_rcc_2026-11-06.zip &&
cd /home/wdenault/slider_online_receipt_patch_2026-11-06 &&
sha256sum -c SHA256SUMS &&
bash start_repair.sh
```

The final command submits one Slurm repair job. You may close the terminal once
its job ID is printed. The repair checks existing receipts, performs the real
batch-3 commit calculation with all writes redirected into memory, installs the
reader fix, archives the matching error, and submits the normal controller.
That controller commits batch 3 using the saved summaries and moves to batch 4.
There is no need to refit the 6,089 already validated genes.

```bash
tail -n 40 /home/wdenault/slider_online_receipt_patch_2026-11-06/logs/repair_*.out
tail -n 40 /home/wdenault/slider_online_receipt_patch_2026-11-06/logs/repair_*.err
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
```

Expected repair output includes:

```
PASS: literal gene T preserved; 2029 summaries validated; batch 3 dry-run commit succeeded.
Receipt repair verified: 6089 validated results preserved; two prior updates unchanged.
Submitted controller JOBID
```

The dry-run calculation may also print `Committed stage 3`, but it explicitly
runs in memory only. The actual on-disk commit occurs in the resumed controller.
Use status.sh to confirm `Committed stages: 3/59` and advancement to batch_04.
The schedule remains six online passes 015-020 and five full steps 021-025.
The seed-retention fix and 120-second submission response deadline are retained.

## Preservation and recovery

No simulation or fitting job is cancelled. The installer refuses active slider
jobs, an ambiguous pending submission, another error, a different stage,
changed code or missing/changed validated outputs. It preserves all task and
attempt records, fits, summaries, completed commits, priors and history. Metadata
are hashed before/after and large fits are checked against receipt file stamps.
The R preflight validates every batch-3 summary checksum and its cohort, then
computes the proposed update without writing analysis files.

Backups, before/after manifests, the preflight prior/history previews, the error
archive and completion receipt are saved under:

```
workflows/slider_online_015_025/repairs/literal_gene_receipts_2026-11-06/
```

The installation is resumable after an interruption. Repeating a completed
repair does not reset later progress or duplicate an already queued controller.
After installation, use the installed workflow's code/start.sh and code/status.sh;
the older HOME launch bundles have a different checksum fingerprint.

## Local validation

The R regression first reproduces the exact original failure with T in batch 3,
then completes all 54 minibatch updates and five full EM steps with T and F in
the cohort. It checks the prior calculations against independently pooled counts,
history numbering, repeated commits and preservation of the source fits. The
production dry-run preflight succeeds without changing fits, receipts or state.
Additional cases preserve literal T, F, TRUE, FALSE, NA and numeric-looking
identifiers, including all-digit checksums. Nine Python repair tests cover state
preservation, failure guards, interrupted installation and a resumed controller
that commits batch 3 without submitting any gene workers.

Only engine.R and its runtime checksum manifest differ from the previously
installed submission patch. The scientific helper and fitting files are identical.
Tests use local fixtures; the production preflight repeats the relevant checks
on RCC before clearing the specific recorded error.
