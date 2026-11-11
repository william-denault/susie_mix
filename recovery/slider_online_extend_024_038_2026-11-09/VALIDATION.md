# Validation record — 2026-11-09

All checks run locally on synthetic fixtures or mocked Slurm. No RCC jobs were
submitted during development. Live RCC evidence was supplied by the user.

- Installed baseline minibatch.py and engine.R hashes exactly match the user’s
  RCC checksums. All vendored scientific R helper files are byte-identical to
  that baseline. Worker requests remain 48 GB; controllers remain 16 GB.
- Full 18,263-gene planning checks: unchanged 93-stage prefix through 023,
  exactly 61 balanced batches for each 024–033 scan (24 x 300, 37 x 299), no
  duplicated/missing genes, fixed random batch membership before runtime ordering,
  rho 0.8/0.5/0.3 at the requested boundaries, and full 034–038, 708 stages total.
- 41 upgrade/receipt/ordering tests: idle and live installation, no cancelled
  workers, interrupted code/config/plan/archive transactions, controller release,
  started-024 refusal, completed-023 handling, resume identity, completion beyond
  the old 95-stage limit, exact accepted-array recovery, invalid or changed
  receipts rejected, delayed accounting never treated as rejection, lost sbatch
  acknowledgement recovered without a second submission, longest-first ordering.
- 26 existing scheduler regressions: queue retries, handoff, ownership, cleanup,
  job-budget accounting, hard-link snapshots, and Linux-readable manifests.
- 96 actual susieRSlidePrior fits over a four-gene synthetic cohort, all iterations
  015–038 with reduced batch counts: exact sufficient-statistic updates, rho
  transitions, online/full transitions, previous-batch prior validation, and
  actual current/previous-result cleanup. These are local model tests, not GTEx
  performance or convergence benchmarks.
- Larger transition test: 61 synthetic genes, 186 updates through 038 including
  one full 61-batch online scan; 1,464 controlled gene fits; literal T/F names;
  zero-mass tissue; forced/inactive component exclusions; upgrade after the first
  023 gene is saved; existing output adoption without refit; identical frozen
  input/state; previous-prior warm starts; exact running-count math; interruption
  after durable commits; interrupted cleanup and orphan temporary-file cleanup;
  retained histories and only current/previous result files through 038.
- History checks: stable old probability columns; iteration/batch/phase/rho and
  all slider values preserved; input versus after_update roles distinguished;
  full-EM batch zero; completed/running/planned iteration summary; repeatable
  publication without duplicate updates. A separate summary test also checks
  original-then-runtime ordering in 023 and longest-first ordering in future scans.
- Final archive: ZIP integrity, all nested SHA256SUMS, LF line endings, Python
  3.6 syntax, and shell syntax. No fitted results or large data files packaged.

Reported RCC failure: controller 50013914 got “Unexpected message received” from
sbatch. Task 001411 was accepted as array 50014103; all 130 workers completed 0:0.
The installer rechecks that evidence and their saved receipts before clearing
only the matching failure journal. The user’s sacct check reported a maximum
20.43 GB peak across the controller/array batch records, with 130 COMPLETED and
one FAILED controller; this is consistent with the observed communication error,
not an OUT_OF_MEMORY stop. It does not certify memory headroom for every gene.
