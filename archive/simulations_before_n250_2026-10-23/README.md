# Simulation archive before the n=250 restart

Created 23 October 2026 before changing any active simulation sources.
The original `simulation results/` directory was moved here, preserving
1,061 files and 35.90 GiB. `results_inventory.csv` records relative paths,
sizes and modification times. Counts and total bytes were verified after
the move. These files are retained locally, not uploaded to RCC.

**A GitHub push is not a backup of the complete simulation archive.**
The raw `simulation results/chunks/` directory and the complete
`simulation results/slide_prior_v1/` directory inside this archive are excluded
by the root `.gitignore`. A fresh clone will not contain those archived files.
Git records the code snapshots, inventory, publication snapshot and legacy
figures, but the full 35.90 GiB result collection requires a separate backup.
Copy this entire archive to independent storage and verify the copy before
removing the last RCC copy. See `../../SIMULATION_ROLLBACK.md`.

- `simulation results/chunks` and `figures`: legacy two-method results matching
  the September screenshots, n=500, 400 requested replicates, PVE 10–40%.
- `simulation results/slide_prior_v1`: newer five-method results and summaries.
- `script/sim`, `job`, shared genotype helpers, gene list and root entry scripts:
  snapshots of the active code immediately before the restart.
- `publication_plots`: a complete copy of the later publication figure workflow
  and its input summaries.
- `legacy_source/` and `legacy_source_2026-09-10.zip`: historical simulation,
  plotting and shared helper sources from commit 0218f96.
- `source_commit.txt`, `source_status.txt`: version-control provenance.
- `coverage_history_audit`: the earlier old/current coverage investigation.

To regenerate the newer archived analysis with its frozen scripts:

```r
source("archive/simulations_before_n250_2026-10-23/regenerate_figures.R")
```

For the original September panel arrangement:

```r
options(susie.archive.dataset = "legacy")
source("archive/simulations_before_n250_2026-10-23/regenerate_figures.R")
options(susie.archive.dataset = NULL)
```

These commands write figures inside this archive. They do not refit models.
Use this wrapper rather than the frozen root entry scripts, which retain
their original absolute-path/entry-point behavior. Exact historical refitting
still requires the original RCC genotype pool and software; old compact
checkpoints did not record package versions or full X/y. No absent inputs
have been reconstructed or claimed to be preserved.

New work lives outside this archive and writes to
`simulation results/n250_r300_v1/` with schema 4.
