# Simulation restart: rollback and data recovery

The version before the n=250 restart is:

`cc02614aa360640cb4c5274f6a2096751cbd4fdf` (22 October 2026).

The restart commit is named **Archive simulations and prepare n250 three-method benchmark**.
It preserves the active source changes, generated jobs, RCC update bundle,
frozen code snapshots, result inventory and tracked figures. Unrelated analysis
work is not part of that commit.

## What a GitHub push preserves

A push of the restart commit preserves its tracked files and reachable Git
history. It does **not** upload ignored files. In particular, these local archive
paths are excluded:

- `archive/simulations_before_n250_2026-10-23/simulation results/chunks/`
- `archive/simulations_before_n250_2026-10-23/simulation results/slide_prior_v1/`

The complete local result collection is 1,061 files, 38,551,712,829 bytes
(35.90 GiB). The root RCC update zip contains scripts and configuration only.
Neither a fresh clone nor that zip is a complete simulation-data backup.
Some September checkpoints exist in commit `b7078e1` on the local branch
`codex/backup-before-push-cleanup-20260914`. That commit is not an ancestor of
the current `main` branch, so pushing `main` does not upload that history.
It is not a backup of the complete current collection.

Before deleting the last RCC copy, copy the **entire** local archive folder
to independent storage and verify its contents against `results_inventory.csv`.
Retain the source data and software needed for refitting separately; compact
simulation records do not save the original X/y or full model fits.

## Roll back the code

Use `git log --oneline` to find the restart commit, then on a working tree with
your later work safely saved:

```sh
git revert <restart-commit-hash>
```

This adds a reversal commit and keeps shared Git history intact. It restores
tracked code and the former locations of tracked figures. It does not recreate
or relocate ignored checkpoint files. Conflicts with later edits must be resolved.

For a separate inspection of the previous tracked version, without changing
the current checkout, use:

```sh
git worktree add --detach ../susie_mix_before_n250 cc02614aa360640cb4c5274f6a2096751cbd4fdf
```

That checkout also needs the separately backed-up data to regenerate analyses.

## Restore the simulation data

Recover `simulation results/chunks/` and `simulation results/slide_prior_v1/`
from the local archive or its independent backup. Place them under the desired
project's `simulation results/` directory for the older scripts. Keep the new
`n250_r300_v1/` results separate. Alternatively, keep the archive intact and use
its `regenerate_figures.R` wrapper to analyze the archived files in place.

An archive stored on the same disk as the working project is not protection
against losing that disk. No independent data backup or GitHub push was performed
when preparing this restart.
