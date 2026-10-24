# Recover the 20 checkpoint-mismatch failures

This utility applies only to original array `49451873`, tasks 16–19, 26–29,
42–47, 64–67 and 70–71. These correspond to simulation conditions 17–20,
27–30, 43–48, 65–68 and 71–72.

The simple recovery is to retain the conflicting checkpoints as backups and
regenerate these 20 conditions using the unchanged simulation code, seeds and
settings. It does not bypass the compatibility check or mix old and new records.
At most 6,000 datasets (18,000 method fits) are repeated. Other conditions keep
their existing results and continue normally.

## Run on RCC

Upload `recover_simulation_checkpoints_2026-10-24.zip` to the project directory:

```sh
cd /project2/mstephens/wdenault/susie_mix
unzip -o recover_simulation_checkpoints_2026-10-24.zip
python3 script/sim/recover_failed_checkpoints.py --apply
```

The archive adds only this recovery utility, its fixture tests and these
instructions. It does not replace frozen simulation sources or the scheduler.
Python 3.6 or later and the existing shared scheduler are required. No R package
installation or job regeneration is needed. For a read-only preview, omit
`--apply`.

The command checks all 20 Slurm accounting states and error logs. Each must
still show FAILED, exit 1:0 and the checkpoint-mismatch error; a different failure
stops the operation before any checkpoint is moved. It also checks frozen source
hashes, scheduler ownership and existing assignments.

Only the pending shared controller is briefly held while its condition list is
updated. Running simulation and EM workers keep running. The controller is then
released and schedules retries within its existing 300-job budget. Do not submit
a second controller or another array. If the controller is already running its
short check, retry the same command in five minutes.

## Preserved state and interruption recovery

Old checkpoints and any matching interrupted `.tmp` saves go into:

`simulation results/n250_r300_v1/checkpoint_backups/array_49451873/`

The original scheduler configuration and move inventory are retained in:

`results_slide_prior_em/shared_scheduler/checkpoint_recovery_49451873.json`

The command can be repeated. It resumes an interrupted archive operation and
does not archive new results after the conditions have been registered. If the
process is forcibly interrupted while holding the controller, rerun the same
command to finish the operation and release its hold. It will not release a
controller hold that was already present before this recovery.

New results return to the usual `n250_r300_v1/chunks/` paths. A `REGISTERED`
message means recovery is queued through the existing scheduler, not that the
20 simulations have finished. The four zero-variance replicate errors in
conditions 21 and 24 are separate and are not changed by this utility. The MAF
rule also remains unchanged: filtering in the full donor pool before n=250
subsampling.

## Local validation

Ten fixture tests passed with Slurm mocked, including scheduler integration,
backup preservation, interruption recovery, repeated commands, rejection of
other failure causes, active-job protection and preservation of user holds:

```sh
python3 -B script/sim/tests/test_checkpoint_recovery.py
```
