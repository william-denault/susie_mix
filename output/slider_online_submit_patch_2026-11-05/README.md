# Recover the accepted Slurm submission for online batch 3

The first two online batches committed successfully. During batch 3, controller
49853687 stopped waiting for sbatch after 10 seconds, leaving an unresolved
submission journal. Slurm nevertheless accepted array 49853903; accounting
reported all 246 workers COMPLETED with exit code 0:0. The stopped successor is
49853902. The last reported state contains 4,837 validated genes: 2,030 + 2,030
in the committed batches, and 777 in batch 3.

This repair is ONLY for this specific stopped state and the already installed
seed-retention patch. It verifies accounting and existing results, restores the
missing submission receipt for task_000038, archives the matching ERROR and
pending submission journal, and resumes the unfinished part of batch 3.

## Run on RCC

Upload slider_online_submit_patch_rcc_2026-11-05.zip to /home/wdenault:

```bash
cd /home/wdenault &&
unzip slider_online_submit_patch_rcc_2026-11-05.zip &&
cd /home/wdenault/slider_online_submit_patch_2026-11-05 &&
sha256sum -c SHA256SUMS &&
bash start_repair.sh
```

The last line prints a Slurm repair-job ID. You may disconnect after submission.
The repair runs on a compute node, then submits the normal controller itself.
No gene workers are cancelled, no completed fit is resubmitted, and neither
committed prior update is replayed. No RDS files are copied or recompressed.

```bash
tail -n 40 /home/wdenault/slider_online_submit_patch_2026-11-05/logs/repair_*.out
tail -n 40 /home/wdenault/slider_online_submit_patch_2026-11-05/logs/repair_*.err
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/code/status.sh
```

Expected repair output:

```
Recovered array 49853903: 246 completed workers; 4837 validated results preserved; two updates unchanged.
Resuming iteration 015, batch 3; no completed gene is resubmitted.
Submitted controller JOBID
```

Do not use the old HOME launch bundles after installation. The patched installed
workflow has a different code fingerprint. Use its installed code/status.sh and,
if a later diagnosed failure needs a restart, its code/start.sh.

## Scope and checks

Only the runtime's minibatch.py and its checksum manifest change. Submission
acknowledgement waits up to 120 seconds instead of 10; the safe window before
another controller becomes eligible grows to 300 seconds. The 23-hour job
limits, one-minute refill eligibility, 298-worker cap and shared 300-job budget
are unchanged. No R code, numerical method, batch membership, gene selection,
rho, priors, counts or history is changed. Online 015-020 and full 021-025 remain
the schedule. The seed-retention fix remains installed.

The repair checks that all 246 accounting records belong to the expected user,
job name and working directory and are COMPLETED/0:0. It also checks their log
filenames and validated receipts. Other active slider jobs, changed inputs,
missing receipts, different errors or a different committed stage block repair.
Large result files are checked using their existing receipt sizes/timestamps;
small analysis metadata, commits, histories and receipts are hashed before and
after. Task histories and retry counts are retained, including the recovered
array. Current successful results are never discarded.

Backups, accounting evidence, the installation transaction and completion
receipt are under the workflow's repairs/submission_ack_2026-11-05 directory.
Installation can resume after an interruption. Repeating a completed repair
does not reset progress or submit a second controller when one is active.

A future submission that remains ambiguous after 120 seconds still stops with
its journal intact; the code does not blindly retry an unacknowledged sbatch.
This avoids duplicate arrays and bypassing the shared slot budget.

The Slurm distinction between submission acknowledgement and eventual job
execution is documented at https://slurm.schedmd.com/sbatch.html . Per-task
accounting is documented at https://slurm.schedmd.com/sacct.html .

## Validation

All 13 new repair tests and the 20 existing runtime Python tests passed locally.
The existing Windows hard-link cleanup test required an unsandboxed rerun.
All R files are byte-for-byte identical to the installed seed-retention repair.

The included tests exercise real receipt validation and controller selection,
preservation of the two committed updates and saved results, recovered attempt
accounting, interrupted installation, repeated installation after progress,
rejection of changed/incomplete/live jobs, the 120-second deadline and the
corresponding safety window. They also verify that another ambiguous timeout
retains the submission journal and is not blindly retried. Slurm calls are
mocked in local tests; the installer repeats the checks against RCC at runtime.
