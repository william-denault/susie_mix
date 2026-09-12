# Run one or several EM iterations

On the cluster, from the project's `job` directory:

```bash
sbatch em_susie_mix
```

The preparation job loads `R/4.2.0`, estimates the tissue priors, saves them,
and submits one Slurm array using the resource settings in `test1`
(`broadwl`, 23 hours, 40 GB, one CPU per chunk). It discovers the existing
`data/temp_index/chunk_*_genes.txt` lists: currently 185 chunks / 18,468 genes.
The original generated `run_chunk_*.R` scripts are not changed or executed.

With no argument, this runs one iteration. To run five consecutive iterations:

```bash
sbatch em_susie_mix 5
```

The count is the number of iterations to run, not the final iteration number.
For example, after iteration 2 completes, `5` runs iterations 3 through 7.
Each iteration estimates fresh tissue priors from all results of the preceding
iteration, appends `prior_history.csv`, and writes its own `iteration_00N/results`.
The next preparation is queued with a Slurm `afterok` dependency on the entire
current array and its preparation job. It starts only after both succeed;
waiting does not occupy a compute node. See the [Slurm dependency documentation](https://slurm.schedmd.com/sbatch.html#OPT_dependency).

The sequence stops after the requested count; it does not test convergence.
The preparation job finishing only means the array was submitted. Its output
gives the array and next preparation IDs. Check the weights after the sequence;
the last history rows are the priors used for the final fit, not a further
update calculated from that final fit.

## What each iteration uses

- No `results_em/prior_history.csv`: pool `susie_mix$pip` from **every gene
  result in `results`**, separately by tissue and coding, to prepare iteration 1.
- Existing history: read its latest iteration, verify all its chunks have
  finished and all gene files exist, then pool `weighted_fit_mix$pip` from that
  iteration to prepare the next one.
- Each coding prior is its pooled PIP sum divided by the sum across all three
  codings. Every predictor PIP contributes, including those outside credible
  sets; there is no PIP or credible-set filter and no gene-level averaging.
- `workhorse_em.R` runs only the weighted mixed-coding fit. It retains the
  original data processing, QC, and SuSiE settings, but omits additive-only,
  unweighted, and permutation fits and marginal association tests.
- Priors are matched to the existing `SMTS` tissue names. A class's mass is
  divided equally among its retained predictors; if a class is absent after
  filtering, weights are renormalized across the remaining classes.

```text
results_em/
  prior_history.csv
  last_array_job_id.txt
  last_continuation_job_id.txt  # Latest automatically queued preparation, if any
  iteration_001/
    priors.csv                 # Frozen priors used to fit this iteration
    manifest.csv               # Frozen chunk/gene assignments
    source_audit.csv            # Errors and nonconvergence in source results
    array_job_ids.txt
    continuation_job_ids.txt   # Next preparation IDs, when a sequence continues
    results/<gene>.rds
    completed/chunk_001.csv     # Per-gene status for this chunk
    completed/chunk_001.done    # Every gene in this chunk was attempted/saved
    logs/susie_mix_<job>_<task>.out
    logs/susie_mix_<job>_<task>.err
  iteration_002/
    ...
```

History rows contain `iteration`, `source_iteration`, `tissue`, `pi_add`,
`pi_rec`, `pi_dom`, the underlying PIP sums, fit/error counts, and a timestamp.
Thus iteration 1's row is estimated from the original scan and is the prior
**used for** iteration 1. Iteration 2's row is estimated from iteration 1.
Each gene RDS contains a list of successful tissues, with `weighted_fit_mix`,
the predictor map, coding priors/weights, and the existing fit metadata.

## Interrupted jobs and exceptional fits

After the array has stopped, retry unfinished chunks with:

```bash
sbatch em_susie_mix resume
```

This keeps the same iteration and priors, submits only unfinished chunks, and
reuses successful gene outputs already saved in those chunks. A failed array
submission can be retried this way too. Concurrent preparation is locked,
and a queued/running array blocks another launch.
An automatically queued preparation also blocks a competing manual launch.

A failed, cancelled, or timed-out array task prevents the next preparation
from running; `--kill-on-invalid-dep=yes` cancels that blocked preparation.
Fix the cause, then resume the interrupted iteration. To finish it and run
two more new iterations:

```bash
sbatch em_susie_mix resume 3
```

Here the count includes the resumed iteration. If the current iteration is
already complete, use the command without `resume`.

To stop automatic continuation while allowing the current array to finish:

```bash
scancel "$(cat ../results_em/last_continuation_job_id.txt)"
```

Use this from `job/` while that preparation is still pending. Its ID is also
printed in the preparation log. If continuation submission itself fails,
the already submitted array keeps running; the log gives the command for
launching the remaining iterations after it finishes.

As in the original scan, gene/tissue errors are saved as explicit records;
they do not stop other genes. Completion means all genes were attempted,
not that all fits succeeded. Error records contribute no PIPs and appear in
the next iteration's `source_audit.csv` and history error counts. These saved
gene/tissue errors do not stop an automatic sequence, just as they do not
block a manual next iteration. Review the CSVs, logs, and source audits when
assessing a sequence. Nonconverged fits contribute their
PIPs and are flagged in the audit. An entirely failed source scan is rejected.

Missing/unreadable gene files, mismatched predictor maps, or invalid PIPs stop
preparation. The initial results must match the full chunk gene list, so a
partially finished original scan cannot initialize the priors. A tissue with
zero total PIP has an undefined initial prior and stops preparation. In later
iterations its previous prior is retained and explicitly marked
`carried_forward_zero_pip`. There is no pseudocount or probability floor.

## Timing

New `completed/chunk_00N.done` RDS markers include `started_at`, `finished_at`
(UTC), and `elapsed_seconds`. The time from the earliest chunk start to the
latest finish measures the fine-mapping phase, including staggered chunk
starts. It excludes prior aggregation and the initial queue wait. Existing
markers without these additional fields remain valid for resuming/advancing.
Slurm accounting can report job queue, start, end, and elapsed times:

```bash
sacct -j "$(cat ../results_em/last_array_job_id.txt)" \
  --format=JobID,State,Submit,Start,End,Elapsed
```

The default cluster project path is `/project2/mstephens/wdenault/susie_mix`.
For a different checkout:

```bash
export SUSIE_MIX_PROJECT_DIR=/path/to/susie_mix
sbatch em_susie_mix
```

The workers need the same GTEx files, PLINK executable, and R packages
(`susieR`, `data.table`, `matrixStats`) as the original scan. No actual scan
is run by the local checks. Test the orchestration and prior calculations with:

```bash
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_iteration.R
bash script/scan_tissue_attempt/tests/test_em_launcher.sh
```
