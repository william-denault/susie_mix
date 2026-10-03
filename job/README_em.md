# Run one or several EM iterations

After fitting, use the [dedicated EM analysis scripts](../script/analysis/README_em_analysis.md)
to generate summaries, descriptive results, and dominant one-CS lead-distance
comparisons from the last completed iteration.

On the cluster, from the project's `job` directory:

```bash
sbatch em_susie_mix
```

The preparation job loads `R/4.2.0`, estimates the tissue priors, saves them,
and submits one Slurm array using the resource settings in `test1`
(`broadwl`, 23 hours, 40 GB, one CPU per chunk). Every **new EM iteration uses
298 chunks**, independently of the baseline chunk layout. For 18,468 genes,
this gives 290 chunks of 62 genes and 8 chunks of 61 genes. The initial gene
universe comes from `data/temp_index/chunk_*_genes.txt`; later iterations use
the preceding iteration's gene manifest. Every gene is assigned exactly once,
in the same order, and the assignment is saved in the new `manifest.csv`.
Small test runs use at most one chunk per gene, without empty jobs. Baseline
gene lists and generated `run_chunk_*.R` scripts are not changed or executed.

Sync `em_utils.R` and `prepare_em_iteration.R` before preparing the next new
iteration; the existing Slurm launcher reads the resulting chunk IDs. Even if
the baseline used 185 chunks, EM will use 298. A resumed iteration retains its
original saved assignments and completion markers; it is never repartitioned
midway. The next new iteration uses the new layout.

To request 20 consecutive iterations after the baseline finishes:

```bash
sbatch em_susie_mix 20
```

Smaller chunks reduce the gene count per job from about 100 to 62. The actual
iteration time still depends on queue waits, available parallel slots, slow
genes, and prior-aggregation time; this does not guarantee 20 iterations per week.

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

- Without `results_em/prior_history.csv`, pool `susie_mix$alpha` from every
  gene result in `results`, separately by tissue and coding, for iteration 1.
- With an existing history, verify the latest iteration is complete, then
  pool `weighted_fit_mix$alpha` from its saved gene fits for the next update.
- Every component with exactly `V > 0` contributes its entire alpha row,
  including predictors outside any reported CS. Positive variances below
  `1e-9` still count; exactly zero-variance assignments are integrated out.
  CS purity, duplicate or overlapping CSs, association P-values, read counts,
  and an additional sample-size threshold do not select EM components.
- With all three coding classes present, normalize their summed active alpha
  masses. Unequal block sizes are handled by dividing each coding prior mass
  across its retained predictors. If a class is absent, optimize the
  conditional-coding objective rather than using naive pooled proportions.
- A fit with no active components contributes zero counts. If the whole
  tissue has no active components, retain its previous prior (or initialize
  uniformly if it has none). PIP sums remain descriptive diagnostics.
- The worker fits only weighted mixed-coding SuSiE and initializes from the
  preceding full posterior and variances using the new coding weights.
  Existing data QC is unchanged. The inner budget is `max_iter=1000`,
  `tol=1e-5`; `min_abs_corr=0.5` controls reported CSs only.

New rows use `update_method=susie_active_component_alpha_v2`, matching the
original unfiltered run. This restores the collapsed variational EM update
in [EM_PRIOR_REVIEW.md](EM_PRIOR_REVIEW.md). No `coding_prior_cs` metadata is
required or used. The former purity-filtered v3 update is no longer active.

```text
results_em/
  prior_history.csv             # Append-only iteration records
  iteration_012/
    priors.csv                 # Frozen weights used for iteration 12
    manifest.csv
    results/<gene>.rds          # Full fits for the next E/M transition
    completed/chunk_001.done
  iteration_013/
    priors.csv                 # New weights calculated from iteration 12
    manifest.csv               # New 298-chunk assignment
    source_audit.csv
    component_counts.rds
    results/<gene>.rds
    completed/chunk_001.csv
    completed/chunk_001.done
    logs/
```

`alpha_total` equals `n_active_components` up to rounding. History retains
`n_components`, `n_zero_variance_components`, fit counts, source errors,
`mstep_q_gain`, and `max_abs_prior_change`. `source_elbo_sum` includes all
valid fits, including all-null fits, when every fit supplies a finite ELBO.
Compare it only for the same data and model settings. Prior history records
weights used for each iteration; iteration 12's row was calculated from
iteration 11, not from the final iteration-12 fits.

## Continue the original unfiltered run after iteration 12

Keep the existing RCC `results/` and `results_em/` directories intact. There
is no baseline rerun, no new results layout, and no need to repeat iterations
1-12. The original complete history, iteration-12 manifest, completion
markers, and full gene `.rds` files must remain on RCC. The partial local
copy of `prior_history.csv` is not a replacement for the RCC history.

Sync these code files, preserving their relative paths:

- `script/scan_tissue_attempt/em_utils.R`
- `script/scan_tissue_attempt/prepare_em_iteration.R`
- `script/scan_tissue_attempt/run_em_chunk.R`
- `script/scan_tissue_attempt/workhorse_em.R`
- `script/scan_tissue_attempt/workhorse_utils.R`
- `script/scan_tissue_attempt/get_gene_annotations.R`
- `job/em_susie_mix`
- `job/em_susie_mix_array`

From the RCC project's `job/` directory, choose one continuation count:

```bash
sbatch em_susie_mix 5
```

This adds iterations **13-17**. Alternatively, `sbatch em_susie_mix 8` adds
iterations **13-20**. Submit only one of these commands. Each new iteration
uses 298 chunks. The existing launcher queues successive preparations after
successful arrays; the count is additional iterations, not the target index.
Do not use `resume` for a completed iteration 12: that option retries an
unfinished iteration with its already frozen priors.

The launcher checks previous completion and active jobs before preparing a
new iteration. If full history, source fits, or completion markers are
missing, reconcile the saved run; do not replace its history or bypass the
checks. Preparation accepts old v2 fits without CS metadata and never
changes their posteriors. Existing snapshots and old numeric history values
are preserved; new diagnostics may be appended as additional columns.

The requested count is a compute budget, not a convergence guarantee.
Inspect prior changes, source-fit ELBOs, and fit errors after the sequence.

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
not that all fits succeeded. Error records contribute no posterior counts and appear in
the next iteration's `source_audit.csv` and history error counts. These saved
gene/tissue errors do not stop an automatic sequence, just as they do not
block a manual next iteration. Review the CSVs, logs, and source audits when
assessing a sequence. Unconverged fits are different: a worker saves their
results and per-gene counts, then exits unsuccessfully without a completion
marker. Resume retries those genes. Preparation also rejects source fits that
are not confirmed converged. An entirely failed source scan is rejected.

Missing/unreadable gene files, mismatched predictor maps, invalid component
posteriors, or failed warm-start/weight checks stop
preparation. The initial results must match the full chunk gene list, so a
partially finished original scan cannot initialize the priors. A tissue with
no active component posteriors retains its previous prior and is marked
`carried_forward_no_active_components`. Without a previous prior, it starts
uniformly and is marked `initialized_uniform_no_active_components`; this is
a fallback, not a learned estimate. A fit with zero reported CSs is valid
and still contributes all of its positive-variance components. Full fits
retain all L rows for subsequent fitting.

All fitting workhorses use `min_abs_corr=0.5` for reporting only. No
association-P, read-count, lead-PIP, or extra sample-size filter is added to
the existing data QC. No pseudocount is added to the M-step. Disconnected
coding groups retain their previous relative total mass because the data
cannot identify it.

## Timing

New `completed/chunk_00N.done` RDS markers include `started_at`, `finished_at`
(UTC), and `elapsed_seconds`. The time from the earliest chunk start to the
latest finish measures the fine-mapping phase, including staggered chunk
starts. It excludes prior aggregation and the initial queue wait. Existing
markers without these additional fields remain valid for resuming/advancing.
For resumed chunks the marker times cover the final attempt; use Slurm
accounting to inspect earlier attempts as well. Slurm accounting can report
job queue, start, end, and elapsed times:

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
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_unfiltered_resume.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_marginal_likelihood.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_susie_smoke.R
bash script/scan_tissue_attempt/tests/test_em_launcher.sh
```
