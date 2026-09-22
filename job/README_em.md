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

- No `results_em/prior_history.csv`: pool `susie_mix$alpha` from **every gene
  result in `results`**, separately by tissue and coding, to prepare iteration 1.
- Existing history: read its latest iteration, verify all its chunks have
  finished and all gene files exist, then pool `weighted_fit_mix$alpha` from that
  iteration to prepare the next one.
- The empirical Bayes update uses the expected coding assignments of SuSiE's
  **active components (`V > 0`)**. Exactly zero-variance components are integrated
  out of the coding-prior M-step. With all three codings present in every fit,
  each prior equals its summed active `alpha` divided by the number of active
  components. All predictors in those rows contribute, even outside credible
  sets. A fit whose variances are all zero contributes no assignment counts.
  Positive variances below `1e-9` still contribute: activity is defined by exact
  zero, not SuSiE's numerical PIP/CS reporting tolerance. Marginal PIP sums are
  saved as full-fit descriptive diagnostics; they do not determine the priors.
- `workhorse_em.R` runs only the weighted mixed-coding fit. It retains the
  original data processing and QC, but omits additive-only, unweighted, and
  permutation fits and marginal association tests. Each fit starts from its
  predecessor's posterior and variance estimates, with the NEW prior weights.
  The inner fitting budget is `max_iter=1000`, `tol=1e-5`.
- Priors are matched to the existing `SMTS` tissue names. A class's mass is
  divided equally among its retained predictors; if a class is absent after
  filtering, weights are renormalized across the remaining classes. In that
  case the M-step optimizes the corresponding conditional-coding objective,
  rather than incorrectly using a global normalized count.

This is collapsed variational EM for the coding-mixture SuSiE model, targeting a lower
bound on the sum of per-gene log marginal likelihoods. It is not the cTWAS
Bernoulli-prior model. The model, derivation, and comparison with the paper
are explained in [EM_PRIOR_REVIEW.md](EM_PRIOR_REVIEW.md).

```text
results_em/
  prior_history.csv
  last_array_job_id.txt
  last_continuation_job_id.txt  # Latest automatically queued preparation, if any
  iteration_001/
    priors.csv                 # Frozen priors used to fit this iteration
    manifest.csv               # Frozen chunk/gene assignments
    source_audit.csv            # Errors and nonconvergence in source results
    component_counts.rds       # Active alpha sums by tissue and available coding classes
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
`pi_rec`, `pi_dom`, active component `alpha` sums, descriptive PIP sums, fit/error
counts, and a timestamp. `n_components` is the total number of component rows;
`n_active_components` counts the rows used, and `n_zero_variance_components`
counts those excluded. `alpha_total` equals `n_active_components` up to rounding.
`n_fits` counts all valid fits, while `n_active_fits` counts fits with at least
one positive-variance component. `mstep_q_gain` checks that the M-step increases its
objective; `max_abs_prior_change` tracks movement of the three coding weights.
`source_elbo_sum` records the summed source-fit ELBO over **all valid fits,
including all-null fits**, only when all included
fits supply a finite value (`n_source_elbo` gives the count). Compare these
values only across iterations with the same gene/tissue fits and model settings.
Thus iteration 1's row is estimated from the original scan and is the prior
**used for** iteration 1. Iteration 2's row is estimated from iteration 1.
Each gene RDS contains a list of successful tissues, with `weighted_fit_mix`,
the predictor map, coding priors/weights, and the existing fit metadata.

## Continuing a run prepared with an earlier update

The next new iteration uses the corrected update automatically, with the
completed fits as initialization. Full fits containing `alpha` and `V` must
be present on the cluster; completion logs or PIPs alone are insufficient.
The old history values and frozen iteration files are preserved. New rows use
`update_method=susie_active_component_alpha_v2`. Existing
`susie_component_alpha_v1` rows retain their tag and values; they used all
component assignments in an uncollapsed variational EM update. Older untagged
rows receive `legacy_pip_share`, identifying the original heuristic PIP-share
update. Newly added diagnostic columns are empty for old rows. Only the
PIP-share method must not be described as marginal-likelihood EM. The v2 change
integrates out exactly inactive assignments and avoids their damping effect.
There is no need to delete completed fits or reset the history. A resumed
iteration keeps its frozen priors; the next newly prepared iteration uses v2.

Sync `em_utils.R`, `prepare_em_iteration.R`, `run_em_chunk.R`, and
`workhorse_em.R` together before launching the next iteration. If an automatic
continuation is already queued, cancel that pending preparation and allow the
current array to finish before replacing worker scripts and submitting again.

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
no active component posteriors in a later iteration retains its previous prior
and is marked `carried_forward_no_active_components`. If no previous prior
exists, it starts uniformly and is marked `initialized_uniform_no_active_components`;
this is a fallback, not a learned estimate. Only exactly zero-variance rows are
excluded from prior learning; full fits retain all L rows for subsequent fitting.
No association-P, read-count, lead-PIP, or CS filter is added to the existing data QC;
`min_abs_corr` remains zero. Positive V means active in the fitted model, not a
guaranteed real signal. No pseudocount is added to the M-step. Disconnected coding groups retain their previous
relative total mass because the data cannot identify it.

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
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_marginal_likelihood.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_susie_smoke.R
bash script/scan_tissue_attempt/tests/test_em_launcher.sh
```
