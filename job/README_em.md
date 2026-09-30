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

- No `results_em/prior_history.csv`: pool `susie_mix$alpha` from **every gene
  result in `results`**, separately by tissue and coding, to prepare iteration 1.
- Existing history: read its latest iteration, verify all its chunks have
  finished and all gene files exist, then pool `weighted_fit_mix$alpha` from that
  iteration to prepare the next one.
- The empirical Bayes update uses the expected coding assignments of components
  with **`V > 0` whose own 95% CS passes `min_abs_corr = 0.5`**. Workers evaluate
  these CSs with `susie_get_cs(..., dedup = FALSE)` while the genotype matrix is
  available and save the eligible component indices in `fit$coding_prior_cs`.
  Identical and partially overlapping CSs count separately: there is no merging,
  overlap threshold, or selection of one representative. The ordinary reported
  `fit$sets` is retained and is not used to reconstruct EB eligibility.
- With all three codings present in every fit, each prior equals its summed
  eligible `alpha` divided by the number of eligible components. All predictors
  in each selected row contribute, including those outside that row's CS.
  Positive variances below `1e-9` are tested too; exactly zero-variance rows are
  excluded. A fit with no eligible components contributes **zero counts** and
  remains a valid fit. Counts are pooled across the other genes in that tissue.
  Marginal PIP sums remain full-fit diagnostics and do not determine the priors.
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

This is a **purity-filtered empirical Bayes update** for the coding-mixture
SuSiE model. It optimizes the coding objective for the selected component rows;
it is not ordinary EM maximizing the full-data marginal likelihood. Eligibility
can change between iterations, so the full-data ELBO need not increase. The
earlier unfiltered derivation and the distinction from the cTWAS model are
recorded in [EM_PRIOR_REVIEW.md](EM_PRIOR_REVIEW.md).

```text
results_em/
  prior_history.csv
  last_array_job_id.txt
  last_continuation_job_id.txt  # Latest automatically queued preparation, if any
  iteration_001/
    priors.csv                 # Frozen priors used to fit this iteration
    manifest.csv               # Frozen chunk/gene assignments
    source_audit.csv            # Errors and nonconvergence in source results
    component_counts.rds       # Eligible alpha sums by tissue and available coding classes
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
`pi_rec`, `pi_dom`, eligible component `alpha` sums, descriptive PIP sums, fit/error
counts, and a timestamp. `n_components` is the total number of component rows;
`n_active_components` counts all `V > 0` rows, and `n_zero_variance_components`
counts the zero-variance rows. **`n_eligible_components` counts the rows used**;
`alpha_total` equals this count up to rounding.
`n_fits` counts all valid fits, while `n_active_fits` counts fits with at least
one positive-variance component. `n_eligible_fits` counts fits contributing to
the update; `n_zero_cs_fits` counts fits without any eligible component.
`cs_min_abs_corr=0.5` and `cs_dedup=FALSE` identify the selection settings.
`mstep_q_gain` checks that the M-step increases its selected-count
objective; `max_abs_prior_change` tracks movement of the three coding weights.
`source_elbo_sum` records the summed source-fit ELBO over **all valid fits,
including all-null fits**, only when all included
fits supply a finite value (`n_source_elbo` gives the count). Compare these
values only across iterations with the same gene/tissue fits and model settings.
Thus iteration 1's row is estimated from the original scan and is the prior
**used for** iteration 1. Iteration 2's row is estimated from iteration 1.
Each gene RDS contains a list of successful tissues, with `weighted_fit_mix`,
the predictor map, coding priors/weights, and the existing fit metadata. The
full `alpha`, `V`, and posterior moments remain available for warm starts.

## Rerunning with purity filtering

New history rows use `update_method=susie_purity_component_alpha_v3`. Preparation
requires `coding_prior_cs` metadata on every valid source fit. Old fitted objects
without it are rejected before an iteration or history row is written: even
full `alpha` and `V` cannot establish genotype correlations without the original
design matrix, and the reported CS list can have removed duplicate components.

For a fresh RCC rerun:

1. Sync the updated `workhorse.R`, `workhorse_em.R`, and `em_utils.R`, together
   with the current preparation and chunk-runner scripts. Finish or stop the
   previous job sequence before replacing scripts used by running workers.
2. Preserve the old `results` and `results_em` directories separately. Rerun the
   baseline scan into a fresh `results` directory using the updated workhorse.
   This regenerates the initial `susie_mix` source fits and also applies purity
   0.5 to additive, mixed, Slide, and permutation fits. Baseline drivers write
   the same gene filenames, so preserve the old outputs before rerunning.
3. With a fresh `results_em` directory and no old `prior_history.csv`, launch
   `sbatch em_susie_mix` (or supply an iteration count) from `job/`.
4. Regenerate the EM and Slide summaries from these new results.

Setting `min_abs_corr=0.5` only in the EM workhorse does not retroactively filter
the baseline source for iteration 1. A `resume` also keeps its frozen priors;
it is not a restart with newly estimated priors. Historical numeric rows and
frozen snapshots are never rewritten; legacy method tags remain distinguishable
if histories are combined. Keep old and new runs separate for this full rerun.

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
eligible component posteriors retains its previous prior and is marked
`carried_forward_no_eligible_components`. If no previous prior exists, it starts
uniformly and is marked `initialized_uniform_no_eligible_components`; this is a
fallback, not a learned estimate. Zero-CS fits are valid results, not failed
genes. Full fits retain all L rows for subsequent fitting.
All fitting workhorses use `min_abs_corr=0.5`; simulations already used this value.
No association-P, read-count, or lead-PIP filter is added to the existing data QC.
No pseudocount is added to the M-step. Disconnected coding groups retain their previous
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
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_purity.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_marginal_likelihood.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_em_susie_smoke.R
bash script/scan_tissue_attempt/tests/test_em_launcher.sh
```
