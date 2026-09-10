# One EM iteration per submission

On the cluster, from the project's `job` directory:

```bash
sbatch em_susie_mix
```

The preparation job loads `R/4.2.0`, estimates the tissue priors, saves them,
and submits one Slurm array using the resource settings in `test1`
(`broadwl`, 23 hours, 40 GB, one CPU per chunk). It discovers the existing
`data/temp_index/chunk_*_genes.txt` lists: currently 185 chunks / 18,468 genes.
The original generated `run_chunk_*.R` scripts are not changed or executed.

When the array finishes, run the same command again for the next iteration.
Each submission performs exactly one iteration; it does not automatically
continue until convergence. The preparation job finishing only means the
array was submitted. Its output gives the array job ID.

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
  iteration_001/
    priors.csv                 # Frozen priors used to fit this iteration
    manifest.csv               # Frozen chunk/gene assignments
    source_audit.csv            # Errors and nonconvergence in source results
    array_job_ids.txt
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

As in the original scan, gene/tissue errors are saved as explicit records;
they do not stop other genes. Completion means all genes were attempted,
not that all fits succeeded. Error records contribute no PIPs and appear in
the next iteration's `source_audit.csv` and history error counts. Review each
chunk's CSV and logs before advancing. Nonconverged fits contribute their
PIPs and are flagged in the audit. An entirely failed source scan is rejected.

Missing/unreadable gene files, mismatched predictor maps, or invalid PIPs stop
preparation. The initial results must match the full chunk gene list, so a
partially finished original scan cannot initialize the priors. A tissue with
zero total PIP has an undefined initial prior and stops preparation. In later
iterations its previous prior is retained and explicitly marked
`carried_forward_zero_pip`. There is no pseudocount or probability floor.

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
