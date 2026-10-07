# Tissue-specific slider-prior EM

For the interrupted October initialization and subsequent bounded scheduling,
use [the recovery workflow](README_slide_prior_recovery.md). It preserves
saved fits, audits errors hidden behind old completion markers, and processes
10 genes per task in isolated R processes under a 300-job submission budget.
The historical launcher described below refuses recovery-managed runs.

This workflow learns 17 slider probabilities separately for each tissue using
`susieRSlidePrior::susie`. The existing `results_slide` scan used the continuous
slider and cannot supply the discrete joint posteriors needed here. A new run
therefore fits uniform probabilities first (`iteration_000`), then learns and
refits at iterations 001, 002, etc. It writes to `results_slide_prior_em`.

Upload/extract `output/slide_prior_em_rcc.zip` into the RCC project, preserving
its `job/` and `script/` paths. The archive contains code, documentation and
tests only. It does not contain fitted results or an R package installation.
The RCC R/4.2.0 environment must have `susieRSlidePrior >= 0.3.0` (the discrete
prior implementation with `delta_prior_counts` and `model_init` support),
`data.table`, and `matrixStats` installed, together with their dependencies.
This was tested locally with susieRSlidePrior 0.3.0 from the `susie_slide_prior`
branch at commit `8f99d4a`.

From the project root, a short package/EM check is:

```bash
module load R/4.2.0
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_em_smoke.R
```

Start initialization plus five EM updates:

```bash
cd /project2/mstephens/wdenault/susie_mix/job
sbatch em_susie_slide_prior 5
```

On a fresh run this requests **six complete fitting passes**: 000 with uniform
weights, followed by 001–005 with learned weights. Use 8 for initialization
plus eight updates. Once iteration 005 is complete, submitting the command
with 5 again continues with 006–010. No mix fits or continuous-slider fits are
used as warm starts. Each later slider iteration uses its own previous full
discrete fits through `model_init`.

The initial gene manifest and tissue names come from the latest local-to-RCC
`results_em/prior_history.csv` and corresponding mix `manifest.csv`. These
files are read only; no mix weights enter the slider analysis. Each new pass
uses up to 298 balanced gene chunks (61–62 genes for the 18,468-gene scan).
The launcher queues the next prior update only after the entire current array
succeeds. This is a fixed number of updates, without an automatic outer-EM
convergence stop.

The grid is `seq(-1, 1, length.out=17)`: -1 is recessive, 0 additive and +1
dominant. The columns `w_01` through `w_17` follow that order; `grid.csv` makes
the mapping explicit. For tissue t and grid point k the M-step is:

```text
N[t,k] = sum over genes g, components l with V[g,l] > 0,
         and free SNPs j of alpha_delta[g,l,j,k]
w_new[t,k] = N[t,k] / sum_k N[t,k]
```

`delta_prior_counts` already excludes SNPs forced to additive coding and
explicit null columns. Its rows retain the actual posterior mass on free
SNPs: a row summing to 0.2 contributes 0.2, not 1. Exactly zero-variance
components are integrated out and contribute zero. If a tissue has no free
mass, its previous prior is retained. No pseudocounts or probability floor are
added. The update uses component probabilities, not PIPs or counts of reported
credible sets. A completely uninformative free assignment reproduces the
current prior and supplies no directional evidence for changing it.

There is **no purity, association P-value, read-depth or N >= 300 screen in
prior learning**. Existing input QC remains: at least 50 matched individuals
to fit a tissue, genotype MAF/HWE/variation checks, and `min_obs=5` forcing
delta=0 where any genotype class has fewer than five individuals. Fits use
L=10, `standardize=FALSE`, Gaussian variance updates by EM, `tol=1e-5` and
`max_iter=1000`. CS purity 0.5 controls reporting only. The stored `min_pv` is
an additive marginal association P-value for later summaries; it never enters
the prior update. Read-count and sample-size metadata are also saved.

Iteration 000 establishes the successful gene/tissue cohort. Initial gene
and tissue failures remain explicit audited records. Subsequent iterations
preserve that cohort; losing a previously successful fit stops the chunk.
Nonconvergence also blocks completion and the next prior update. Initial
failure records can be inspected in `completed/chunk_*.csv`, the gene RDS
files and, after pooling, `source_audit.csv`.

After stopped or failed jobs are no longer active, retry unfinished chunks:

```bash
sbatch em_susie_slide_prior resume 5
```

For iteration 000, this means finish initialization and then run five
updates. For a later iteration, the resumed update counts as the first of the
five. Completed chunks and successful saved genes are reused, and frozen
priors are not recomputed on resume. A queued automatic continuation must
finish or be cancelled before a manual launch; the launcher prints its ID.

If a fit reached its numerical iteration limit, inspect the log and increase
only that budget, for example:

```bash
SUSIE_SLIDE_MAX_ITER=5000 sbatch em_susie_slide_prior resume 5
```

This does not change the model, tolerance, prior or cohort. The larger budget
is inherited by the chained jobs and recorded on newly fitted gene results.
Other fitting settings and the package version are frozen; a change requires
inspection rather than silently mixing runs.

Each iteration saves `priors.csv`, `settings.rds`, `grid.csv`, `manifest.csv`,
`pooled_counts.csv`, `source_audit.csv`, full gene `results/*.rds`, completion
markers and logs. Each gene result is a named list of tissues, with the full
fit in `fit_slide_prior`. `prior_history.csv` includes changes in weights,
free posterior mass, active-component counts, M-step Q gains, and source-fit
ELBO sums. `source_elbo` on row m describes the fits from iteration m-1; these
are variational lower bounds, not exact marginal log likelihoods. Grid
probabilities are fixed inside each call; `estimate_prior_method="EM"` in
the fitter updates Gaussian effect variances, whereas this outer workflow
updates the slider probabilities.

Local verification covers fractional counts and forced/null/zero-V cases,
uninformative assignments, immutable snapshots, failed-fit handling and
restart behavior, mocked Slurm dependencies, five real-package prior updates,
and a synthetic GTEx worker run including data import and warm starts. Full
GTEx execution and cluster resource usage still require the RCC run.

```bash
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_em.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_em_smoke.R
Rscript --vanilla script/scan_tissue_attempt/tests/test_slide_prior_em_worker.R
bash script/scan_tissue_attempt/tests/test_slide_prior_em_launcher.sh
```
