# Slider EM handoff: online 014–019, full EM 020–024

This revision replaces the separate-folder minibatch bundle. It continues in
`/project2/mstephens/wdenault/susie_mix/results_slide_prior_em` and uses completed
iteration 013 as its fixed starting checkpoint. It requires iteration 014 and
all later iteration folders to be absent before handoff. It never overwrites
an existing iteration or silently changes the requested starting number.

## Launch on RCC

Upload `slider_online_014_rcc_2026-11-04.zip` to `/home/wdenault`, then:

```bash
cd /home/wdenault &&
unzip slider_online_014_rcc_2026-11-04.zip &&
cd /home/wdenault/slider_online_014_2026-11-04 &&
sha256sum -c SHA256SUMS &&
bash start.sh --handoff
```

The new code, jobs, plan, provenance, state and controller logs are installed in:

```
susie_mix/workflows/slider_online_014_024/
    code/          # immutable, checksummed R/Python sources and Slurm jobs
    config.json
    HANDOFF.json
    seed/          # checkpoint certificate, counts, cohorts and original history
    plan/
    stages.csv
    logs/
    tasks/         # submitted arrays and per-worker logs
```

Once installed, the workflow is entirely inside `susie_mix`; no continuing job
depends on the extracted HOME bundle. Raw GTEx inputs and installed R packages
remain external dependencies, as they are for the existing analysis. Their
paths/settings and R package versions are captured and checked at startup.

Check progress:

```bash
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_014_024/code/status.sh
squeue --user=wdenault -o "%.18i %.28j %.12T %.12l %R"
```

After inspecting an interruption, resume using:

```bash
bash /project2/mstephens/wdenault/susie_mix/workflows/slider_online_014_024/code/start.sh --resume --handoff
```

`--handoff` is idempotent after `HANDOFF.json` is committed. If the initial
launch could not obtain a safe pending-controller window, the workflow files
remain installed; retry with `--resume --handoff`. Do not use the previous
separate-folder bundle or the old `resume_slide_prior_em latest` to continue
this new history format. The launcher refuses an earlier separate minibatch
experiment still in the queue and refuses duplicate new controllers.

## What happens to the current jobs

Two different controllers cannot own iteration 014 in the same directory.
The handoff therefore waits until the existing controller is pending for
`BeginTime`, at least 90 seconds before its earliest start. It verifies the
job ID against the saved continuation ID and checks the command, owner and
state with `scontrol`. It then cancels **only that pending controller**, using
an exact job ID plus user/name/state filters. The state filter is documented
in [Slurm's scancel reference](https://slurm.schedmd.com/scancel.html).

All running and pending gene workers remain untouched. The same existing,
unmodified controller script is resubmitted with target iteration 013 and a
23-hour limit, so it can refill remaining genes and certify 013 without
preparing 014. Its descendants inherit target 013. The new online controller
waits for every old slider job to leave the queue and for 013's completion
certificate before writing new results or extending root prior history.

The pending-controller replacement and submission IDs are journaled. An
ambiguous sbatch response leaves `pending_submission.json` and stops automatic
retry to avoid duplicate jobs. If the old finisher fails before certifying 013,
the online branch stops and reports that failure; inspect `logs/finish_013_*`
and the subsequent original runtime controller logs before retrying it.

If 014 already exists, this bundle stops before cancelling any job. Keep those
results and select a later handoff boundary in a subsequent revision; do not
delete a started iteration just to satisfy this launcher's check.

## Numbering and result layout

| Scans | Method | Updates per scan |
|---|---|---|
| 000–013 | Existing full EM | Existing history preserved |
| 014–019 | Online/minibatch EM | Nine updates, rho = 0.3 |
| 020–024 | Full-data EM | One update |

For example:

```
results_slide_prior_em/
    prior_history.csv
    prior_history_long.csv
    online_latest_priors.csv
    iteration_013/
    iteration_014/
        ONLINE_OWNER.json
        priors.csv                 # prior at the start of the scan
        prior_history.csv          # accumulated history through this scan
        prior_history_long.csv     # this scan's grid values/probabilities
        settings.rds
        grid.csv
        manifest.csv
        results/GENE.rds           # one output per gene, as before
        batches/batch_01/
            priors.csv             # prior actually used for this batch
            prior_after.csv        # prior after its update
            counts.csv
            running_counts_per_gene.rds
            history.csv
            summaries/
            success/
            receipts/
        ... batch_09/
    ... iteration_019/
    iteration_020/batches/full_step/
    ... iteration_024/
```

Every saved gene has `em_iteration` equal to the visible scan number (014,
015, etc.), plus `online_stage` and `online_batch` attributes. Each online
scan fits every gene once, under its own batch's frozen prior. Therefore the
genes within an online scan need not all have the same fitted prior; their
batch metadata records the exact one. Every full-data scan uses one prior
for all genes.

The final results are fitted under `iteration_024/batches/full_step/priors.csv`.
`prior_after.csv` and root `online_latest_priors.csv` are the resulting final
M-step, one update ahead of those fits. No extra unrequested full scan is run.

## Prior history at every within-scan update

The existing columns and their order, including `w_01`–`w_17`, remain at the
front of root `prior_history.csv`. Old rows retain their original values.
A byte-for-byte copy of the original history is kept at
`workflows/slider_online_014_024/original_prior_history.csv`, independently of
the preserved source history in `seed/`.

New identifiers are:

- `iteration`: the visible scan number.
- `batch`: 0 for the scan's input prior, 1–9 for online updates, or 1 for the
  update after a full-data scan.
- `prior_role`: `input` or `after_update`.
- `phase`: `minibatch` or `full`.
- `update` / `stage`: global online-workflow update 1–59; an input row refers
  to the preceding committed stage. Legacy rows have no online stage.
- `rho`: 0.3 online, 1 for full-data updates.
- `slider_01`–`slider_17`: the actual grid values paired with `w_01`–`w_17`.
- `checkpoint_iteration`: 013 for newly generated rows.

The grid is `seq(-1,1,length.out=17)`: -1 recessive, 0 additive, +1 dominant.
Weights are probabilities at these grid values. Each update records every
tissue's complete 17-value distribution, not just a mean or the last batch.

`prior_history_long.csv` supplies one row per tissue/update/grid point:

```
iteration,batch,phase,prior_role,update,tissue,rho,slider_index,slider_value,probability
```

For example, select all post-update priors within scan 014 in R:

```r
root <- "/project2/mstephens/wdenault/susie_mix/results_slide_prior_em"
h <- read.csv(file.path(root, "prior_history_long.csv"))
h14 <- subset(h, iteration == 14 & prior_role == "after_update")
sort(unique(h14$batch))  # 1 through 9 after the scan completes
```

Histories are updated only after an entire batch is validated, never from a
partial batch or completion-order subset. An atomic publication journal allows
restart after a CSV write and refuses unrelated changes by another writer.
The scan-input row is published after that scan has been prepared.

## Statistics, scheduling and storage

The 18,263-ID expanded autosomal cohort, including the 646 recovered IDs, is
fixed from certified 013; the 34 deferred mappings stay excluded. The fitter,
convergence tolerance, QC, active-component and free-SNP count rules are
unchanged. Source counts C initialize S=C/18,263 and the first 014 prior is
normalize(C), the ordinary next M-step from 013.

Each of six passes reshuffles all genes into two batches of 2,030 and seven
of 2,029 (fixed seed 20261104 plus pass, saved assignments). Batches are balanced
by gene count, without runtime/credible-set selection. Each batch updates:

```
S_new[t,k] = 0.7*S_old[t,k] + 0.3*batch_counts[t,k]/batch_gene_count
prior_new[t,k] = S_new[t,k]/sum_k S_new[t,k]
```

Zero-mass tissue rows retain their preceding prior. Five ordinary full-data
E/M steps follow, using the established pooled-count M-step and Q check.
The schedule has 59 updates and stops at 024; convergence must still be
assessed from the prior history. Counts are not normalized per gene before
pooling and no pseudocounts or credible-set filters are introduced.

The default cap is 298 workers (48 GB, 23 hours, one gene per array element),
refilled approximately every five minutes. Controllers use 16 GB and 23 hours.
All user jobs share the 300-slot budget; space is reserved for continuation
controllers. A successor is queued before lengthy metadata operations. Batch
completion, frozen priors, hashes and cohort checks prevent partial updates.
Saved validated outputs can be adopted after an interruption without refitting;
the commit transaction prevents applying rho twice.

Existing 013 fits are pinned with hard links while they are needed as warm
starts. During online 014, only 013 and 014 retain full results; old 000–012
result files are retired after the seed and first batch preparation validate.
Before 015, all 014 genes have validated branch warm starts, so the seed pins
and 013 full results can be removed. This continues until only 023 and 024
full results remain. Every prior, count summary, history, audit and provenance
file remains. Unexpected objects in a results directory stop cleanup for
inspection; no arbitrary folder deletion is performed.

The existing compact-result fix is retained: the captured R `$call` is removed,
while fitted posteriors remain available for warm starts. There is no serial
bulk recompression pass.

## Reproducibility and later repository cleanup

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the current entry points and a
proposed organization for simulations, baseline scans, mix EM, slider EM and
online EM. This bundle does not relocate active scripts or delete old source
bundles. It does provide one immutable local source/job snapshot for the new
workflow, frozen plans, original priors, package versions, exact gene mapping,
configuration, checkpoint identity and every subsequent prior update.

Local checks cover all 59 synthetic stages, history preservation and duplicate
prevention, all nine within-scan updates, displayed iteration numbers, shared
result folders, real slider warm starts, ownership-aware cleanup, queue
reservation, and the controller-only handoff. No production RCC job was run
during these tests. Development tests use local fixture/library paths and are
not required on the RCC login node.
