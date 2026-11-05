# Reproduction entry points and planned repository cleanup

The new online workflow is self-contained inside the project at
`workflows/slider_online_014_024`. Its `code/` includes the statistical helpers,
fitter, scheduler, Slurm controller/worker scripts and frozen renamed-gene map.
Its own jobs do not depend on a bundle extracted in RCC HOME.

The existing project has several historical launchers. The current source map,
based on the files inspected for this bundle, is:

| Analysis | R sources / driver | Jobs / documentation |
|---|---|---|
| Simulations (SuSiE, mix, discrete slider) | `script/sim/`; `write_jobs.R`, `preflight.R`, `sim_workhorse.R` | `job/launch_simulation_slide`, generated `run_simulation_slide_batch_*`; `script/sim/README_simulations.md`, `RCC_SIMULATION_UPDATE.md` |
| Baseline gene scan (including SuSiE) | `script/scan_tissue_attempt/workhorse.R`, `workhorse_utils.R`, `write_job.R`; generated `script/scan_job/run_chunk_*.R` | `job/test1` is the chunk-array launcher; audit generated drivers before reproducing a historical scan |
| SuSiE-mix EM | `workhorse_em.R`, `em_utils.R`, `prepare_em_iteration.R` under `script/scan_tissue_attempt/` | `job/em_susie_mix`, `job/em_susie_mix_array`; `job/README_em.md` |
| Full discrete-slider EM through 013 | **Production:** `slider_inplace_runtime_2026-10-29/script/scan_tissue_attempt/`; frozen runtime includes the expanded-cohort mapping | `slider_inplace_runtime_2026-10-29/job/resume_slide_prior_em`, `slide_prior_parallel_worker`; runtime `frozen.json` and expansion provenance |
| Online slider 014–019, full slider 020–024 | `workflows/slider_online_014_024/code/engine.R`, `minibatch.py`, vendored `script/scan_tissue_attempt/` | Same `code/` folder: `start.sh`, `status.sh`, `controller_job`, `worker_job`, `handoff.py`; current README |

Historical `job/em_susie_slide_prior` describes the original initialization
workflow. It is not the entry point for the expanded, recovery-managed current
run. Likewise, the continuous-slider `results_slide` scan is a different model
from discrete-slider-prior EM; its script/job chain needs a separate provenance
audit before it is assigned a canonical reproduction command.

The later cleanup should produce one documented workflow directory for each
analysis, with a clear launch/resume/status command, and separate code from
runtime state, logs, temporary PLINK output, fitted results and old bundles.
A suitable target is:

```
workflows/
    simulations/
    susie_scan/
    susie_mix_em/
    susie_slider_em/
    susie_slider_online_em/
shared/                 # verified common helpers only
config/                 # raw-input locations and frozen analysis settings
provenance/             # package/source versions, manifests, seeds, checksums
job/                    # documented entry points or workflow-specific jobs
archive/                # superseded source bundles and launchers
```

Before moving anything, inventory which deployed scripts each active Slurm job
reads, compare generated drivers with their source generator, and preserve
immutable versions and manifests for historical results. Absolute paths in
saved workflow state and warm-start checkpoints must be migrated deliberately;
moving a directory alone can break resumption. Validate one small fixture for
each analysis from the cleaned tree before declaring the old copies redundant.

For exact reproduction, scripts alone are not sufficient. Preserve the R
package source/version (including the discrete-slider implementation), R and
PLINK versions, GTEx/annotation input identities, original and renamed gene
maps, eligibility/exclusion rules, random seeds, model settings, prior
checkpoints and the corresponding source hashes. Protected raw GTEx data need
not be copied into a code archive; document its location, access requirements
and identity. The current online bundle records its inputs and package versions
but is not a replacement for a complete environment/data provenance audit of
all historical analyses.

No broad project reorganization is performed by this bundle. The full cleanup
can be done as a separate change once active paths have been inventoried.
