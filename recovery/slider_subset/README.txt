Resume the stopped expanded slider run with 646 verified additions
================================================================

Diagnosis: iteration 005 completed and validated all 473,380 tissue fits.
Gate 49661983 correctly submitted setup 49662008. Setup stopped at the
mapping preflight: 646 of 680 IDs passed; 34 had differing RCC annotation
start/end coordinates. All 34 had unique matching expression/annotation
IDs and matching names/chromosomes. They were not PLINK or fitting failures.

At the user's request, this patch defers exactly those 34 candidates and
retains the 646 passing IDs unchanged. It does not relax mapping validation.
The cohort is fixed from 006 through 040: 18,263 autosomal expression IDs
(17,617 original + 646 added), or 19,114 manifest entries including the
851 previously excluded non-autosomal entries.

Upload slider_subset_patch_rcc_2026-10-30.zip to /home/wdenault, then:

cd "$HOME" &&
unzip slider_subset_patch_rcc_2026-10-30.zip &&
cd "$HOME/slider_subset_patch_2026-10-30" &&
sha256sum -c SHA256SUMS &&
bash apply_and_resume.sh

The patch targets the existing HOME/slider_inplace_2026-10-29 bundle. It
verifies its original checksums and the actual RCC mapping report, refuses
active slider jobs or an already initialized expanded run, backs up changed
files and the full 680-row report, then installs the fixed selection. A
durable plan permits an interrupted file replacement to be resumed.

It then invokes the existing --resume command. The gate will verify 005,
setup will repeat the strict metadata and genotype checks for all 646 IDs,
and the controller will start 006 and continue toward 040. Completed 005
genes are not refitted by setup: they are reused as warm starts for 006.
New genes start cold at 006; their counts first update the prior for 007.

Results remain in the same results_slide_prior_em folder. The existing
two-checkpoint cleanup policy is unchanged. Prior/count history is retained.
The 34 deferred IDs and reasons are recorded in excluded_candidates.csv
inside the HOME bundle and in results_slide_prior_em/
inplace_expansion_2026-10-29/excluded_candidates.csv after setup.
The full failed RCC report is preserved in
$HOME/slider_subset_patch_state_2026-10-30/mapping_checks.before_subset.csv.

Monitor using:
squeue -u "$USER" -o "%.18i %.28j %.12T %R"
cat "$HOME/slider_inplace_launch_2026-10-29/last_submission.json"

The last receipt identifies the current gate/setup log under
$HOME/slider_inplace_launch_2026-10-29/logs. Once fitting starts, controller
logs are in the project/slider_inplace_runtime_2026-10-29/job directory.

Do not re-extract the original 2026-10-29 ZIP over the patched bundle.
If new input checks fail, inspect the reported error; additional genes are
never silently discarded. No existing EM result is deleted by this patch.
