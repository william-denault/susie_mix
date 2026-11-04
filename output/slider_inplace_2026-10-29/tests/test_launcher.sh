#!/bin/bash
# Mock Slurm and R: no cluster access, fits, or real result changes.
set -euo pipefail
repo=$(pwd -P)
mkdir -p "$repo/tmp"
fixture=$(mktemp -d "$repo/tmp/parallel_launcher.XXXXXX")
export PAR_TEST="$fixture" SUSIE_MIX_PROJECT_DIR="$fixture/project with spaces" SLURM_JOB_ID=9000 USER=test_user
mkdir -p "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em"
controller="$SLIDER_TEST_RUNTIME/job/resume_slide_prior_em"
worker="$SLIDER_TEST_RUNTIME/job/slide_prior_parallel_worker"
export PAR_WORKERS=0 PAR_OTHER=0 PAR_SIM=0 PAR_CONTROL=0 PAR_FAIL_QUEUE=0 PAR_SELF=1 PAR_STATUS=batch
export PAR_FAIL_R=0 PAR_FAIL_ARRAY=0 PAR_FAIL_NEXT=0 PAR_AMBIGUOUS=0 PAR_FAIL_GENE=0 PAR_LOCK=0
export PAR_ACTION=fit
export PAR_UNIT=gene
unset SUSIE_SLIDE_WORKERS SUSIE_SLIDE_POLL_MINUTES SUSIE_SLIDE_TARGET_ITERATION SUSIE_SLIDE_WAIT_FOR_SIMULATIONS
module() { :; }
flock() { printf '%s\n' "$*" >> "$PAR_TEST/lock_calls"; [[ "$PAR_LOCK" == 0 ]]; }
squeue() {
    [[ "$PAR_FAIL_QUEUE" == 0 ]] || return 1
    [[ "$*" == *--array* && "$*" == *--states=all* && "$*" == *%200j* ]] || return 1
    [[ "$PAR_SELF" == 0 ]] || printf '9000|9000|slide_inplace_control|RUNNING\n'
    local i
    for ((i=1;i<=PAR_WORKERS;i++)); do printf '9100_%s|9100|slide_inplace_fit|PENDING\n' "$i"; done
    for ((i=1;i<=PAR_OTHER;i++)); do printf '8000_%s|8000|unrelated|RUNNING\n' "$i"; done
    for ((i=1;i<=PAR_SIM;i++)); do printf '7000_%s|7000|susie_slide_sim    |PENDING\n' "$i"; done
    [[ "$PAR_CONTROL" == 0 ]] || printf '9200|9200|slide_recover_control|RUNNING\n'
    printf '9999|9999|slide_inplace_control|COMPLETED\n'
    return 0
}
Rscript() {
    [[ "$2" == "$SUSIE_MIX_PROJECT_DIR/slider_inplace_runtime_2026-10-29/script/scan_tissue_attempt/run_slide_prior_parallel.R" ]] || return 1
    local mode=$3
    printf '%s\n' "$mode" >> "$PAR_TEST/r_calls"
    [[ "$PAR_FAIL_R" == 0 ]] || return 1
    case "$mode" in
        plan)
            printf '%s\n' "$6" > "$PAR_TEST/capacity"
            [[ -f "$7" ]] || return 1
            if [[ "$PAR_STATUS" != batch || "$6" == 0 ]]; then
                printf '%s\n0\n' "${PAR_STATUS/batch/waiting}"; return
            fi
            local b="$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/iteration_000/recovery/parallel/batches/test"
            mkdir -p "$b"
            printf 'batch\n0\n%s\n%s\n%s\n' "$b" "$6" "$PAR_ACTION" ;;
        task)
            if [[ "$PAR_UNIT" == gene ]]; then printf 'G001\tvalidate\nG002\tvalidate\n'
            else printf 'G001\tfit\t%s\n' "$PAR_UNIT"; fi ;;
        gene)
            printf '%s\n' "$8" >> "$PAR_TEST/gene_calls"
            printf '%s\n' "${9:-gene}" >> "$PAR_TEST/unit_calls"
            [[ "$PAR_FAIL_GENE" == 0 || "$8" != G001 ]] ;;
        *) return 1 ;;
    esac
}
sbatch() {
    printf '%s\n' "$*" >> "$PAR_TEST/submissions"
    if [[ "$*" == *--array=* ]]; then
        printf '%s\n' "$@" > "$PAR_TEST/array_args"
        [[ "$PAR_FAIL_ARRAY" == 0 ]] || return 1
        [[ "$PAR_AMBIGUOUS" == 0 ]] || { printf 'unexpected response\n'; return; }
        printf '9100;cluster\n'
    else
        printf '%s\n' "$@" > "$PAR_TEST/next_args"
        printf '%s\n' "$SUSIE_SLIDE_TARGET_ITERATION" "$SUSIE_SLIDE_POLL_MINUTES" "$SUSIE_SLIDE_WORKERS" > "$PAR_TEST/exported"
        [[ "$PAR_FAIL_NEXT" == 0 ]] || return 1
        printf '9101;cluster\n'
    fi
}
python3() { [[ "$1" == *slider_inplace_runtime_2026-10-29/verify_inplace.py ]]; }
export -f module flock squeue Rscript sbatch python3
reset_calls() { rm -f -- "$PAR_TEST/submissions" "$PAR_TEST/r_calls" "$PAR_TEST/array_args" "$PAR_TEST/next_args"; }
run() { "$BASH" "$controller" latest >"$PAR_TEST/out" 2>"$PAR_TEST/err"; }
reject() { if run; then echo 'Expected rejection' >&2; exit 1; fi; }
run
[[ $(cat "$PAR_TEST/capacity") == 298 ]]
grep -q -- '--array=1-298' "$PAR_TEST/array_args"
grep -q -- '--begin=now+5minutes' "$PAR_TEST/next_args"
grep -q -- '--dependency=afterany:9000' "$PAR_TEST/next_args"
[[ $(cat "$PAR_TEST/exported") == $'40\n5\n298' ]]
[[ ! -f "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/parallel_pending_submission.txt" ]]
# A single completion is replenished: no threshold at 200 workers.
reset_calls; export PAR_WORKERS=297; run
[[ $(cat "$PAR_TEST/capacity") == 1 ]]
grep -q -- '--array=1-1' "$PAR_TEST/array_args"
# Other jobs count, including pending elements; never exceed the shared budget.
reset_calls; export PAR_WORKERS=100 PAR_OTHER=150; run
[[ $(cat "$PAR_TEST/capacity") == 48 ]]
reset_calls; export PAR_OTHER=198; run
[[ $(cat "$PAR_TEST/capacity") == 0 && ! -f "$PAR_TEST/array_args" && -f "$PAR_TEST/next_args" ]]
# Existing simulation chain keeps priority across its 300-task batch boundaries.
reset_calls; export PAR_WORKERS=0 PAR_OTHER=0 PAR_SIM=298; run
[[ ! -f "$PAR_TEST/r_calls" && ! -f "$PAR_TEST/array_args" && -f "$PAR_TEST/next_args" ]]
grep -q 'Waiting for the existing simulation chain' "$PAR_TEST/out"
export PAR_SIM=0
reset_calls; export PAR_STATUS=complete; run
[[ ! -f "$PAR_TEST/submissions" ]]
export PAR_STATUS=batch
for flag in PAR_FAIL_QUEUE PAR_CONTROL PAR_FAIL_R PAR_LOCK; do
    reset_calls; export "$flag=1"; reject
    [[ ! -f "$PAR_TEST/submissions" ]]; export "$flag=0"
done
reset_calls; export PAR_SELF=0; reject; export PAR_SELF=1
[[ ! -f "$PAR_TEST/submissions" ]]
reset_calls; export PAR_FAIL_ARRAY=1; reject; export PAR_FAIL_ARRAY=0
[[ ! -f "$PAR_TEST/next_args" && ! -f "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/parallel_pending_submission.txt" ]]
reset_calls; export PAR_FAIL_NEXT=1; reject; export PAR_FAIL_NEXT=0
[[ -f "$PAR_TEST/array_args" && ! -f "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/parallel_pending_submission.txt" ]]
reset_calls; export PAR_AMBIGUOUS=1; reject; export PAR_AMBIGUOUS=0
[[ -f "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/parallel_pending_submission.txt" ]]
reset_calls; reject; [[ ! -f "$PAR_TEST/submissions" ]]
rm -f -- "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/parallel_pending_submission.txt"
# A failed validation gene does not discard successful neighbors in its task.
export SLURM_ARRAY_TASK_ID=1 PAR_FAIL_GENE=1
batch="$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/iteration_000/recovery/parallel/batches/test"
if "$BASH" "$worker" "$SUSIE_MIX_PROJECT_DIR" 0 "$batch" >"$PAR_TEST/worker.out" 2>"$PAR_TEST/worker.err"; then exit 1; fi
[[ $(cat "$PAR_TEST/gene_calls") == $'G001\nG002' ]]
export PAR_FAIL_GENE=0 PAR_STATUS=batch PAR_OTHER=0 PAR_SIM=0 PAR_WORKERS=0
unset SLURM_ARRAY_TASK_ID
reset_calls; export PAR_OTHER=299; reject
[[ ! -f "$PAR_TEST/submissions" ]]
export PAR_OTHER=0
# The old shared-scheduler flag must not dispatch this private controller.
mkdir -p "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/shared_scheduler"
printf '{}\n' > "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/shared_scheduler/enabled.json"
export SUSIE_SLIDE_WAIT_FOR_SIMULATIONS=0 PAR_SIM=150
reset_calls; run
[[ $(cat "$PAR_TEST/capacity") == 148 ]]
grep -Fq -- "$SUSIE_MIX_PROJECT_DIR/slider_inplace_runtime_2026-10-29/job/slide_prior_parallel_worker" "$PAR_TEST/array_args"
grep -Fq -- "$SUSIE_MIX_PROJECT_DIR/slider_inplace_runtime_2026-10-29/job/resume_slide_prior_em" "$PAR_TEST/next_args"
echo 'PASS: private routing, target 040, 298 cap, 5-minute refill, simulations share 300 slots, duplicate/failure guards, and validation reuse.'
