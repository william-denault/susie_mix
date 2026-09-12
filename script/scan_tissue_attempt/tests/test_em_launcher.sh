#!/bin/bash
# Run from the repository root with bash. All cluster commands are mocked;
# this test never submits real jobs or reads actual scan results.
set -euo pipefail
repo_dir=$(pwd -P)
mkdir -p "$repo_dir/tmp"
test_base=$(cd "$repo_dir/tmp" && pwd -P)
test_dir=$(mktemp -d "$test_base/em_launcher_test.XXXXXX")
cleanup() {
    local resolved
    resolved=$(cd "$test_dir" && pwd -P)
    if [[ "$resolved" == "$test_base"/em_launcher_test.* ]]; then
        rm -rf -- "$resolved"
    fi
}
trap cleanup EXIT
export SUSIE_MIX_PROJECT_DIR="${test_dir}/project with spaces"
export EM_TEST_ARGS="${test_dir}/submission.txt"
export EM_TEST_CONT_ARGS="${test_dir}/continuation.txt"
export EM_TEST_CALLS="${test_dir}/calls.txt"
export EM_TEST_R_ARGS="${test_dir}/r_args.txt"
export EM_TEST_FAIL_R=0 EM_TEST_ACTIVE=0 EM_TEST_LOCKED=0 EM_TEST_FAIL_SUBMIT=0
export EM_TEST_FAIL_CONTINUATION=0 EM_TEST_ITERATION=1 EM_TEST_ACTIVE_CONTINUATION=0
export SLURM_JOB_ID=4000
export USER="${USER:-em_test}"

module() { :; }
flock() { [[ "$EM_TEST_LOCKED" == 0 ]]; }
squeue() {
    if [[ "$EM_TEST_ACTIVE" == 1 ]]; then printf '4242\n'; fi
    if [[ "$EM_TEST_ACTIVE_CONTINUATION" != 0 ]]; then printf '%s\n' "$EM_TEST_ACTIVE_CONTINUATION"; fi
}
Rscript() {
    printf '%s\n' "$@" > "$EM_TEST_R_ARGS"
    if [[ "$EM_TEST_FAIL_R" == 1 ]]; then return 1; fi
    if [[ "$*" == *prepare_em_iteration.R* ]]; then
        local iteration_dir
        printf -v iteration_dir '%s/results_em/iteration_%03d' "$SUSIE_MIX_PROJECT_DIR" "$EM_TEST_ITERATION"
        mkdir -p "$iteration_dir"
        printf '%s\n' "$iteration_dir"
        local i sep=""
        for ((i=1; i<=185; i++)); do printf '%s%d' "$sep" "$i"; sep=","; done
        printf '\n'
    fi
}
sbatch() {
    local arg is_array=0
    for arg in "$@"; do if [[ "$arg" == --array=* ]]; then is_array=1; fi; done
    if [[ "$is_array" == 1 ]]; then
        printf '%s\n' "$@" > "$EM_TEST_ARGS"
        printf 'array\n' >> "$EM_TEST_CALLS"
        if [[ "$EM_TEST_FAIL_SUBMIT" == 1 ]]; then return 1; fi
        printf '%d;test_cluster\n' "$((4242 + 2 * (EM_TEST_ITERATION - 1)))"
    else
        printf '%s\n' "$@" > "$EM_TEST_CONT_ARGS"
        printf 'continuation\n' >> "$EM_TEST_CALLS"
        if [[ "$EM_TEST_FAIL_CONTINUATION" == 1 ]]; then return 1; fi
        printf '%d;test_cluster\n' "$((4243 + 2 * (EM_TEST_ITERATION - 1)))"
    fi
}
export -f module flock squeue Rscript sbatch

"$BASH" "$repo_dir/job/em_susie_mix" > "${test_dir}/out" 2>&1
grep -q -- '--array=1,2,3,.*184,185' "$EM_TEST_ARGS"
grep -Fxq "${SUSIE_MIX_PROJECT_DIR}/job/em_susie_mix_array" "$EM_TEST_ARGS"
grep -Fxq "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001" "$EM_TEST_ARGS"
grep -Fxq '4242' "${SUSIE_MIX_PROJECT_DIR}/results_em/last_array_job_id.txt"
grep -Fxq 'new' "$EM_TEST_R_ARGS"
[[ ! -e "$EM_TEST_CONT_ARGS" ]]

# Reject invalid counts before preparation or submission.
for count in 0 -1 1.5 abc 02 999999999999999999999; do
    rm -f -- "$EM_TEST_ARGS" "$EM_TEST_R_ARGS"
    if "$BASH" "$repo_dir/job/em_susie_mix" "$count" > "${test_dir}/out" 2>&1; then
        echo "Expected invalid count $count to fail" >&2; exit 1
    fi
    [[ ! -e "$EM_TEST_ARGS" && ! -e "$EM_TEST_R_ARGS" ]]
done
if "$BASH" "$repo_dir/job/em_susie_mix" 2 extra > "${test_dir}/out" 2>&1; then
    echo 'Expected excess arguments to fail' >&2; exit 1
fi
if env -u SLURM_JOB_ID "$BASH" "$repo_dir/job/em_susie_mix" 2 > "${test_dir}/out" 2>&1; then
    echo 'Expected multiple iterations outside sbatch to fail' >&2; exit 1
fi

# An active array, failed prior preparation, or held lock must not submit.
for failure in EM_TEST_ACTIVE EM_TEST_FAIL_R EM_TEST_LOCKED; do
    rm -f -- "$EM_TEST_ARGS"
    export "$failure=1"
    if "$BASH" "$repo_dir/job/em_susie_mix" > "${test_dir}/out" 2>&1; then
        echo "Expected launch to fail with $failure" >&2; exit 1
    fi
    [[ ! -e "$EM_TEST_ARGS" ]]
    export "$failure=0"
done

"$BASH" "$repo_dir/job/em_susie_mix" resume > "${test_dir}/out" 2>&1
grep -Fxq 'resume' "$EM_TEST_R_ARGS"
export EM_TEST_FAIL_SUBMIT=1
if "$BASH" "$repo_dir/job/em_susie_mix" resume > "${test_dir}/out" 2>&1; then
    echo 'Expected failed array submission to stop the launcher' >&2; exit 1
fi
grep -q 'Retry: sbatch em_susie_mix resume' "${test_dir}/out"
[[ ! -e "$EM_TEST_CONT_ARGS" ]]
export EM_TEST_FAIL_SUBMIT=0

# Three iterations submit three arrays and exactly two preparations. A next
# preparation waits for the whole array and the current lock-holding job.
rm -f -- "$EM_TEST_CALLS"
"$BASH" "$repo_dir/job/em_susie_mix" 3 > "${test_dir}/out" 2>&1
grep -Fxq -- '--dependency=afterok:4242:4000' "$EM_TEST_CONT_ARGS"
grep -Fxq -- '--kill-on-invalid-dep=yes' "$EM_TEST_CONT_ARGS"
grep -Fxq -- '--export=ALL' "$EM_TEST_CONT_ARGS"
grep -Fxq -- "--chdir=${SUSIE_MIX_PROJECT_DIR}/job" "$EM_TEST_CONT_ARGS"
grep -Fxq "${SUSIE_MIX_PROJECT_DIR}/job/em_susie_mix" "$EM_TEST_CONT_ARGS"
[[ "$(tail -n 1 "$EM_TEST_CONT_ARGS")" == 2 ]]
grep -Fxq '4243' "${SUSIE_MIX_PROJECT_DIR}/results_em/last_continuation_job_id.txt"
grep -Fxq '4243' "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001/continuation_job_ids.txt"

# A second manual submission cannot jump ahead of a queued continuation.
export EM_TEST_ACTIVE_CONTINUATION=4243
rm -f -- "$EM_TEST_ARGS" "$EM_TEST_CONT_ARGS" "$EM_TEST_R_ARGS"
if "$BASH" "$repo_dir/job/em_susie_mix" > "${test_dir}/out" 2>&1; then
    echo 'Expected queued continuation to block manual submission' >&2; exit 1
fi
[[ ! -e "$EM_TEST_ARGS" && ! -e "$EM_TEST_CONT_ARGS" && ! -e "$EM_TEST_R_ARGS" ]]

# The continuation recognizes its own Slurm ID and starts the next iteration.
export SLURM_JOB_ID=4243 EM_TEST_ITERATION=2
"$BASH" "$repo_dir/job/em_susie_mix" 2 > "${test_dir}/out" 2>&1
grep -Fxq -- '--dependency=afterok:4244:4243' "$EM_TEST_CONT_ARGS"
[[ "$(tail -n 1 "$EM_TEST_CONT_ARGS")" == 1 ]]
grep -Fxq 'new' "$EM_TEST_R_ARGS"
export SLURM_JOB_ID=4245 EM_TEST_ITERATION=3 EM_TEST_ACTIVE_CONTINUATION=4245
rm -f -- "$EM_TEST_CONT_ARGS"
"$BASH" "$repo_dir/job/em_susie_mix" 1 > "${test_dir}/out" 2>&1
[[ ! -e "$EM_TEST_CONT_ARGS" ]]
[[ "$(grep -c '^array$' "$EM_TEST_CALLS")" == 3 ]]
[[ "$(grep -c '^continuation$' "$EM_TEST_CALLS")" == 2 ]]
grep -Fxq '4246' "${SUSIE_MIX_PROJECT_DIR}/results_em/last_array_job_id.txt"

# A resume count includes the resumed iteration; its successor is always new.
export EM_TEST_ACTIVE_CONTINUATION=0
"$BASH" "$repo_dir/job/em_susie_mix" resume 2 > "${test_dir}/out" 2>&1
grep -Fxq 'resume' "$EM_TEST_R_ARGS"
[[ "$(tail -n 1 "$EM_TEST_CONT_ARGS")" == 1 ]]
if grep -Fxq 'resume' "$EM_TEST_CONT_ARGS"; then
    echo 'Continuation must prepare a new iteration' >&2; exit 1
fi

# An array submission failure never queues a continuation.
export EM_TEST_FAIL_SUBMIT=1
rm -f -- "$EM_TEST_CONT_ARGS"
if "$BASH" "$repo_dir/job/em_susie_mix" 2 > "${test_dir}/out" 2>&1; then
    echo 'Expected array submission failure to stop the chain' >&2; exit 1
fi
[[ ! -e "$EM_TEST_CONT_ARGS" ]]
export EM_TEST_FAIL_SUBMIT=0 EM_TEST_FAIL_CONTINUATION=1
if "$BASH" "$repo_dir/job/em_susie_mix" 3 > "${test_dir}/out" 2>&1; then
    echo 'Expected continuation submission failure to be reported' >&2; exit 1
fi
grep -q 'After the array finishes, run: sbatch em_susie_mix 2' "${test_dir}/out"
grep -Fxq '4246' "${SUSIE_MIX_PROJECT_DIR}/results_em/last_array_job_id.txt"
export EM_TEST_FAIL_CONTINUATION=0

# A failed prior update must stop a multi-iteration request before either job.
export EM_TEST_FAIL_R=1
rm -f -- "$EM_TEST_ARGS" "$EM_TEST_CONT_ARGS"
if "$BASH" "$repo_dir/job/em_susie_mix" 2 > "${test_dir}/out" 2>&1; then
    echo 'Expected prior preparation failure to stop the chain' >&2; exit 1
fi
[[ ! -e "$EM_TEST_ARGS" && ! -e "$EM_TEST_CONT_ARGS" ]]

# The worker passes the project, exact iteration directory, and chunk ID.
export SLURM_ARRAY_TASK_ID=185
if "$BASH" "$repo_dir/job/em_susie_mix_array" "$SUSIE_MIX_PROJECT_DIR" \
    "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001"; then
    echo 'Expected worker R failure to propagate to Slurm' >&2; exit 1
fi
export EM_TEST_FAIL_R=0
"$BASH" "$repo_dir/job/em_susie_mix_array" "$SUSIE_MIX_PROJECT_DIR" \
    "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001"
grep -Fxq '185' "$EM_TEST_R_ARGS"
grep -Fxq "${SUSIE_MIX_PROJECT_DIR}/script/scan_tissue_attempt/run_em_chunk.R" "$EM_TEST_R_ARGS"
echo 'All EM launcher tests passed (mock Slurm/R commands).'
