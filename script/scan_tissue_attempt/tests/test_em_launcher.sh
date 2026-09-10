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
export EM_TEST_R_ARGS="${test_dir}/r_args.txt"
export EM_TEST_FAIL_R=0 EM_TEST_ACTIVE=0 EM_TEST_LOCKED=0 EM_TEST_FAIL_SUBMIT=0
export USER="${USER:-em_test}"

module() { :; }
flock() { [[ "$EM_TEST_LOCKED" == 0 ]]; }
squeue() { if [[ "$EM_TEST_ACTIVE" == 1 ]]; then printf '4242\n'; fi; }
Rscript() {
    printf '%s\n' "$@" > "$EM_TEST_R_ARGS"
    if [[ "$EM_TEST_FAIL_R" == 1 ]]; then return 1; fi
    if [[ "$*" == *prepare_em_iteration.R* ]]; then
        mkdir -p "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001"
        printf '%s\n' "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001"
        local i sep=""
        for ((i=1; i<=185; i++)); do printf '%s%d' "$sep" "$i"; sep=","; done
        printf '\n'
    fi
}
sbatch() {
    printf '%s\n' "$@" > "$EM_TEST_ARGS"
    if [[ "$EM_TEST_FAIL_SUBMIT" == 1 ]]; then return 1; fi
    printf '4242;test_cluster\n'
}
export -f module flock squeue Rscript sbatch

"$BASH" "$repo_dir/job/em_susie_mix" > "${test_dir}/out" 2>&1
grep -q -- '--array=1,2,3,.*184,185' "$EM_TEST_ARGS"
grep -Fxq "${SUSIE_MIX_PROJECT_DIR}/job/em_susie_mix_array" "$EM_TEST_ARGS"
grep -Fxq "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001" "$EM_TEST_ARGS"
grep -Fxq '4242' "${SUSIE_MIX_PROJECT_DIR}/results_em/last_array_job_id.txt"
grep -Fxq 'new' "$EM_TEST_R_ARGS"

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

# The worker passes the project, exact iteration directory, and chunk ID.
export SLURM_ARRAY_TASK_ID=185
"$BASH" "$repo_dir/job/em_susie_mix_array" "$SUSIE_MIX_PROJECT_DIR" \
    "${SUSIE_MIX_PROJECT_DIR}/results_em/iteration_001"
grep -Fxq '185' "$EM_TEST_R_ARGS"
grep -Fxq "${SUSIE_MIX_PROJECT_DIR}/script/scan_tissue_attempt/run_em_chunk.R" "$EM_TEST_R_ARGS"
echo 'All EM launcher tests passed (mock Slurm/R commands).'
