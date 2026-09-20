#!/bin/bash
# Run from the project root. Slurm, R, and flock are mocked; no jobs are sent.
set -euo pipefail
repo_dir=$(pwd -P)
test_base="$repo_dir/tmp/slide_simulation_validation"
mkdir -p "$test_base"
test_dir=$(mktemp -d "$test_base/launchers.XXXXXX")
export SUSIE_MIX_PROJECT_DIR="$test_dir/project with spaces"
mkdir -p "$SUSIE_MIX_PROJECT_DIR/script/sim/jobs_slide" "$SUSIE_MIX_PROJECT_DIR/job"
cp "$repo_dir/script/sim/jobs_slide/submission_batches.csv" "$SUSIE_MIX_PROJECT_DIR/script/sim/jobs_slide/"
cp "$repo_dir"/job/run_simulation_slide_batch_* "$SUSIE_MIX_PROJECT_DIR/job/"
cp "$repo_dir/job/launch_simulation_slide" "$SUSIE_MIX_PROJECT_DIR/job/"
for id in 1 300 301 600 601 900 901 1125; do
    touch "$SUSIE_MIX_PROJECT_DIR/script/sim/jobs_slide/sim_job_${id}.R"
done
export SLIDE_TEST_DIR="$test_dir" SLIDE_FAIL_R=0 SLIDE_FAIL_ARRAY=0 SLIDE_FAIL_CONT=0
export SLIDE_LOCKED=0 SLIDE_ACTIVE="" SLIDE_BATCH=1
export SLURM_JOB_ID=6000 USER="${USER:-slide_test}"
module() { :; }
flock() { [[ "$SLIDE_LOCKED" == 0 ]]; }
squeue() { printf '%s\n' "$SLIDE_ACTIVE"; }
Rscript() {
    printf '%s\n' "$@" > "$SLIDE_TEST_DIR/r_args"
    [[ "$SLIDE_FAIL_R" == 0 ]]
}
sbatch() {
    if [[ "$*" == *run_simulation_slide_batch_* ]]; then
        printf '%s\n' "$@" > "$SLIDE_TEST_DIR/array_args"
        printf 'array\n' >> "$SLIDE_TEST_DIR/calls"
        [[ "$SLIDE_FAIL_ARRAY" == 0 ]] || return 1
        printf '%s;test_cluster\n' "$((7000 + 2 * SLIDE_BATCH))"
    else
        printf '%s\n' "$@" > "$SLIDE_TEST_DIR/cont_args"
        printf 'continuation\n' >> "$SLIDE_TEST_DIR/calls"
        [[ "$SLIDE_FAIL_CONT" == 0 ]] || return 1
        printf '%s;test_cluster\n' "$((7001 + 2 * SLIDE_BATCH))"
    fi
}
export -f module flock squeue Rscript sbatch

# Every first/last task selects the right numbered R file, including batch 4.
for batch in 1 2 3 4; do
    last=299
    if (( batch == 4 )); then last=224; fi
    grep -Fxq "#SBATCH --array=0-$last" "$repo_dir/job/run_simulation_slide_batch_$batch"
    for task in 0 "$last"; do
        export SLURM_ARRAY_TASK_ID=$task
        "$BASH" "$repo_dir/job/run_simulation_slide_batch_$batch" > "$test_dir/out" 2>&1
        expected=$(( (batch - 1) * 300 + task + 1 ))
        grep -Fxq -- '--vanilla' "$test_dir/r_args"
        grep -Fxq "$SUSIE_MIX_PROJECT_DIR/script/sim/jobs_slide/sim_job_${expected}.R" "$test_dir/r_args"
    done
done
export SLURM_ARRAY_TASK_ID=225
if "$BASH" "$repo_dir/job/run_simulation_slide_batch_4" > "$test_dir/out" 2>&1; then
    echo 'Expected out-of-range job to fail' >&2; exit 1
fi
for invalid in -1 abc 01 400; do
    export SLURM_ARRAY_TASK_ID=$invalid
    if "$BASH" "$repo_dir/job/run_simulation_slide" > "$test_dir/out" 2>&1; then
        echo 'Expected invalid array task to fail' >&2; exit 1
    fi
done
export SLURM_ARRAY_TASK_ID=1
if "$BASH" "$repo_dir/job/run_simulation_slide" > "$test_dir/out" 2>&1; then
    echo 'Expected missing R script to fail' >&2; exit 1
fi
export SLURM_ARRAY_TASK_ID=0 SLIDE_FAIL_R=1
if "$BASH" "$repo_dir/job/run_simulation_slide" > "$test_dir/out" 2>&1; then
    echo 'Expected R failure to propagate to Slurm' >&2; exit 1
fi
export SLIDE_FAIL_R=0

# Four stages submit four arrays and only three small continuation jobs. Each
# continuation must depend on the whole preceding array AND its launcher.
for batch in 1 2 3 4; do
    export SLIDE_BATCH=$batch
    rm -f "$test_dir/cont_args"
    "$BASH" "$repo_dir/job/launch_simulation_slide" "$batch" > "$test_dir/out" 2>&1
    grep -Fxq "$SUSIE_MIX_PROJECT_DIR/job/run_simulation_slide_batch_$batch" "$test_dir/array_args"
    if (( batch < 4 )); then
        grep -Fxq -- "--dependency=afterok:$((7000 + 2 * batch)):$SLURM_JOB_ID" "$test_dir/cont_args"
        grep -Fxq -- '--kill-on-invalid-dep=yes' "$test_dir/cont_args"
        [[ "$(tail -n 1 "$test_dir/cont_args")" == "$((batch + 1))" ]]
        export SLURM_JOB_ID=$((7001 + 2 * batch)) SLIDE_ACTIVE=$((7001 + 2 * batch))
    else
        [[ ! -e "$test_dir/cont_args" ]]
    fi
done
[[ "$(grep -c '^array$' "$test_dir/calls")" == 4 ]]
[[ "$(grep -c '^continuation$' "$test_dir/calls")" == 3 ]]

# An active array or queued continuation blocks a duplicate manual launch.
export SLURM_JOB_ID=9000
for active in 7008 7007; do
    export SLIDE_ACTIVE=$active
    rm -f "$test_dir/array_args"
    if "$BASH" "$repo_dir/job/launch_simulation_slide" > "$test_dir/out" 2>&1; then
        echo 'Expected active work to block a duplicate launch' >&2; exit 1
    fi
    [[ ! -e "$test_dir/array_args" ]]
done
export SLIDE_ACTIVE="" SLIDE_LOCKED=1
if "$BASH" "$repo_dir/job/launch_simulation_slide" > "$test_dir/out" 2>&1; then
    echo 'Expected lock contention to stop submission' >&2; exit 1
fi
[[ ! -e "$test_dir/array_args" ]]
export SLIDE_LOCKED=0
for invalid in 0 -1 abc 01 5; do
    if "$BASH" "$repo_dir/job/launch_simulation_slide" "$invalid" > "$test_dir/out" 2>&1; then
        echo 'Expected invalid batch to fail' >&2; exit 1
    fi
done
if env -u SLURM_JOB_ID "$BASH" "$repo_dir/job/launch_simulation_slide" > "$test_dir/out" 2>&1; then
    echo 'Expected direct execution outside sbatch to fail' >&2; exit 1
fi

# Submission failures stop the chain and give an explicit recovery command.
export SLIDE_FAIL_ARRAY=1
rm -f "$test_dir/cont_args"
if "$BASH" "$repo_dir/job/launch_simulation_slide" 2 > "$test_dir/out" 2>&1; then
    echo 'Expected failed array submission to stop the chain' >&2; exit 1
fi
[[ ! -e "$test_dir/cont_args" ]]
grep -q 'sbatch job/launch_simulation_slide 2' "$test_dir/out"
export SLIDE_FAIL_ARRAY=0 SLIDE_FAIL_CONT=1 SLIDE_BATCH=2
if "$BASH" "$repo_dir/job/launch_simulation_slide" 2 > "$test_dir/out" 2>&1; then
    echo 'Expected continuation failure to be reported' >&2; exit 1
fi
grep -q 'sbatch job/launch_simulation_slide 3' "$test_dir/out"
grep -Fxq '7004' "$SUSIE_MIX_PROJECT_DIR/simulation results/slide_v1/launcher/last_array_job_id.txt"
echo 'PASS: batch boundaries, R failure propagation, four-stage chain, duplicate protection and submission failures (mock Slurm).'
