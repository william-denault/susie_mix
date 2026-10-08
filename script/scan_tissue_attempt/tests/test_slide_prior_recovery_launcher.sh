#!/bin/bash
# All Slurm and R commands are mocked; no network, jobs or actual fits are touched.
set -euo pipefail
repo=$(pwd -P)
mkdir -p "$repo/tmp"
fixture=$(mktemp -d "$repo/tmp/slide_recovery_launcher.XXXXXX")
cleanup() {
    local resolved
    resolved=$(cd "$fixture" && pwd -P)
    if [[ "$resolved" == "$repo"/tmp/slide_recovery_launcher.* ]]; then rm -rf -- "$resolved"; fi
}
trap cleanup EXIT
export SUSIE_MIX_PROJECT_DIR="$fixture/project with spaces"
export REC_TEST="$fixture" SLURM_JOB_ID=9000 USER="${USER:-test_user}"
export REC_OTHERS=0 REC_ACTIVE=0 REC_ACTIVE_NAME=0 REC_STATUS=batch REC_LOCK=0
export REC_FAIL_ARRAY=0 REC_FAIL_NEXT=0 REC_FAIL_GENES=0 REC_FAIL_GENE=0 REC_GENES_COUNT=2
export REC_PLAN_CALLS=0
unset SUSIE_SLIDE_GENES_PER_TASK
mkdir -p "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em"
controller="$repo/job/recover_slide_prior_em"
worker="$repo/job/recover_slide_prior_em_array"
module() { :; }
flock() { [[ "$REC_LOCK" == 0 ]]; }
squeue() {
    [[ "$*" == *--array* && "$*" == *"%i|%F|%j"* ]] || return 1
    printf '9000|9000|slide_recover_control\n'
    local i
    for ((i=1; i<=REC_OTHERS; i++)); do printf '%s|7000|unrelated\n' "7000_$i"; done
    if [[ "$REC_ACTIVE" == 1 ]]; then printf '9100_1|9100|worker\n'; fi
    if [[ "$REC_ACTIVE_NAME" == 1 ]]; then printf '9200_1|9200|slide_recover_worker\n'; fi
}
Rscript() {
    local mode=$3
    printf '%s\n' "$mode" >> "$REC_TEST/r_calls"
    case "$mode" in
        audit|verify|exclude-sex-mt) return 0 ;;
        new) printf '1\n' ;;
        plan)
            printf '%s\n' "$6" > "$REC_TEST/capacity"
            printf '%s\n' "$SUSIE_SLIDE_GENES_PER_TASK" > "$REC_TEST/gene_cap"
            printf '%s\n' "$SUSIE_SLIDE_AUTOMATIC" > "$REC_TEST/automatic"
            if [[ "$REC_STATUS" != batch ]]; then printf '%s\n' "$REC_STATUS"; return; fi
            local batch="$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/iteration_000/recovery/batches/batch_test"
            mkdir -p "$batch"
            printf 'batch\n%s\n%s\n' "$batch" "$6" ;;
        genes)
            [[ "$REC_FAIL_GENES" == 0 ]] || return 1
            local i
            for ((i=1; i<=REC_GENES_COUNT; i++)); do printf 'G%03d\n' "$i"; done ;;
        gene)
            printf '%s\n' "$8" >> "$REC_TEST/gene_calls"
            [[ "$REC_FAIL_GENE" == 0 || "$8" != G001 ]] ;;
        *) return 1 ;;
    esac
}
sbatch() {
    if [[ "$*" == *--array=* ]]; then
        printf '%s\n' "$@" > "$REC_TEST/array_args"
        printf '%s\n' "$SUSIE_SLIDE_GENES_PER_TASK" > "$REC_TEST/array_gene_cap"
        [[ "$REC_FAIL_ARRAY" == 0 ]] || return 1
        printf '9100;cluster\n'
    else
        printf '%s\n' "$@" > "$REC_TEST/next_args"
        printf '%s\n' "$SUSIE_SLIDE_GENES_PER_TASK" > "$REC_TEST/next_gene_cap"
        [[ "$REC_FAIL_NEXT" == 0 ]] || return 1
        printf '9101;cluster\n'
    fi
}
export -f module flock squeue Rscript sbatch
reset_calls() { : > "$REC_TEST/r_calls"; : > "$REC_TEST/array_args"; : > "$REC_TEST/next_args"; }
expect_fail() { if "$@" > "$REC_TEST/output" 2>&1; then echo 'Expected failure' >&2; exit 1; fi; }
reset_calls
"$BASH" "$controller" audit 0 > "$REC_TEST/output" 2>&1
[[ ! -s "$REC_TEST/array_args" ]]
grep -Fxq audit "$REC_TEST/r_calls"
reset_calls
"$BASH" "$controller" exclude-sex-mt 0 > "$REC_TEST/output" 2>&1
grep -Fxq exclude-sex-mt "$REC_TEST/r_calls"
[[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/next_args" ]]

"$BASH" "$controller" run 0 > "$REC_TEST/output" 2>&1
grep -Fxq '298' "$REC_TEST/capacity"
grep -Fxq -- '--array=1-298' "$REC_TEST/array_args"
grep -Fxq -- '--dependency=afterany:9100:9000' "$REC_TEST/next_args"
grep -Fxq continue "$REC_TEST/next_args"
grep -Fxq '9100' "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/last_array_job_id.txt"
grep -Fxq 0 "$REC_TEST/automatic"
grep -Fxq 4 "$REC_TEST/gene_cap"
grep -Fxq 4 "$REC_TEST/array_gene_cap"
grep -Fxq 4 "$REC_TEST/next_gene_cap"
"$BASH" "$controller" continue 0 > "$REC_TEST/output" 2>&1
grep -Fxq 1 "$REC_TEST/automatic"
# The optional five-gene cap reaches both workers and subsequent coordinators.
export SUSIE_SLIDE_GENES_PER_TASK=5
"$BASH" "$controller" resume 0 > "$REC_TEST/output" 2>&1
grep -Fxq 5 "$REC_TEST/gene_cap"
grep -Fxq 5 "$REC_TEST/array_gene_cap"
grep -Fxq 5 "$REC_TEST/next_gene_cap"
grep -Fq 'at most 5 genes/task' "$REC_TEST/output"
for cap in 0 16 4.5 -1 bad 04; do
    reset_calls
    expect_fail env SUSIE_SLIDE_GENES_PER_TASK="$cap" "$BASH" "$controller" resume 0
    [[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/r_calls" ]]
done
unset SUSIE_SLIDE_GENES_PER_TASK

# Count every other submitted array element, plus this and the next controller.
export REC_OTHERS=295
"$BASH" "$controller" resume 0 > "$REC_TEST/output" 2>&1
grep -Fxq '3' "$REC_TEST/capacity"
export REC_OTHERS=298
reset_calls
expect_fail "$BASH" "$controller" resume 0
[[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/r_calls" ]]
export REC_OTHERS=0

# Exact IDs and job names both protect against duplicate submissions.
for flag in REC_ACTIVE REC_ACTIVE_NAME REC_LOCK; do
    export "$flag=1"
    reset_calls
    expect_fail "$BASH" "$controller" resume 0
    [[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/r_calls" ]]
    expect_fail "$BASH" "$controller" exclude-sex-mt 0
    [[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/r_calls" ]]
    export "$flag=0"
done
for args in 'run 0 299' 'run 0 0' 'run -1' 'bad 0' 'run 01'; do
    expect_fail "$BASH" "$controller" $args
done
expect_fail env -u SLURM_JOB_ID "$BASH" "$controller" run 0

# Worker failures cannot cancel inspection; blocked retries cannot advance EM.
export REC_STATUS=blocked
reset_calls
expect_fail "$BASH" "$controller" resume 0
[[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/next_args" ]]
export REC_STATUS=verify
reset_calls
"$BASH" "$controller" resume 0 > "$REC_TEST/output" 2>&1
grep -Fxq verify "$REC_TEST/r_calls"
[[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/next_args" ]]
export REC_STATUS=batch

export REC_FAIL_ARRAY=1
reset_calls
expect_fail "$BASH" "$controller" run 0
[[ ! -s "$REC_TEST/next_args" ]]
export REC_FAIL_ARRAY=0 REC_FAIL_NEXT=1
reset_calls
expect_fail "$BASH" "$controller" run 0
grep -q 'Array 9100 remains submitted' "$REC_TEST/output"
export REC_FAIL_NEXT=0 SUSIE_SLIDE_PILOT=1
reset_calls
"$BASH" "$controller" run 0 5 > "$REC_TEST/output" 2>&1
grep -Fxq -- '--array=1-5' "$REC_TEST/array_args"
[[ ! -s "$REC_TEST/next_args" ]]
unset SUSIE_SLIDE_PILOT

# Every gene uses a separate Rscript invocation. A failed gene does not erase
# another gene's successful work, and enumeration failures cannot look complete.
export SLURM_ARRAY_TASK_ID=1 REC_FAIL_GENE=1
batch="$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/iteration_000/recovery/batches/batch_test"
expect_fail "$BASH" "$worker" "$SUSIE_MIX_PROJECT_DIR" 0 "$batch"
[[ "$(cat "$REC_TEST/gene_calls")" == $'G001\nG002' ]]
export REC_FAIL_GENE=0 REC_FAIL_GENES=1
: > "$REC_TEST/gene_calls"
expect_fail "$BASH" "$worker" "$SUSIE_MIX_PROJECT_DIR" 0 "$batch"
[[ ! -s "$REC_TEST/gene_calls" ]]
export REC_FAIL_GENES=0 REC_GENES_COUNT=16
expect_fail "$BASH" "$worker" "$SUSIE_MIX_PROJECT_DIR" 0 "$batch"
[[ ! -s "$REC_TEST/gene_calls" ]]
# Adopting recovery prevents accidentally resubmitting the old large chunks.
reset_calls
touch "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/iteration_000/recovery/state.rds"
expect_fail "$BASH" "$repo/job/em_susie_slide_prior" resume
[[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/r_calls" ]]
echo 'PASS: account capacity, configurable gene caps and propagation, array bounds, duplicate protection, afterany inspection, pilot, submission failures and isolated gene commands.'
