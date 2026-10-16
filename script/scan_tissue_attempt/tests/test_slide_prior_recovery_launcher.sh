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
export REC_FAIL_VERIFY=0
export REC_FAIL_GRANT=0
export REC_ROLLING_WORKERS=0 REC_CONTROL=0 REC_QUEUE_FAIL=0 REC_QUEUE_EMPTY=0 REC_DONE_CONTROL=0
unset SUSIE_SLIDE_GENES_PER_TASK SUSIE_SLIDE_TARGET_ITERATION SUSIE_SLIDE_PILOT SUSIE_SLIDE_EXCLUDE_NODES SLURM_RESTART_COUNT SUSIE_SLIDE_ROLLING SUSIE_SLIDE_REFILL_BELOW SUSIE_SLIDE_POLL_MINUTES
mkdir -p "$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em"
controller="$repo/job/recover_slide_prior_em"
worker="$repo/job/recover_slide_prior_em_array"
module() { :; }
flock() { [[ "$REC_LOCK" == 0 ]]; }
squeue() {
    [[ "$*" == *--array* && "$*" == *"%i|%F|%j"* ]] || return 1
    printf '%s\n' "$*" > "$REC_TEST/queue_args"
    [[ "$REC_QUEUE_FAIL" == 0 ]] || return 1
    [[ "$REC_QUEUE_EMPTY" == 0 ]] || return 0
    qrow() {
        if [[ "${SUSIE_SLIDE_ROLLING:-0}" == 1 ]]; then printf '%s|%s\n' "$1" "${2:-RUNNING}";
        else printf '%s\n' "$1"; fi
    }
    qrow '9000|9000|slide_recover_control'
    local i
    for ((i=1; i<=REC_OTHERS; i++)); do qrow "7000_$i|7000|unrelated"; done
    if [[ "$REC_ACTIVE" == 1 ]]; then qrow '9100_1|9100|worker'; fi
    if [[ "$REC_ACTIVE_NAME" == 1 ]]; then qrow '9200_1|9200|slide_recover_worker'; fi
    for ((i=1; i<=REC_ROLLING_WORKERS; i++)); do qrow "9100_$i|9100|slide_recover_worker"; done
    if [[ "$REC_CONTROL" == 1 ]]; then qrow '9300|9300|slide_recover_control'; fi
    if [[ "$REC_DONE_CONTROL" == 1 ]]; then qrow '9400|9400|slide_recover_control' COMPLETED; fi
}
Rscript() {
    local mode=$3
    printf '%s\n' "$mode" >> "$REC_TEST/r_calls"
    case "$mode" in
        audit|exclude-sex-mt) return 0 ;;
        grant-retry)
            printf '%s\n' "$6" "$7" "$8" > "$REC_TEST/grant_args"
            [[ "$REC_FAIL_GRANT" == 0 ]] ;;
        verify) [[ "$REC_FAIL_VERIFY" == 0 ]] ;;
        new) printf '1\n' ;;
        plan)
            if [[ "${SUSIE_SLIDE_ROLLING:-0}" == 1 ]]; then
                [[ $# == 7 && -f "$7" ]] || return 1
                cp "$7" "$REC_TEST/queue_snapshot"
            else
                [[ $# == 6 ]] || return 1
            fi
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
            printf 'fixture gene output\n'
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
        printf '%s\n' "${SUSIE_SLIDE_TARGET_ITERATION:-}" > "$REC_TEST/next_target"
        printf '%s\n' "${SUSIE_SLIDE_ROLLING:-}" "${SUSIE_SLIDE_REFILL_BELOW:-}" "${SUSIE_SLIDE_POLL_MINUTES:-}" > "$REC_TEST/next_rolling"
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

# Reviewed retries run under the same lock/queue guards, before planning.
export SUSIE_SLIDE_EXCLUDE_NODES=midway2-0003
reset_calls
"$BASH" "$controller" retry 0 G001 3 'Three reviewed NODE_FAIL records' > "$REC_TEST/output" 2>&1
[[ "$(cat "$REC_TEST/r_calls")" == $'grant-retry\nplan' ]]
[[ "$(cat "$REC_TEST/grant_args")" == $'G001\n3\nThree reviewed NODE_FAIL records' ]]
grep -Fxq -- '--exclude=midway2-0003' "$REC_TEST/array_args"
grep -Fxq -- '--exclude=midway2-0003' "$REC_TEST/next_args"
export REC_FAIL_GRANT=1
reset_calls
expect_fail "$BASH" "$controller" retry 0 G001 3 'Three reviewed NODE_FAIL records'
[[ "$(cat "$REC_TEST/r_calls")" == grant-retry && ! -s "$REC_TEST/array_args" ]]
export REC_FAIL_GRANT=0 REC_ACTIVE_NAME=1
reset_calls
expect_fail "$BASH" "$controller" retry 0 G001 3 'Three reviewed NODE_FAIL records'
[[ ! -s "$REC_TEST/r_calls" && ! -s "$REC_TEST/array_args" ]]
export REC_ACTIVE_NAME=0
for args in 'retry 0' 'retry 0 G001 3' 'retry 0 G001 bad reason'; do
    expect_fail "$BASH" "$controller" $args
done
expect_fail env SUSIE_SLIDE_EXCLUDE_NODES='node one' "$BASH" "$controller" resume 0
unset SUSIE_SLIDE_EXCLUDE_NODES

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
# A verified iteration queues exactly one next controller up to the target.
export SUSIE_SLIDE_TARGET_ITERATION=20
export SUSIE_SLIDE_EXCLUDE_NODES=midway2-0003
reset_calls
"$BASH" "$controller" resume 0 > "$REC_TEST/output" 2>&1
grep -Fxq verify "$REC_TEST/r_calls"
grep -Fxq -- '--dependency=afterok:9000' "$REC_TEST/next_args"
grep -Fxq -- '--exclude=midway2-0003' "$REC_TEST/next_args"
grep -Fxq new "$REC_TEST/next_args"
grep -Fxq 1 "$REC_TEST/next_args"
grep -Fxq 20 "$REC_TEST/next_target"
[[ ! -s "$REC_TEST/array_args" ]]
reset_calls
"$BASH" "$controller" continue 20 > "$REC_TEST/output" 2>&1
[[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/next_args" ]]
export REC_FAIL_VERIFY=1
reset_calls
expect_fail "$BASH" "$controller" continue 19
[[ ! -s "$REC_TEST/next_args" ]]
export REC_FAIL_VERIFY=0 SUSIE_SLIDE_PILOT=1
reset_calls
"$BASH" "$controller" resume 0 > "$REC_TEST/output" 2>&1
[[ ! -s "$REC_TEST/next_args" ]]
unset SUSIE_SLIDE_PILOT
export REC_FAIL_NEXT=1
reset_calls
expect_fail "$BASH" "$controller" continue 19
grep -q 'next-iteration submission failed' "$REC_TEST/output"
export REC_FAIL_NEXT=0
for target in -1 01 bad 1.5 999999999; do
    reset_calls
    expect_fail env SUSIE_SLIDE_TARGET_ITERATION="$target" "$BASH" "$controller" resume 0
    [[ ! -s "$REC_TEST/r_calls" && ! -s "$REC_TEST/next_args" ]]
done
expect_fail env SUSIE_SLIDE_TARGET_ITERATION=19 "$BASH" "$controller" continue 20
# Exercise every completed boundary: 000 -> 001 through 019 -> 020, then stop.
for ((i=0; i<=20; i++)); do
    reset_calls
    "$BASH" "$controller" continue "$i" > "$REC_TEST/output" 2>&1
    [[ ! -s "$REC_TEST/array_args" ]]
    if (( i < 20 )); then
        grep -Fxq new "$REC_TEST/next_args"
        grep -Fxq "$((i + 1))" "$REC_TEST/next_args"
    else
        [[ ! -s "$REC_TEST/next_args" ]]
    fi
done
export REC_STATUS=batch
reset_calls
"$BASH" "$controller" continue 7 > "$REC_TEST/output" 2>&1
grep -Fxq continue "$REC_TEST/next_args"
grep -Fxq 20 "$REC_TEST/next_target"
grep -Fxq 4 "$REC_TEST/next_gene_cap"
unset SUSIE_SLIDE_TARGET_ITERATION
unset SUSIE_SLIDE_EXCLUDE_NODES
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

# Rolling recovery adopts live workers, counts them against the global budget,
# and schedules its next check without an array-wide dependency.
export SUSIE_SLIDE_ROLLING=1 SUSIE_SLIDE_TARGET_ITERATION=20
export SUSIE_SLIDE_EXCLUDE_NODES=midway2-0003 REC_ROLLING_WORKERS=2
export REC_DONE_CONTROL=1
reset_calls
"$BASH" "$controller" resume 0 > "$REC_TEST/output" 2>&1
grep -Fxq 296 "$REC_TEST/capacity"
grep -Fxq -- '--array=1-296' "$REC_TEST/array_args"
grep -Fxq -- '--dependency=afterany:9000' "$REC_TEST/next_args"
grep -Fxq -- '--begin=now+5minutes' "$REC_TEST/next_args"
grep -Fq -- '--states=all' "$REC_TEST/queue_args"
grep -Fxq '9100_1|9100|slide_recover_worker' "$REC_TEST/queue_snapshot"
! grep -Fq '9400' "$REC_TEST/queue_snapshot"
[[ "$(cat "$REC_TEST/next_rolling")" == $'1\n200\n5' ]]
grep -Fxq 20 "$REC_TEST/next_target"
grep -Fxq -- '--exclude=midway2-0003' "$REC_TEST/array_args"
grep -Fxq -- '--exclude=midway2-0003' "$REC_TEST/next_args"
export REC_OTHERS=10
"$BASH" "$controller" continue 0 > "$REC_TEST/output" 2>&1
grep -Fxq 286 "$REC_TEST/capacity"
export REC_OTHERS=0
"$BASH" "$controller" continue 0 100 > "$REC_TEST/output" 2>&1
grep -Fxq 98 "$REC_TEST/capacity"
# Full queues and unfinished genes are waiting states, not failed runs.
export REC_STATUS=waiting REC_OTHERS=296
reset_calls
"$BASH" "$controller" continue 0 > "$REC_TEST/output" 2>&1
grep -Fxq 0 "$REC_TEST/capacity"
[[ ! -s "$REC_TEST/array_args" && -s "$REC_TEST/next_args" ]]
grep -Fxq -- '--dependency=afterany:9000' "$REC_TEST/next_args"
export REC_OTHERS=0 REC_ROLLING_WORKERS=210
reset_calls
"$BASH" "$controller" continue 0 > "$REC_TEST/output" 2>&1
[[ ! -s "$REC_TEST/array_args" && -s "$REC_TEST/next_args" ]]
export REC_ROLLING_WORKERS=2 REC_STATUS=verify
reset_calls
expect_fail "$BASH" "$controller" continue 0
[[ ! -s "$REC_TEST/next_args" ]]
! grep -Fxq verify "$REC_TEST/r_calls"
# A second coordinator, a failed/empty query, or a lock must stop all submissions.
export REC_STATUS=batch
for flag in REC_CONTROL REC_QUEUE_FAIL REC_QUEUE_EMPTY REC_LOCK; do
    export "$flag=1"
    reset_calls
    expect_fail "$BASH" "$controller" resume 0
    [[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/next_args" && ! -s "$REC_TEST/r_calls" ]]
    export "$flag=0"
done
for mode in audit verify new exclude-sex-mt; do
    reset_calls
    expect_fail "$BASH" "$controller" "$mode" 0
    [[ ! -s "$REC_TEST/r_calls" ]]
done
expect_fail "$BASH" "$controller" retry 0 G001 3 'Node failure'
for setting in SUSIE_SLIDE_ROLLING=2 SUSIE_SLIDE_REFILL_BELOW=0 SUSIE_SLIDE_REFILL_BELOW=299 SUSIE_SLIDE_POLL_MINUTES=0 SUSIE_SLIDE_POLL_MINUTES=61 SUSIE_SLIDE_PILOT=1; do
    reset_calls
    expect_fail env "$setting" "$BASH" "$controller" resume 0
    [[ ! -s "$REC_TEST/r_calls" ]]
done
# Rolling configuration survives an EM boundary, which retains afterok.
export REC_STATUS=verify REC_ROLLING_WORKERS=0
reset_calls
"$BASH" "$controller" continue 0 > "$REC_TEST/output" 2>&1
grep -Fxq -- '--dependency=afterok:9000' "$REC_TEST/next_args"
grep -Fxq new "$REC_TEST/next_args"
[[ "$(cat "$REC_TEST/next_rolling")" == $'1\n200\n5' ]]
export REC_STATUS=blocked
reset_calls
expect_fail "$BASH" "$controller" continue 0
[[ ! -s "$REC_TEST/array_args" && ! -s "$REC_TEST/next_args" ]]
export REC_STATUS=batch
export REC_DONE_CONTROL=0
unset SUSIE_SLIDE_ROLLING SUSIE_SLIDE_TARGET_ITERATION SUSIE_SLIDE_EXCLUDE_NODES

# Every gene uses a separate Rscript invocation. A failed gene does not erase
# another gene's successful work, and enumeration failures cannot look complete.
export SLURM_ARRAY_TASK_ID=1 REC_FAIL_GENE=1
batch="$SUSIE_MIX_PROJECT_DIR/results_slide_prior_em/iteration_000/recovery/batches/batch_test"
expect_fail "$BASH" "$worker" "$SUSIE_MIX_PROJECT_DIR" 0 "$batch"
[[ "$(cat "$REC_TEST/gene_calls")" == $'G001\nG002' ]]
# A requeue must retain the previous gene logs, including when run repeatedly.
log="$batch/gene_logs/G001_9000_restart_0.out"
[[ "$(cat "$log")" == 'fixture gene output' ]]
export SLURM_RESTART_COUNT=1 REC_FAIL_GENE=0
"$BASH" "$worker" "$SUSIE_MIX_PROJECT_DIR" 0 "$batch" > "$REC_TEST/output" 2>&1
[[ "$(cat "$log")" == 'fixture gene output' ]]
[[ "$(cat "$batch/gene_logs/G001_9000_restart_1.out")" == 'fixture gene output' ]]
"$BASH" "$worker" "$SUSIE_MIX_PROJECT_DIR" 0 "$batch" > "$REC_TEST/output" 2>&1
[[ "$(cat "$batch/gene_logs/G001_9000_restart_1.out")" == $'fixture gene output\nfixture gene output' ]]
grep -Fxq '#SBATCH --open-mode=append' "$worker"
grep -Fxq '#SBATCH --open-mode=append' "$controller"
unset SLURM_RESTART_COUNT
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
echo 'PASS: reviewed retry guards, node exclusion propagation, restart log preservation, capacity, gene caps, 20-update chain, duplicate protection, pilot and submission failures.'
