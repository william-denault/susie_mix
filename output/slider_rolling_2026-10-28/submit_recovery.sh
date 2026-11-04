#!/bin/bash
# Submit retirement (optional), compaction, then a resume-to-005 gate.
set -euo pipefail
bundle=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
project=${1:?Usage: submit_recovery.sh PROJECT [--retire-through-003]}
mode=${2:-keep-all}
[[ "$mode" == keep-all || "$mode" == --retire-through-003 ]] || exit 2
state="${SLIDER_STORAGE_STATE:-$HOME/slider_storage_state_2026-10-28}"
retired="${SLIDER_RETIRE_STATE:-$HOME/slider_retire_000_003_2026-10-28}"
scratch="${SLIDER_STORAGE_SCRATCH:-/scratch/midway2/$USER/slider_storage_2026-10-28}"
mkdir -p "$state/logs"
[[ ! -e "$state/queued_jobs.tsv" ]] || { echo "Submission receipt exists: $state/queued_jobs.tsv. See README for phase-specific resume commands." >&2; exit 1; }
(cd "$bundle" && sha256sum -c SHA256SUMS >/dev/null)
common=(--parsable --kill-on-invalid-dep=yes --exclude=midway2-0003,midway2-0011)
dependency=()
if [[ "$mode" == --retire-through-003 ]]; then
  python3 "$bundle/storage.py" retire "$project" "$retired" --through 3
  job=$(sbatch "${common[@]}" --cpus-per-task=1 --mem=16G --job-name=slide_storage_retire \
    --output="$state/logs/retire_%j.out" --error="$state/logs/retire_%j.err" \
    "$bundle/storage_job" "$bundle" retire "$project" "$retired" --through 3 --apply)
  job=${job%%;*}; [[ "$job" =~ ^[0-9]+$ ]] || exit 1
  printf 'retire\t%s\n' "$job" >> "$state/queued_jobs.tsv"
  dependency=(--dependency="afterok:$job")
fi
job=$(sbatch "${common[@]}" "${dependency[@]}" --job-name=slide_storage_compact \
  --output="$state/logs/compact_%j.out" --error="$state/logs/compact_%j.err" \
  "$bundle/storage_job" "$bundle" compact "$project" "$state" "$scratch" --workers 2)
job=${job%%;*}; [[ "$job" =~ ^[0-9]+$ ]] || exit 1
printf 'compact\t%s\n' "$job" >> "$state/queued_jobs.tsv"
resume=$(sbatch "${common[@]}" --dependency="afterok:$job" --cpus-per-task=1 --mem=2G --time=00:10:00 \
  --job-name=slide_storage_resume --output="$state/logs/resume_%j.out" --error="$state/logs/resume_%j.err" \
  "$bundle/storage_job" "$bundle" resume "$project" "$state" --target 5)
resume=${resume%%;*}; [[ "$resume" =~ ^[0-9]+$ ]] || exit 1
printf 'resume\t%s\n' "$resume" >> "$state/queued_jobs.tsv"
cat "$state/queued_jobs.tsv"
printf 'State: %s\nScratch escrow: %s\n' "$state" "$scratch"
echo 'Iteration 005 will finish and stop before 006. Expansion is a separate gated step.'
