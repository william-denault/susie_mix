#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
project=${1:-/project2/mstephens/wdenault/susie_mix}
original_state=${2:-$HOME/slider_storage_state_2026-10-28}
scratch=${3:-/scratch/midway2/$USER/slider_storage_2026-10-28}
state=${SLIDER_PARTIAL_STATE:-$HOME/slider_partial_finalize_2026-10-29}
[[ $# -le 3 ]] || { echo 'Usage: stop_compaction_resume_005.sh [PROJECT [ORIGINAL_STATE [SCRATCH]]]' >&2; exit 1; }
(cd "$bundle" && sha256sum -c SHA256SUMS >/dev/null)
exec python3 "$bundle/handoff.py" "$project" "$original_state" "$scratch" "$state"
