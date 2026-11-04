#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
project=${SUSIE_MIX_PROJECT_DIR:-/project2/mstephens/wdenault/susie_mix}
original=${SLIDER_STORAGE_STATE:-$HOME/slider_storage_state_2026-10-28}
scratch=${SLIDER_STORAGE_SCRATCH:-/scratch/midway2/$USER/slider_storage_2026-10-28}
partial=${SLIDER_PARTIAL_STATE:-$HOME/slider_partial_finalize_2026-10-29}
launch=${SLIDER_INPLACE_STATE:-$HOME/slider_inplace_launch_2026-10-29}
[[ $# -eq 0 || ( $# -eq 1 && "$1" == --resume ) ]] || { echo 'Usage: start_005_to_040.sh [--resume]' >&2; exit 1; }
(cd "$bundle" && sha256sum -c SHA256SUMS >/dev/null)
exec python3 "$bundle/orchestrate.py" "$project" "$original" "$scratch" "$partial" "$launch" "$@"
