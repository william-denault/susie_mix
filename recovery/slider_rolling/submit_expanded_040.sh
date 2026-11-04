#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
project=${1:-/project2/mstephens/wdenault/susie_mix}
storage_state=${2:-$HOME/slider_storage_state_2026-10-28}
branch="$project/slider_expanded_040_2026-10-28"
launch_state="$HOME/slider_expanded_040_state_2026-10-28"
[[ $# -le 3 ]] || { echo 'Usage: submit_expanded_040.sh [PROJECT [STORAGE_STATE [--resume]]]' >&2; exit 1; }
extra=()
if [[ $# == 3 ]]; then
  [[ "$3" == --resume ]] || exit 1
  extra=(--resume)
fi
(cd "$bundle" && sha256sum -c SHA256SUMS >/dev/null)
exec python3 "$bundle/launch.py" submit "$project" "$storage_state" "$branch" "$launch_state" "${extra[@]}"
