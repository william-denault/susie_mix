#!/bin/bash
set -euo pipefail
patch=$(cd "$(dirname "$0")" && pwd -P)
target=${SLIDER_INPLACE_BUNDLE:-$HOME/slider_inplace_2026-10-29}
project=${SUSIE_MIX_PROJECT_DIR:-/project2/mstephens/wdenault/susie_mix}
state=${SLIDER_SUBSET_PATCH_STATE:-$HOME/slider_subset_patch_state_2026-10-30}
(cd "$patch" && sha256sum -c SHA256SUMS >/dev/null)
python3 "$patch/install_subset.py" "$project" "$target" "$state"
exec bash "$target/start_005_to_040.sh" --resume
