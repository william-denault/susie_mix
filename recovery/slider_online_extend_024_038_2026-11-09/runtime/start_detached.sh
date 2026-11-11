#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
cd "$bundle"
sha256sum -c SHA256SUMS >/dev/null
args=(--handoff)
if [[ -e /project2/mstephens/wdenault/susie_mix/workflows/slider_online_015_025/config.json ]]; then
    args+=(--resume)
fi
nohup bash "$bundle/start.sh" "${args[@]}" >> "$bundle/handoff_launch.log" 2>&1 < /dev/null &
echo "Started detached handoff launcher PID $!."
echo "This is a launcher PID, not a Slurm job ID."
echo "Log: $bundle/handoff_launch.log"
