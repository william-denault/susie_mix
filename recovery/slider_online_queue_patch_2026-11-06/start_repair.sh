#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
cd "$bundle"
sha256sum -c SHA256SUMS >/dev/null
mkdir -p logs
sbatch --parsable --exclude=midway2-0003,midway2-0011 \
  --chdir="$bundle" --output="$bundle/logs/repair_%j.out" \
  --error="$bundle/logs/repair_%j.err" "$bundle/repair_job" "$bundle"
