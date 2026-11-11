#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
cd "$bundle"
sha256sum -c SHA256SUMS >/dev/null
python3 -u upgrade.py launch
