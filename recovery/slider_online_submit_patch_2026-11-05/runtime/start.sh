#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
cd "$bundle"
sha256sum -c SHA256SUMS >/dev/null
exec python3 "$bundle/minibatch.py" launch "$@"
