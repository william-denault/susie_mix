#!/bin/bash
set -euo pipefail
bundle=$(cd "$(dirname "$0")" && pwd -P)
exec python3 "$bundle/minibatch.py" status "$@"
