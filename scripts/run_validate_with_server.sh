#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

port="${PORT:-4173}"
log_file="${TMPDIR:-/tmp}/lookup-page-http.log"

python -m http.server "$port" --bind 127.0.0.1 >"$log_file" 2>&1 &
server_pid="$!"

cleanup() {
  kill "$server_pid" >/dev/null 2>&1 || true
}
trap cleanup EXIT

sleep 1
python scripts/validate_page.py
