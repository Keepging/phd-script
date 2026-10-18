#!/usr/bin/env bash
# Usage: bash run.sh [out_dir]   (needs network access to rest.uniprot.org)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$HERE/output}"; mkdir -p "$OUT"
if [ ! -d "$HERE/.venv" ]; then
  uv venv -p python3.12 "$HERE/.venv" && uv pip install -p "$HERE/.venv/bin/python" -r "$HERE/requirements.txt"
fi
PY="$HERE/.venv/bin/python"
"$PY" "$HERE/fetch_uniprot.py" "$HERE/input/site_meta.tsv" "$OUT"
"$PY" "$HERE/score_sites.py" "$HERE/input/site_meta.tsv" "$OUT/uniprot_sequences.tsv" "$OUT"
