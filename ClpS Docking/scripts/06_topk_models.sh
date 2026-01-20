#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/06_topk_models.sh runs 5
RUNS_DIR="${1:?Need runs directory}"
K="${2:-5}"

for d in "$RUNS_DIR"/*; do
  [[ -d "$d" ]] || continue
  sc="$(ls "$d"/*.sc 2>/dev/null | head -n 1)"
  [[ -f "$sc" ]] || continue

  echo "## $(basename "$d")"
  # Print: Isc decoy_tag (assuming tag is column 2 on SCORE lines)
  awk '$1=="SCORE:"{print $NF,$2}' "$sc" | sort -n | head -n "$K"
  echo
done
