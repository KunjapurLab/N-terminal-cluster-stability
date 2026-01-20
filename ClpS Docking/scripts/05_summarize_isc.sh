#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/05_summarize_isc.sh runs 10
RUNS_DIR="${1:?Need runs directory}"
K="${2:-10}"

for d in "$RUNS_DIR"/*; do
  [[ -d "$d" ]] || continue
  sc="$(ls "$d"/*.sc 2>/dev/null | head -n 1)"
  [[ -f "$sc" ]] || continue

  # Assumes I_sc is last column in SCORE lines (matches your workflow).
  topk_mean="$(awk '$1=="SCORE:"{print $NF}' "$sc" | sort -n | head -n "$K" | awk '{s+=$1} END{if(NR>0) print s/NR; else print "NA"}')"
  minv="$(awk '$1=="SCORE:"{print $NF}' "$sc" | sort -n | head -n 1)"
  maxv="$(awk '$1=="SCORE:"{print $NF}' "$sc" | sort -n | tail -n 1)"

  echo "$(basename "$d")  top${K}_mean_Isc=${topk_mean}  min=${minv}  max=${maxv}"
done
