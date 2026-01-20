#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/03_relax_template.sh inputs/processed/3O2H_AB_pep5.pdb
source config/rosetta_paths.sh

PDB="${1:?Need input PDB}"
OUT="runs/template"
mkdir -p "$OUT"

"$ROSETTA_BIN/relax.static.linuxgccrelease" \
  -s "$PDB" \
  -database "$ROSETTA_DB" \
  @config/flags/relax.flags \
  -out:path:all "$OUT"

echo "Relax outputs in: $OUT"
