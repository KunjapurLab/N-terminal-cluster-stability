#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/00_clean_and_concat.sh inputs/3O2H.pdb A B
source config/rosetta_paths.sh

PDB="${1:?Need input PDB path}"
CHAIN_A="${2:?Need receptor chain (e.g., A)}"
CHAIN_B="${3:?Need peptide chain (e.g., B)}"

NAME="$(basename "$PDB" .pdb)"
OUT="inputs/processed"
mkdir -p "$OUT"

python "$ROSETTA_PY/clean_pdb.py" "$PDB" "$CHAIN_A"
python "$ROSETTA_PY/clean_pdb.py" "$PDB" "$CHAIN_B"

cat "${NAME}_${CHAIN_A}.pdb" "${NAME}_${CHAIN_B}.pdb" > "$OUT/${NAME}_${CHAIN_A}${CHAIN_B}.pdb"

echo "Wrote: $OUT/${NAME}_${CHAIN_A}${CHAIN_B}.pdb"
