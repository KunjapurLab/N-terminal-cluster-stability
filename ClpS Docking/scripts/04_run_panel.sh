#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/04_run_panel.sh runs/template/<relaxed_pdb> config/peptides.txt 200
source config/rosetta_paths.sh

TEMPLATE="${1:?Need template PDB (relaxed/trimmed)}"
PEPTIDES_FILE="${2:?Need peptides list file}"
NSTRUCT="${3:-200}"

if [[ ! -f "$TEMPLATE" ]]; then
  echo "Template not found: $TEMPLATE" >&2
  exit 1
fi

while IFS= read -r pep; do
  pep="$(echo "$pep" | tr -d '\r' | xargs)"
  [[ -z "$pep" ]] && continue

  RESFILE="inputs/resfiles/${pep}.resfile"
  if [[ ! -f "$RESFILE" ]]; then
    echo "Missing resfile: $RESFILE (run scripts/02_make_resfiles.py)" >&2
    exit 1
  fi

  OUT="runs/${pep}"
  mkdir -p "$OUT"

  echo "=== ${pep} ==="

  # 1) Mutate peptide with FixBB
  "$ROSETTA_BIN/fixbb.static.linuxgccrelease" \
    -s "$TEMPLATE" \
    -resfile "$RESFILE" \
    -database "$ROSETTA_DB" \
    @config/flags/fixbb.flags \
    -out:path:all "$OUT"

  FIXED="$(ls "$OUT"/*_0001.pdb 2>/dev/null | head -n 1)"
  if [[ -z "${FIXED:-}" ]]; then
    echo "FixBB did not produce *_0001.pdb in $OUT" >&2
    exit 1
  fi

  # 2) Prepack
  "$ROSETTA_BIN/FlexPepDocking.static.linuxgccrelease" \
    -s "$FIXED" \
    -database "$ROSETTA_DB" \
    @config/flags/prepack.flags \
    -out:path:all "$OUT"

  PREPACK="$(ls "$OUT"/*_prepack_0001.pdb 2>/dev/null | head -n 1)"
  if [[ -z "${PREPACK:-}" ]]; then
    echo "Prepack did not produce *_prepack_0001.pdb in $OUT" >&2
    exit 1
  fi

  # 3) Refine
  "$ROSETTA_BIN/FlexPepDocking.static.linuxgccrelease" \
    -s "$PREPACK" \
    -database "$ROSETTA_DB" \
    @config/flags/refine.flags \
    -nstruct "$NSTRUCT" \
    -out:path:all "$OUT"

  echo "Done: ${pep} (outputs in $OUT)"
done < "$PEPTIDES_FILE"
