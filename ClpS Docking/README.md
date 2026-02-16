# FlexPepDock peptide panel pipeline

Scripts to run a peptide-variant panel through Rosetta:
clean → trim peptide → generate resfiles → relax template → FixBB → FlexPepDock prepack → FlexPepDock refine → summarize I_sc.

## Requirements
- Rosetta installed locally/HPC (this repo does not include Rosetta)
- Python 3
- Bash

## Setup
```bash
cd "ClpS Docking"
cp config/rosetta_paths.example.sh config/rosetta_paths.sh
# edit config/rosetta_paths.sh for your installation
```

Put your PDB in inputs/ (e.g. inputs/3O2H.pdb), then:

```bash
bash scripts/00_clean_and_concat.sh inputs/3O2H.pdb A B
python3 scripts/01_trim_peptide.py --in inputs/processed/3O2H_AB.pdb --out inputs/processed/3O2H_AB_pep5.pdb --pep-chain B --start 1 --end 5
python3 scripts/02_make_resfiles.py --peptides config/peptides.txt --outdir inputs/resfiles --chain B
bash scripts/03_relax_template.sh inputs/processed/3O2H_AB_pep5.pdb
bash scripts/04_run_panel.sh runs/template/3O2H_AB_pep5_0001.pdb config/peptides.txt 200
bash scripts/05_summarize_isc.sh runs 10
bash scripts/06_topk_models.sh runs 5
```

Chain letters and peptide residue numbering can vary by PDB. If your peptide chain is not B or residues aren’t numbered from 1, adjust:

01_trim_peptide.py --pep-chain ... --start ... --end ...

02_make_resfiles.py --chain ...

Outputs go under runs/<PEPTIDE>/.


---

## `config/rosetta_paths.example.sh`
```bash
#!/usr/bin/env bash
# Copy to config/rosetta_paths.sh and edit for your machine/HPC.

export ROSETTA_ROOT="/path/to/rosetta"
export ROSETTA_BIN="$ROSETTA_ROOT/main/source/bin"
export ROSETTA_DB="$ROSETTA_ROOT/main/database"
export ROSETTA_PY="$ROSETTA_ROOT/main/source/scripts/python"
```
