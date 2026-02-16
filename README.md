# N-terminal Cluster Stability Toolkit

Dear curious scientist,

Thank you for considering our bioinformatics software. This software suite is designed to take you from raw FASTQ.gz/FASTA data through to data visualization of multiple tiered bins of FACS data. This suite is introduced in "Combinatorial mutagenesis of N-terminal sequences reveals unexpected and expanded stability determinants of the Escherichia coli N-degron pathway" by Sen et al., Biorxiv (2025) and "Developing a flow cytometric method to evaluate the stability of protein N-termini" by Sen et al., accepted at Methods In Enzymology (2025). As such, we suggest you use that text to accompany this read me. We utilize protein stability index (PSI) as our metric for evaluating the stability of sequences. We direct you to work from the Elledge Lab, particularly Yen et al. (2008) and Timms et al. (2019) for strong references regarding the development and utilization of this weighted average metric. Furthermore, we have included a small test dataset to help you pilot out this workflow.

## Repository Layout

- `Scripts/`: FASTQ/FASTA processing and database-driven visualization utilities
- `N-FIVE/`: PSI prediction GUI scripts and SHAP analysis scripts
- `Databases/`: example SQLite databases used by analysis scripts
- `Sample Data/`: sample sequencing files

## Quick Start

1. Create and activate a virtual environment (Windows).
2. Install dependencies:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

If you use Command Prompt instead of PowerShell:

```bat
python -m venv .venv
.venv\Scripts\activate.bat
pip install -r requirements.txt
```

3. Run scripts from the repository root (recommended):

```powershell
python "N-FIVE/SS 13March2025 Predict PSI with model v7.py"
python "N-FIVE/SS 21Feb2025 SHAP Interaction Visualizer.py"
python "Scripts/4.) 5 Position Heatmap Generator.py"
```

## Set up Instructions
1.) Clone/download the full repository and run scripts from the repository root. Put your input `.FASTQ/.FASTQ.gz/.FASTA` files where convenient (for example in `Sample Data/`).

2.) You will have to install a series of python libraries. You can optionally choose to do this in a new virtual environment. Within your target environment, install the libraries using "pip install -r requirements.txt"

3.) If necessary, begin by converting your FASTQ or FASTQ.gz files to a FASTA format using the "FASTQ.gz to FASTA.py" or "FASTQ to FASTA.py" scripts. Note that for large compressed file sizes you can expect a 3-5x increase in file size, so have the appropriate hard-drive space ready before hand.

4.) If necessary, demultiplex samples pooled in the same run using "FASTA Demultiplexer.py"

5.) Optionally, to analyze amino acid and codon distribution and bias, run "FASTA Codon and AA Distributions.py"

6.) Next, run "Sequence Dictionary Generator.py" to generate a sequence-stability database. By default, output databases are written to `Databases/`.

7.) From this database, utilize any visualization scripts (demarked with `4.)`). These scripts now auto-locate databases in `Databases/`, and most also accept a database path as the first command-line argument.

## ClpS Docking

If you are interested in running our ClpS docking workflows, please refer to the subfolder for ClpS Docking and its associated ReadMe.

## File Compatibility

The N-FIVE and analysis scripts auto-locate required files in standard repo locations.

- Model files:
  - `N-FIVE/N-FIVE WT model.pkl` (preferred)
  - scripts look for database files in `Databases/` automatically
- SHAP cache:
  - stored as `N-FIVE/shap_interaction_values.pkl`

Most scripts also accept a database path from the command line as the first argument.

Example:

```powershell
python "Scripts/4.) GraphPad Heatmap Export.py" "Databases/Feb2025 NGS for Dec2024 WT resort.db"
```

For the logo script (three databases):

```powershell
python "Scripts/4.) Least Stable Sequences Logomaker.py" "Databases/Feb2025 NGS for Dec2024 WT resort.db" "Databases/10M LFTR- consensus 22Nov.db" "Databases/10M ClpS- consensus 24Nov.db"
```

## Citation

If you utilize these scripts, please cite: Sen et al. Biorxiv (2025). "Combinatorial mutagenesis of N-terminal sequences reveals unexpected and expanded stability determinants of the Escherichia coli N-degron pathway." 

## Contact

For support, contact the Kunjapur Lab: `kunjapurlab [at] udel (dot) edu` with subject prefix `[N-degron methods]`.
