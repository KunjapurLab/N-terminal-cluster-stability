# NGS Software Suite Read Me


## Introduction
Dear curious scientist,

Thank you for considering our bioinformatics software for your studies. This software suite is designed to take you from raw FASTQ.gz/FASTA data through to data visualization of multiple tiered bins of FACS data. This suite is introduced in "Developing a flow cytometric method to evaluate the stability of protein N-termini" by Sen et al., Methods In Enzymology (2025). As such, we suggest you use that text to accompany this read me. We utilize protein stability index (PSI) as our metric for evaluating the stability of sequences. We direct you to work from the Elledge Lab, particularly Yen et al. (2008) and Timms et al. (2019) for strong references regarding the development and utilization of this weighted average metric. Furthermore, we have included a small test dataset to help you pilot out this workflow.


## Features
- Convert raw FASTQ(.gz) → FASTA  
- Demultiplex pooled FASTA files  
- Compute Protein Stability Index (PSI)  
- Generate amino-acid & codon usage reports  
- Visualize multi-tiered FACS binned data as heatmaps, WebLogos, and degron lists  
- …and more!

## Installation
1. Prerequisites
   - Windows 10 or later  
   - Python 3.8–3.13  
   - Sufficient disk space for your expanded FASTA files

2. Setup  
   ```bash
   # (optional) create a virtualenv:
   python -m venv ngs-env
   source ngs-env/Scripts/activate

   # install dependencies
   pip install -r requirements.txt

## Set up Instructions
1.) Download and extract all scripts into the same folder. Add any .FASTQ/.FASTA files you would like to analyze to this folder.

2.) You will have to install a series of python libraries. You can optionally choose to do this in a new virtual environment. Within your target environment, install the libraries using "pip install -r requirements.txt"

3.) If necessary, begin by converting your FASTQ or FASTQ.gz files to a FASTA format using the "FASTQ.gz to FASTA.py" or "FASTQ to FASTA.py" scripts. Note that for large 
compressed file sizes you can expect a 3-5x increase in file size, so have the appropriate hard-drive space ready before hand.

4.) If necessary, demultiplex samples pooled in the same run using "FASTA Demultiplexer.py"

5.) Optionally, to analyze amino acid and codon distribution and bias, run "FASTA Codon and AA Distributions.py" 

6.) Next, run "Sequence Dictionary Generator.py" to generate a sequence-stability database. This will be the foundation for future analysis.

7.) From this database, feel free to utilize any visualization scripts (demarked with 4.) in the folder). Adjust the line in the code to account for your database.

## Citations
If you utilize these scripts, please cite:
Sen et al. Biorxiv (2025) "Combinatorial mutagenesis of N-terminal sequences in E. coli reveals unexpected N-degron stability determinants." 
Sen et al., Methods In Enzymology (2025). "Developing a flow cytometric method to evaluate the stability of protein N-termini."

At the moment, these publications are in the preparation & submission phase.

## Contact
If you need assistance, please contact kunjapurlab [at] udel (dot) edu and put [N-degron methods] at the beginning of your email subject line. Learn more about the lab here: https://www.kunjapurlab.org/ 

You may also email the first author and primary developer of this code at sabyasachi.sen96@gmail.com. Feel free to reach out on LinkedIn as well: https://www.linkedin.com/in/sabyasachi-sen/ 

## FAQS 
Q: What versions of Python has this code been tested?
A: The code has been developed and tested on Python 3.13

## History
This code was last updated in April 2025
