-----NGS Software Suite Read Me------



Dear curious scientist,

Thank you for considering our bioinformatics software for your studies. This software suite is designed to take you from raw FASTQ.gz/FASTA data through to data visualization of multiple tiered bins of FACS data. This suite is introduced in "Developing a flow cytometric method to evaluate the stability of protein N-termini" by Sen et al., Methods In Enzymology(2025). As such, we suggest you use that text to accompany this read me. We also utilize protein stability index (PSI) as our method for evaluating the stability of sequences. We direct you to work from the Elledge Lab, particularly Yen et al. (2008) and Timms et al. (2019) for strong references regarding PSI. Furthermore, we have included a small test dataset to help you pilot out this workflow.

This code was developed using Python 3.13 and was last updated in May 2025. A working version of Python as well as an IDE for simple editing are required. All instructions are written for usage in a Windows operating system. Instructions below are written to match the sequences and aims set in Sen et al. 2025.

If you need assistance, please contact kunjapurlab@udel.edu and put [N-degron methods] at the beginning of your email subject line.

-----Set up-----

1.) Download and extract all scripts into the same folder. Add any .FASTQ/.FASTA files you would like to analyze to this folder.
2.) You will have to install a series of python libraries. You can optionally choose to do this in a new virtual environment. Type CMD+R and then type:
	pip install [placeholder]
3.) If necessary, begin by converting your FASTQ or FASTQ.gz files to a FASTA format using the "FASTQ.gz to FASTA.py" or "FASTQ to FASTA.py" scripts. Note that for large 
compressed file sizes you can expect a 3-5x increase in file size, so have the appropriate hard-drive space ready before hand.
4.) If necessary, demultiplex samples pooled in the same run using "FASTA Demultiplexer.py"
5.) To analyze the composition of your mutagenized region to evaluate base pair and amino acid distributions, optionally run the "FASTA Codon and AA Distributions.py" script.
6.) Next, run "Sequence Dictionary Generator.py" to generate a sequence-stability database. This will be the foundation for future analysis.
7.) From this database, feel free to utilize any visualization scripts (demarked with 4.) in the folder). Adjust the line in the code to account for your database.

-----FAQs-----

Q: What kind of data is best used for these scripts?

A: This code has been developed to analyze 15 consecutive mutagenized base pairs across four sorted bins. The code would have to be adapted to alternative data sets.