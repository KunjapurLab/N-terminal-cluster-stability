import os
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import numpy as np

# Codon table for translation
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def matches_consensus(sequence, consensus):
    for seq_base, cons_base in zip(sequence, consensus):
        if cons_base != 'X' and seq_base != cons_base:
            return False
    return True

def find_consensus_matches(read, consensus):
    """Find if the consensus motif appears anywhere in the read."""
    consensus_length = len(consensus)
    for i in range(len(read) - consensus_length + 1):
        substring = read[i:i + consensus_length]
        if matches_consensus(substring, consensus):
            return True
    return False

def translate_sequence(seq):
    """Translates a nucleotide sequence into an amino acid sequence."""
    protein = ''.join(codon_table.get(seq[i:i+3], 'X') for i in range(0, len(seq), 3))
    return protein

def analyze_fasta_upstream(FASTA_file, common, consensus, bases_upstream=15):
    """Analyze per-position base and amino acid distributions upstream of a common sequence."""
    base_position_counts = defaultdict(Counter)  # key = position, value = base counts
    aa_position_counts = defaultdict(Counter)    # key = position, value = AA counts
    match_count = 0
    total_upstream_seqs = []

    with open(FASTA_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            read = line.strip()

            if not find_consensus_matches(read, consensus):
                continue
            match_count += 1

            common_index = read.find(common)
            if common_index >= bases_upstream:
                upstream_seq = read[common_index - bases_upstream:common_index]
                total_upstream_seqs.append(upstream_seq)

                for pos, base in enumerate(upstream_seq):
                    base_position_counts[pos][base] += 1

    if bases_upstream % 3 == 0:
        for seq in total_upstream_seqs:
            aa_seq = translate_sequence(seq)
            if '*' not in aa_seq:
                for pos, aa in enumerate(aa_seq):
                    aa_position_counts[pos][aa] += 1

    return base_position_counts, aa_position_counts, match_count

def process_fasta_files(fasta_files, common, consensus, bases_upstream=15):
    total_base_counts = defaultdict(Counter)
    total_aa_counts = defaultdict(Counter)
    total_matches = 0

    for fasta_file in fasta_files:
        if os.path.exists(fasta_file):
            base_counts, aa_counts, match_count = analyze_fasta_upstream(fasta_file, common, consensus, bases_upstream)
            for pos in base_counts:
                total_base_counts[pos].update(base_counts[pos])
            for pos in aa_counts:
                total_aa_counts[pos].update(aa_counts[pos])
            total_matches += match_count
        else:
            print(f"Warning: {fasta_file} not found!")

    return total_base_counts, total_aa_counts, total_matches

def print_and_save_distributions(base_counts, aa_counts, bases_upstream, output_file):
    with open(output_file, 'w') as f:
        f.write("Base Pair Distributions (per position):\n")
        print("\nBase Pair Distributions (per position):")
        for pos in range(bases_upstream):
            total = sum(base_counts[pos].values())
            if total == 0:
                continue
            f.write(f"Position {pos+1}:\n")
            print(f"Position {pos+1}:")
            for base in 'ACGT':
                percent = (base_counts[pos][base] / total) * 100 if total > 0 else 0
                f.write(f"  {base}: {percent:.2f}%\n")
                print(f"  {base}: {percent:.2f}%")
            f.write("\n")
            print()

        if aa_counts:
            f.write("\nAmino Acid Distributions (per position):\n")
            print("\nAmino Acid Distributions (per position):")
            for pos in range(bases_upstream // 3):
                total = sum(aa_counts[pos].values())
                if total == 0:
                    continue
                f.write(f"AA Position {pos+1}:\n")
                print(f"AA Position {pos+1}:")
                for aa, count in aa_counts[pos].most_common():
                    percent = (count / total) * 100 if total > 0 else 0
                    f.write(f"  {aa}: {percent:.2f}%\n")
                    print(f"  {aa}: {percent:.2f}%")
                f.write("\n")
                print()


def plot_stacked_bar(data_counts, labels, title, ylabel, output_file, base_plot=False):
    """Plot stacked bar chart for base or AA distributions."""
    import numpy as np

    positions = sorted(data_counts.keys())

    if base_plot:
        # For base pairs: only A, C, G, T
        categories = ['A', 'C', 'G', 'T']
    else:
        # For amino acids: all observed AAs
        categories = sorted({cat for counts in data_counts.values() for cat in counts})

    # Build a matrix for plotting
    data_matrix = []
    for cat in categories:
        row = [(data_counts[pos][cat] if cat in data_counts[pos] else 0) for pos in positions]
        data_matrix.append(row)

    # Convert counts to percentages
    data_matrix = np.array(data_matrix)
    column_sums = np.sum(data_matrix, axis=0)
    percent_matrix = (data_matrix / column_sums) * 100

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    bottom = np.zeros(len(positions))

    for i, cat in enumerate(categories):
        bars = ax.bar(positions, percent_matrix[i], bottom=bottom, label=cat)
        if base_plot:  # Only annotate inside bars if base_plot
            for bar, percent in zip(bars, percent_matrix[i]):
                if percent > 2:  # Only label if big enough
                    height = bar.get_height()
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bottom[bars.index(bar)] + height / 2,
                        f"{percent:.1f}%",
                        ha='center', va='center', fontsize=8, color='black'
                    )
        bottom += percent_matrix[i]

    ax.set_xlabel('Position', fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.set_xticks(positions)
    ax.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


# Main
if __name__ == "__main__":
    #If analyzing less than four files, placeholder names can be listed in "empty" slots
    fasta_files = ["21Jan2025 B1.fasta", "21Jan2025 B2.fasta", "21Jan2025 B3.fasta", "21Jan2025 B4.fasta"]
    common = "ATCAGTGACTTC"
    consensus = "XXXACTGGAGGATGGTCGCACGCTTTCGGACTACAACATCCAGAAAGAATCTACCCTTCATTTGGTTCTGCGTCTGCGTGGAGGAXXXXXXXXXXXXXXXATCAGTGACTTCATCGCATCCAAGGGCGAGGAGCTCTTTACTGGCGTAGTACCAATT"

    bases_upstream = 15  # Change to whatever you want
    output_txt = "upstream_distributions.txt"
    base_plot = "basepair_distribution.png"
    aa_plot = "aminoacid_distribution.png"

    base_counts, aa_counts, total_matches = process_fasta_files(fasta_files, common, consensus, bases_upstream)

    print(f"\nTotal sequences matching the consensus: {total_matches}\n")
    print_and_save_distributions(base_counts, aa_counts, bases_upstream, output_txt)
    print(f"\nDistributions saved to '{output_txt}'.")

    # Plot base pairs (A, C, G, T only, with % labels inside bars)
    plot_stacked_bar(base_counts, labels=['A', 'C', 'G', 'T'], 
                 title="Base Pair Distribution per Position",
                 ylabel="Percentage (%)",
                 output_file=base_plot,
                 base_plot=True)

# Plot amino acids (all observed AAs, no inside labels)
    if aa_counts:
        plot_stacked_bar(aa_counts, labels=None, 
                     title="Amino Acid Distribution per Position",
                     ylabel="Percentage (%)",
                     output_file=aa_plot,
                     base_plot=False)