import os
import sqlite3
from collections import Counter

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

# Add PSI calculations to result list
def PSI_calc(result_list):
    for entry in result_list:
        total_reads = sum(entry[1:5])
        if total_reads > 0:
            PSI = (4 * entry[4] + 3 * entry[3] + 2 * entry[2] + 1 * entry[1]) / total_reads
            entry[5] = round(PSI, 2)
    return result_list

# Save result list to SQLite database
def save_result_list_to_db(result_list, db_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ResultList (
            Sequence TEXT PRIMARY KEY,
            B1 INTEGER,
            B2 INTEGER,
            B3 INTEGER,
            B4 INTEGER,
            PSI REAL
        )
    ''')
    cursor.executemany('''
        INSERT OR REPLACE INTO ResultList (Sequence, B1, B2, B3, B4, PSI)
        VALUES (?, ?, ?, ?, ?, ?)
    ''', result_list)
    conn.commit()
    conn.close()
    print(f"Result list saved to {db_file}.")

# Function to check if a sequence matches a consensus sequence
def matches_consensus(sequence, consensus):
    for seq_base, cons_base in zip(sequence, consensus):
        if cons_base != 'X' and seq_base != cons_base:
            return False
    return True

# Function to find matches to the consensus using a sliding window
def find_consensus_matches(read, consensus):
    """Finds if the consensus motif appears anywhere in the read."""
    consensus_length = len(consensus)
    for i in range(len(read) - consensus_length + 1):
        substring = read[i:i + consensus_length]
        if matches_consensus(substring, consensus):
            return True  # A match is found
    return False

# Process a single FASTA file to generate sequence data and check against consensus
def fasta_to_deg_dict(FASTA_file, common, consensus):
    """Processes a FASTA file to count amino acid sequences and consensus matches."""
    temp_codons_list = []
    AA_list = []
    match_count = 0  # Track how many reads match the consensus

    with open(FASTA_file, 'r') as f:
        for line in f:
            if line[0] != ">":
                read = line.strip()  # Read the entire sequence

                # Check for consensus match anywhere in the read
                if find_consensus_matches(read, consensus):
                    match_count += 1  # Increment match count
                else:
                    continue  # Skip reads that don't match the consensus

                # Find the `common` sequence and extract 15 bases before it
                common_index = read.find(common)
                if common_index >= 15:  # Ensure there are enough bases in front of the `common` sequence
                    preceding_sequence = read[common_index - 15:common_index]
                    temp_codons_list.append(preceding_sequence)

    # Translate extracted 15-base sequences to amino acids
    for entry in temp_codons_list:
        protein = ''.join(codon_table.get(entry[i:i + 3], "X") for i in range(0, 15, 3))
        if '*' not in protein:  # Skip sequences with stop codons
            AA_list.append(protein)

    # Count the occurrences of each amino acid sequence
    AA_count = Counter(AA_list)
    return dict(AA_count), match_count



# Process multiple FASTA files
def process_fasta_files(fasta_files, common, consensus):
    """Processes multiple FASTA files and counts consensus matches."""
    result_dict = {}
    total_match_count = 0  # Track total matches across all files

    for num, fasta_file in enumerate(fasta_files):
        if os.path.exists(fasta_file):
            AA_count, match_count = fasta_to_deg_dict(fasta_file, common, consensus)
            total_match_count += match_count  # Accumulate match counts
            for sequence, count in AA_count.items():
                if sequence not in result_dict:
                    result_dict[sequence] = [sequence] + [0] * 4 + [0.0]
                result_dict[sequence][num + 1] = count

    # Convert dict to list for PSI calculation
    result_list = list(result_dict.values())
    return PSI_calc(result_list), total_match_count


# Main execution
if __name__ == '__main__':
    fasta_files = ["10M ClpS- B1.fasta", "10M ClpS- B2.fasta", "10M ClpS- B3.fasta", "10M ClpS- B4.fasta"]
    common = "ATCAGTGACTTC"
    # Updated consensus sequence with X as wildcard
    consensus = "XXXACTGGAGGATGGTCGCACGCTTTCGGACTACAACATCCAGAAAGAATCTACCCTTCATTTGGTTCTGCGTCTGCGTGGAGGAXXXXXXXXXXXXXXXATCAGTGACTTCATCGCATCCAAGGGCGAGGAGCTCTTTACTGGCGTAGTACCAATT"
    db_file = "10M ClpS- consensus 24Nov.db"

    # Process the FASTA files
    result_list, total_match_count = process_fasta_files(fasta_files, common, consensus)

    # Save the results to the database
    save_result_list_to_db(result_list, db_file)

    # Print matching statistics
    total_sequences = sum(sum(entry[1:5]) for entry in result_list)
    print(f"Total sequences: {total_sequences}")
    print(f"Sequences matching the consensus: {total_match_count}")
    print(f"Difference: {total_sequences - total_match_count}")

