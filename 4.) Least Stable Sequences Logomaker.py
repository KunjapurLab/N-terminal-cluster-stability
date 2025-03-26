import sqlite3
import time
import pandas as pd
import matplotlib.pyplot as plt
import logomaker  # pip install logomaker

def load_data_from_db(db_file):
    """
    Load data from SQLite database.
    """
    start_time = time.time()
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT Sequence, B1, B2, B3, B4, PSI FROM ResultList")
    result_list = cursor.fetchall()
    conn.close()
    end_time = time.time()
    print(f"[{db_file}] Database load time: {end_time - start_time:.2f} seconds")
    return result_list

def recalculate_psi(result_list):
    """
    Recalculate PSI values for each entry in the result_list.
    """
    recalculated_results = []
    for entry in result_list:
        sequence = entry[0]
        b1, b2, b3, b4 = entry[1:5]
        total_reads = b1 + b2 + b3 + b4
        if total_reads > 0:
            # Recalculate PSI using the weighted sum formula
            psi = (b4 + 2 * b3 + 3 * b2 + 4 * b1) / total_reads
        else:
            psi = 0  # Default to 0 if there are no reads
        recalculated_results.append((sequence, b1, b2, b3, b4, round(psi, 2)))
    return recalculated_results

def filter_sequences_by_read_count(result_list, threshold=4):
    """
    Filter sequences with total read counts greater than a threshold.
    """
    start_time = time.time()
    filtered_result_list = [entry for entry in result_list if sum(entry[1:5]) > threshold]
    end_time = time.time()
    print(f"Filtering time: {end_time - start_time:.2f} seconds")
    return filtered_result_list

def process_db_full_sequence(db_file, recalc=False, threshold=4):
    """
    Process a database and return full entries (including the entire motif and PSI)
    after loading, optionally recalculating PSI, and filtering by read count.
    """
    result_list = load_data_from_db(db_file)
    if recalc:
        result_list = recalculate_psi(result_list)
    result_list = filter_sequences_by_read_count(result_list, threshold=threshold)
    return result_list

def generate_weblogo(sequences, title="WebLogo of Sequences", output_file=None):
    """
    Generate a sequence logo from a list of motif sequences using Logomaker.
    The logo will display the deviation from the baseline occurrence rate (1/20)
    for each amino acid at each motif position. Amino acids with frequencies above
    this mean appear above the y-axis and those below appear below.
    
    Custom colors are applied:
      - F, L, W, Y: blue
      - R, K: orange
      - All others: black
    """
    # Ensure all sequences have the same length.
    seq_lengths = {len(seq) for seq in sequences}
    if len(seq_lengths) != 1:
        raise ValueError("Not all sequences are the same length. Cannot generate a proper weblogo.")
    motif_length = seq_lengths.pop()
    
    # Define the 20 canonical amino acids in a fixed order.
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    
    # Build a count matrix: each row represents a motif position.
    counts_list = [{} for _ in range(motif_length)]
    for seq in sequences:
        for pos, aa in enumerate(seq):
            counts_list[pos][aa] = counts_list[pos].get(aa, 0) + 1

    # Convert the list of dictionaries into a DataFrame.
    count_df = pd.DataFrame(counts_list)
    # Ensure all canonical amino acids are included (fill missing ones with 0).
    count_df = count_df.reindex(columns=amino_acids, fill_value=0)
    count_df.index.name = "Position"
    
    # Convert counts to frequencies (per position).
    freq_df = count_df.div(count_df.sum(axis=1), axis=0)
    
    # Compute deviation from the baseline frequency.
    # Baseline is assumed to be 1/20 for the 20 canonical amino acids.
    baseline = 1 / len(amino_acids)
    diff_df = freq_df - baseline

    # Define custom color scheme.
    custom_color_scheme = {}
    for aa in amino_acids:
        if aa in ['F', 'L', 'W', 'Y']:
            custom_color_scheme[aa] = 'blue'
        elif aa in ['R', 'K']:
            custom_color_scheme[aa] = 'orange'
        else:
            custom_color_scheme[aa] = 'black'
    
    # Create the logo using Logomaker.
    plt.figure(figsize=(motif_length / 1.5, 4))
    logo = logomaker.Logo(diff_df, shade_below=0.5, fade_below=0.5,
                           font_name='Arial Rounded MT Bold',
                           color_scheme=custom_color_scheme)
    # Draw a horizontal line at 0 to indicate the baseline.
    plt.axhline(0, color='black', linewidth=1)
    plt.title(title)
    plt.xlabel("Position in Motif")
    plt.ylabel("Deviation from Mean Occurrence Rate")
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file)
        print(f"WebLogo saved to {output_file}")
    else:
        plt.show()

# --- Main Execution ---

if __name__ == '__main__':
    # Specify your database files.
    db_file1 = "Feb2025 NGS for Dec2024 WT resort.db"   # Use original PSI values
    db_file2 = "10M LFTR- consensus 22Nov.db"            # Recalculate PSI
    db_file3 = "10M ClpS- consensus 24Nov.db"            # Recalculate PSI
    
    # Gather sequences from all databases.
    all_sequences = []
    for db_file, recalc in [(db_file1, False), (db_file2, True), (db_file3, True)]:
        entries = process_db_full_sequence(db_file, recalc=recalc, threshold=4)
        all_sequences.extend(entries)
    
    # --- Analysis for the Lowest 10,000 PSI Sequences ---
    # Sort all entries by PSI (ascending) and select the lowest 10,000 sequences.
    all_sequences_sorted = sorted(all_sequences, key=lambda x: x[5])
    lowest_entries = all_sequences_sorted[:100000]
    lowest_sequences = [entry[0] for entry in lowest_entries if entry[0]]
    
    print(f"Selected {len(lowest_sequences)} sequences with the lowest PSI for logo generation.")
    
    try:
        generate_weblogo(lowest_sequences, title="WebLogo of 10,000 Lowest PSI Sequences")
    except ValueError as e:
        print(f"Error generating WebLogo for lowest sequences: {e}")
    
    # --- Analysis for the Highest 100,000 PSI Sequences ---
    # Sort all entries by PSI (descending) and select the highest 100,000 sequences.
    all_sequences_sorted_desc = sorted(all_sequences, key=lambda x: x[5], reverse=True)
    highest_entries = all_sequences_sorted_desc[:100000]
    highest_sequences = [entry[0] for entry in highest_entries if entry[0]]
    
    print(f"Selected {len(highest_sequences)} sequences with the highest PSI for logo generation.")
    
    try:
        generate_weblogo(highest_sequences, title="WebLogo of 100,000 Highest PSI Sequences")
    except ValueError as e:
        print(f"Error generating WebLogo for highest sequences: {e}")
