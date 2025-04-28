import sqlite3
import time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
    print(f"Database load time: {end_time - start_time:.2f} seconds")
    return result_list

def filter_sequences_by_read_count(result_list):
    """
    Filter sequences with total read counts greater than 10 across bins.
    """
    start_time = time.time()
    filtered_result_list = [entry for entry in result_list if sum(entry[1:5]) > 20]
    end_time = time.time()
    print(f"Filtering time: {end_time - start_time:.2f} seconds")
    return filtered_result_list


def calculate_mean_psi_all_positions(result_list):
    """
    Calculate the mean PSI for each amino acid at each position, excluding sequences containing 'X'.
    """
    start_time = time.time()
    psi_data = []

    for entry in result_list:
        sequence = entry[0]
        psi = entry[5]

        # Exclude sequences containing 'X'
        if 'X' in sequence:
            continue

        # Collect PSI data for all positions in the sequence
        for position, amino_acid in enumerate(sequence, start=1):
            psi_data.append({'Position': position, 'AminoAcid': amino_acid, 'PSI': psi})

    # Create a DataFrame
    df = pd.DataFrame(psi_data)

    # Group by amino acid and position, then calculate the mean PSI
    mean_psi = df.groupby(['AminoAcid', 'Position'])['PSI'].mean().unstack()
    end_time = time.time()
    print(f"Mean PSI calculation time: {end_time - start_time:.2f} seconds")
    return mean_psi

def plot_aa_position_heatmap(mean_psi):
    """
    Plot a heatmap for mean PSI values grouped by amino acids and positions.
    """
    start_time = time.time()
    plt.figure(figsize=(14, 12))
    sns.heatmap(mean_psi, annot=True, fmt=".2f", cmap="Blues", cbar=True, vmin=1.5, vmax=4.0)
    plt.title("WT BL21 Mean PSI Heatmap (Amino Acids vs. Positions)")
    plt.xlabel("Position")
    plt.ylabel("Amino Acid")
    plt.tight_layout()
    plt.show()
    end_time = time.time()
    print(f"Heatmap plotting time: {end_time - start_time:.2f} seconds")
    
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
            psi = (4*b4 + 3 * b3 + 2 * b2 + 1 *b1) / total_reads
        else:
            psi = 0  # Default to 0 if there are no reads

        # Round PSI to two decimal places
        recalculated_entry = (sequence, b1, b2, b3, b4, round(psi, 2))
        recalculated_results.append(recalculated_entry)

    return recalculated_results


# Main part of the code
if __name__ == '__main__':
    db_file = "Feb2025 NGS for Dec2024 WT resort.db"
    output_file = "LFTR- sequences_with_x.txt"

    # Load data from the database
    result_list = load_data_from_db(db_file)

    # Recalculate PSI values
    recalculated_result_list = recalculate_psi(result_list)

    # Filter sequences
    filtered_result_list = filter_sequences_by_read_count(recalculated_result_list)

    # Calculate the mean PSI for all positions, excluding 'X'
    mean_psi_all_positions = calculate_mean_psi_all_positions(filtered_result_list)

    # Plot the heatmap
    plot_aa_position_heatmap(mean_psi_all_positions)

