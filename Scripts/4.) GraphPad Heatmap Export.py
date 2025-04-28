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
    Filter sequences with total read counts greater than 20 across bins.
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

def output_data_for_graphpad(mean_psi, output_file):
    """
    Output the mean PSI data in CSV format suitable for GraphPad Prism.
    The output has an 'AminoAcid' column and one column per position.
    """
    # Reset index so that the amino acids become a column in the CSV file.
    df_to_output = mean_psi.reset_index()
    df_to_output.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")

def plot_aa_position_heatmap(mean_psi):
    """
    Plot a heatmap for mean PSI values grouped by amino acids and positions.
    """
    start_time = time.time()
    plt.figure(figsize=(14, 12))
    sns.heatmap(mean_psi, annot=True, fmt=".2f", cmap="Blues", cbar=True)
    plt.title("ClpS- Mean PSI Heatmap (Amino Acids vs. Positions)")
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
            psi = (4*b4 + 3*b3 + 2*b2 + 1*b1) / total_reads
        else:
            psi = 0  # Default to 0 if there are no reads

        # Round PSI to two decimal places
        recalculated_entry = (sequence, b1, b2, b3, b4, round(psi, 2))
        recalculated_results.append(recalculated_entry)

    return recalculated_results

# Main part of the code
if __name__ == '__main__':
    db_file = "Feb2025 NGS for Dec2024 WT resort.db"
    graphpad_output_file = "WTBL21Feb2025AVITI_graphpad_heatmap_data.csv"

    # Load data from the database
    result_list = load_data_from_db(db_file)

    # Uncomment the next line if you wish to recalculate PSI values using your weighted formula
    # result_list = recalculate_psi(result_list)

    # Filter sequences based on total read counts
    filtered_result_list = filter_sequences_by_read_count(result_list)

    # Calculate the mean PSI for all positions (excluding sequences with 'X')
    mean_psi_all_positions = calculate_mean_psi_all_positions(filtered_result_list)

    # Output the data to a CSV file formatted for GraphPad Prism
    output_data_for_graphpad(mean_psi_all_positions, graphpad_output_file)

    # Plot the heatmap (optional)
    plot_aa_position_heatmap(mean_psi_all_positions)
