import os
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

def calculate_mean_psi_two_positions(result_list, pos1, pos2):
    """
    Calculate the mean PSI for each combination of two positions.

    Parameters:
        result_list (list): The filtered result list.
        pos1 (int): The first position.
        pos2 (int): The second position.

    Returns:
        pd.DataFrame: Mean PSI grouped by the two positions.
    """
    start_time = time.time()
    psi_data = []

    for entry in result_list:
        sequence = entry[0]
        psi = entry[5]

        # Exclude sequences containing 'X'
        if 'X' in sequence or len(sequence) < max(pos1, pos2):
            continue

        # Extract the amino acids at the specified positions
        p1 = sequence[pos1 - 1]
        p2 = sequence[pos2 - 1]
        psi_data.append({f'P{pos1}': p1, f'P{pos2}': p2, 'PSI': psi})

    # Create a DataFrame
    df = pd.DataFrame(psi_data)

    # Group by the two positions, then calculate the mean PSI
    mean_psi = df.groupby([f'P{pos1}', f'P{pos2}'])['PSI'].mean().unstack()
    end_time = time.time()
    print(f"Mean PSI calculation for P{pos1}-P{pos2} time: {end_time - start_time:.2f} seconds")
    return mean_psi

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

        # Round PSI to two decimal places
        recalculated_entry = (sequence, b1, b2, b3, b4, round(psi, 2))
        recalculated_results.append(recalculated_entry)

    return recalculated_results

def plot_heatmap(mean_psi, pos1, pos2, save_dir="heatmaps"):
    """
    Plot and save a heatmap for mean PSI values grouped by two specified positions.

    Parameters:
        mean_psi (pd.DataFrame): Data for the heatmap.
        pos1 (int): The first position.
        pos2 (int): The second position.
        save_dir (str): Directory where heatmaps will be saved.
    """
    start_time = time.time()
    # Ensure the save directory exists
    os.makedirs(save_dir, exist_ok=True)
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(mean_psi, annot=True, fmt=".2f", cmap="Blues", cbar=True, vmin=1.5, vmax=3.5)
    plt.title(f"WT BL21 Mean PSI Heatmap (P{pos1}-P{pos2} Amino Acid Combinations)")
    plt.xlabel(f"P{pos2}")
    plt.ylabel(f"P{pos1}")
    plt.tight_layout()
    
    # Save the heatmap to the specified folder
    file_name = os.path.join(save_dir, f"heatmap_P{pos1}_P{pos2}.png")
    plt.savefig(file_name, dpi=300)
    plt.close()  # Close the plot to prevent it from popping up
    end_time = time.time()
    print(f"Heatmap saved for P{pos1}-P{pos2}: {file_name}. Plotting time: {end_time - start_time:.2f} seconds")

# Main part of the code
if __name__ == '__main__':
    db_file = "Feb2025 NGS for Dec2024 WT resort.db"

    # Load data from the database
    result_list = load_data_from_db(db_file)

    # Optionally recalculate PSI values
    # recalculated_result_list = recalculate_psi(result_list)

    # Filter sequences
    filtered_result_list = filter_sequences_by_read_count(result_list)

    # Generate and save heatmaps for all unique two-position combinations
    positions = [1,2, 3, 4, 5]
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            pos1, pos2 = positions[i], positions[j]
            mean_psi = calculate_mean_psi_two_positions(filtered_result_list, pos1, pos2)
            plot_heatmap(mean_psi, pos1, pos2)
