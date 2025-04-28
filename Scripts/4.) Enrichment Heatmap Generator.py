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

def filter_by_psi(result_list, psi_threshold, comparison):
    """
    Filter sequences based on PSI threshold.
    """
    if comparison == "greater":
        return [entry for entry in result_list if entry[5] > psi_threshold]
    elif comparison == "less":
        return [entry for entry in result_list if entry[5] < psi_threshold]
    else:
        raise ValueError("Comparison must be 'greater' or 'less'.")

def count_amino_acids_at_positions(result_list, max_positions=5):
    """
    Count occurrences of each amino acid at each position.
    """
    counts = [{} for _ in range(max_positions)]
    for entry in result_list:
        sequence = entry[0]
        for i in range(min(len(sequence), max_positions)):
            aa = sequence[i]
            counts[i][aa] = counts[i].get(aa, 0) + 1
    return counts

def calculate_normalized_counts(counts):
    """
    Normalize counts at each position.
    """
    normalized_counts = []
    for position_counts in counts:
        total = sum(position_counts.values())
        normalized_counts.append({aa: count / total for aa, count in position_counts.items()})
    return normalized_counts

def calculate_enrichment(norm_counts_pop1, norm_counts_pop2):
    """
    Calculate enrichment/depletion between two populations.
    """
    enrichment = []
    for pos1, pos2 in zip(norm_counts_pop1, norm_counts_pop2):
        enrichment_pos = {}
        all_aas = set(pos1.keys()).union(set(pos2.keys()))
        for aa in all_aas:
            freq1 = pos1.get(aa, 0)
            freq2 = pos2.get(aa, 0)
            enrichment_pos[aa] = freq1 / freq2 if freq2 > 0 else float('inf') if freq1 > 0 else 0
        enrichment.append(enrichment_pos)
    return enrichment

def plot_enrichment_heatmap(enrichment, max_positions=5):
    """
    Plot a heatmap for enrichment/depletion values.
    """
    for i, enrichment_pos in enumerate(enrichment[:max_positions]):
        df = pd.DataFrame.from_dict(enrichment_pos, orient="index", columns=["Enrichment"]).reset_index()
        df.columns = ["Amino Acid", "Enrichment"]
        df_pivot = df.pivot_table(index="Amino Acid", values="Enrichment")

        plt.figure(figsize=(10, 8))
        sns.heatmap(df_pivot, annot=True, fmt=".2f", cmap="coolwarm", cbar=True)
        plt.title(f"Enrichment Heatmap for Position {i + 1}")
        plt.xlabel("Amino Acids")
        plt.ylabel("Enrichment/Depletion")
        plt.show()
        
def plot_combined_heatmap(enrichment, max_positions=5):
    """
    Plot a combined heatmap for all positions in a single plot.
    Each column corresponds to a position, and rows represent amino acids.
    """
    # Combine enrichment data for all positions into a single DataFrame
    combined_data = {}
    for i, enrichment_pos in enumerate(enrichment[:max_positions]):
        for aa, value in enrichment_pos.items():
            if aa not in combined_data:
                combined_data[aa] = {}
            combined_data[aa][f"P {i + 1}"] = value

    # Convert to DataFrame
    df_combined = pd.DataFrame(combined_data).T  # Transpose to have AAs as rows
    df_combined = df_combined.fillna(0)  # Replace missing values with 0

    # Exclude amino acid 'X'
    if 'X' in df_combined.index:
        df_combined = df_combined.drop(index='X')

    # Sort amino acids alphabetically
    df_combined = df_combined.sort_index()

    # Plot the heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(df_combined, annot=True, fmt=".2f", cmap="coolwarm", cbar=True)
    plt.title("Enrichment for PSI>3 compared to PSI<2")
    plt.xlabel("Positions")
    plt.ylabel("Amino Acids")
    plt.show()



# Main part of the code
if __name__ == '__main__':
    db_file = "Feb2025 NGS for Dec2024 WT resort.db"

    # Load data from the database
    result_list = load_data_from_db(db_file)

    # Filter sequences by PSI thresholds
    population1 = filter_by_psi(result_list, 3, "greater")
    population2 = filter_by_psi(result_list, 2, "less") 

    # Count amino acid occurrences at each position
    counts_pop1 = count_amino_acids_at_positions(population1)
    counts_pop2 = count_amino_acids_at_positions(population2)

    # Normalize the counts
    norm_counts_pop1 = calculate_normalized_counts(counts_pop1)
    norm_counts_pop2 = calculate_normalized_counts(counts_pop2)

    # Calculate enrichment/depletion
    enrichment = calculate_enrichment(norm_counts_pop1, norm_counts_pop2)

    # Plot heatmaps for enrichment/depletion
    plot_combined_heatmap(enrichment)
