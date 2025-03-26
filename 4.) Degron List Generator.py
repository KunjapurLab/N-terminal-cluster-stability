import sqlite3
import os

def fetch_result_list_from_db(db_file):
    """Fetches the result list from the SQLite database."""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM ResultList")
    result_list = cursor.fetchall()

    conn.close()
    return result_list

def write_output(output_file, common_counts, num_reads_list, result_list, avg_reads_per_bin):
    """Writes the result to a text file, including sequence counts and average reads per bin."""
    with open(output_file, 'w') as output_file:
        output_file.write('Across 4 bins of data, here are the read counts and PSI for each 5 AA sequence.\n\n')
        for entry in result_list:
            output_file.write(f'Sequence: {entry[0]}, B1: {entry[1]}, B2: {entry[2]}, B3: {entry[3]}, B4: {entry[4]}, PSI: {entry[5]}\n')

        output_file.write("\nSummary of Common Counts and Reads:\n")
        for i in range(4):
            output_file.write(f'B{i+1}: Common Count = {common_counts[i]}, Total Reads = {num_reads_list[i]}\n')

        output_file.write("\nNumber of unique sequences per bin:\n")
        for i in range(4):
            unique_sequences = sum(1 for entry in result_list if entry[i+1] > 0)  # Count sequences with non-zero read count in bin
            output_file.write(f'B{i+1}: {unique_sequences} unique sequences\n')

        output_file.write("\nAverage number of reads per sequence per bin:\n")
        for i in range(4):
            unique_sequences = sum(1 for entry in result_list if entry[i+1] > 0)  # Number of unique sequences in bin
            total_reads = sum(entry[i+1] for entry in result_list)  # Total reads in the bin
            
            # Avoid division by zero if there are no sequences in the bin
            if unique_sequences > 0:
                avg_reads = total_reads / unique_sequences
            else:
                avg_reads = 0
            
            output_file.write(f'B{i+1}: Average Reads per Sequence = {avg_reads:.2f}\n')

def fetch_common_counts_and_reads(db_file):
    """Fetches common counts and reads from the database to summarize."""
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    cursor.execute("SELECT B1, B2, B3, B4 FROM ResultList")
    counts = list(zip(*cursor.fetchall()))
    common_counts = [sum(counts[i]) for i in range(4)]
    num_reads_list = [sum(counts[i]) for i in range(4)]

    conn.close()
    return common_counts, num_reads_list

# Main execution
if __name__ == '__main__':
    db_file = "10M ClpS- consensus 22Nov.db"  # Path to your SQLite database file
    output_file = " 22Nov2024 10M ClpS- Degron List"

    # Fetch result list from the database (no recalculation of PSI)
    result_list = fetch_result_list_from_db(db_file)

    # Fetch common counts and read counts
    common_counts, num_reads_list = fetch_common_counts_and_reads(db_file)

    # Write the results to a text file (no PSI recalculation)
    write_output(output_file + ".txt", common_counts, num_reads_list, result_list, None)

    print("Run completed.")
