import tkinter as tk
from tkinter import filedialog, messagebox

def split_fasta_by_barcode(input_fasta, barcode_to_filename):
    """
    Splits a FASTA file into separate files based on a specific 3-letter barcode
    at the beginning of each sequence line.

    Parameters:
    - input_fasta (str): Path to the input FASTA file.
    - barcode_to_filename (dict): Dictionary where keys are 3-letter barcodes
      and values are output filenames.
    """
    # Open each output file for writing
    output_files = {barcode: open(filename, 'w') for barcode, filename in barcode_to_filename.items()}

    # Read the input FASTA file
    with open(input_fasta, 'r') as f:
        while True:
            header = f.readline().strip()
            sequence = f.readline().strip()
            
            # Break the loop if we reach the end of the file
            if not header or not sequence:
                break

            # Extract the 3-letter barcode from the beginning of the sequence
            barcode = sequence[:3]

            # Write the header and sequence to the appropriate file
            if barcode in output_files:
                output_files[barcode].write(f"{header}\n{sequence}\n")

    # Close all output files
    for file in output_files.values():
        file.close()

    print("FASTA file splitting completed.")
    messagebox.showinfo("Completed", "FASTA file splitting completed.")

# GUI function
def run_gui():
    def add_entry():
        barcode = barcode_entry.get().strip().upper()
        filename = filename_entry.get().strip()
        if len(barcode) == 3 and filename:
            barcode_to_filename[barcode] = filename
            barcode_entry.delete(0, tk.END)
            filename_entry.delete(0, tk.END)
            update_barcode_list()
        else:
            messagebox.showwarning("Input Error", "Please enter a 3-letter barcode and a filename.")

    def update_barcode_list():
        barcode_listbox.delete(0, tk.END)
        for barcode, filename in barcode_to_filename.items():
            barcode_listbox.insert(tk.END, f"{barcode}: {filename}")

    def select_input_file():
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        if file_path:
            input_file_entry.delete(0, tk.END)
            input_file_entry.insert(0, file_path)

    def start_splitting():
        input_fasta = input_file_entry.get().strip()
        if input_fasta and barcode_to_filename:
            split_fasta_by_barcode(input_fasta, barcode_to_filename)
        else:
            messagebox.showwarning("Input Error", "Please select an input file and enter at least one barcode and filename.")

    # Create the main window
    root = tk.Tk()
    root.title("FASTA Splitter")

    # Input FASTA file section
    tk.Label(root, text="Input FASTA file:").grid(row=0, column=0, padx=5, pady=5)
    input_file_entry = tk.Entry(root, width=40)
    input_file_entry.grid(row=0, column=1, padx=5, pady=5)
    tk.Button(root, text="Browse", command=select_input_file).grid(row=0, column=2, padx=5, pady=5)

    # Barcode and filename entry section
    tk.Label(root, text="Barcode (3 letters):").grid(row=1, column=0, padx=5, pady=5)
    barcode_entry = tk.Entry(root, width=10)
    barcode_entry.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    
    tk.Label(root, text="Output Filename:").grid(row=2, column=0, padx=5, pady=5)
    filename_entry = tk.Entry(root, width=30)
    filename_entry.grid(row=2, column=1, padx=5, pady=5)

    tk.Button(root, text="Add Barcode", command=add_entry).grid(row=2, column=2, padx=5, pady=5)

    # List of added barcodes and filenames
    tk.Label(root, text="Barcodes and Filenames:").grid(row=3, column=0, columnspan=3, padx=5, pady=5)
    barcode_listbox = tk.Listbox(root, width=50)
    barcode_listbox.grid(row=4, column=0, columnspan=3, padx=5, pady=5)

    # Start button
    tk.Button(root, text="Start Splitting", command=start_splitting).grid(row=5, column=0, columnspan=3, padx=5, pady=5)

    # Initialize the barcode to filename mapping
    barcode_to_filename = {}

    root.mainloop()

# Run the GUI
run_gui()