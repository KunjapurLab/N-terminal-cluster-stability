#!/usr/bin/env python3

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import random
import numpy as np
import joblib
import sys
import os
import xgboost as xgb

# Helper to bundle resources with PyInstaller
def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

def find_model_file():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    candidates = [
        resource_path("20Feb2025 WT BL21 model.pkl"),
        resource_path("N-FIVE WT model.pkl"),
        os.path.join(script_dir, "20Feb2025 WT BL21 model.pkl"),
        os.path.join(script_dir, "N-FIVE WT model.pkl"),
        os.path.join(repo_root, "N-FIVE", "N-FIVE WT model.pkl"),
    ]

    seen = set()
    for path in candidates:
        normalized = os.path.normpath(path)
        if normalized in seen:
            continue
        seen.add(normalized)
        if os.path.isfile(normalized):
            return normalized
    return None, candidates


# Load XGBoost model
model_path = find_model_file()
if isinstance(model_path, tuple):
    _, searched_paths = model_path
    raise FileNotFoundError(
        "Model file not found. Searched:\n" + "\n".join(f"- {p}" for p in searched_paths)
    )

try:
    xgb_model = joblib.load(model_path)
except Exception as e:
    raise Exception(f"Model loading failed from {model_path}: {e}")

# Amino acids (fixed order)
AA_LIST = list('ARNDCEQGHILKMFPSTWYV')
AA_TO_IDX = {aa: idx for idx, aa in enumerate(AA_LIST)}
NUM_POSITIONS = 5
NUM_AA = len(AA_LIST)

# ------------------------------ Encoding ------------------------------ #

def encode_single_sequence(seq):
    if len(seq) != NUM_POSITIONS:
        raise ValueError("Sequence must be exactly 5 amino acids long.")
    encoded = np.zeros(NUM_POSITIONS * NUM_AA, dtype=np.float32)
    for pos, aa in enumerate(seq.upper()):
        if aa in AA_TO_IDX:
            encoded[pos * NUM_AA + AA_TO_IDX[aa]] = 1
        else:
            raise ValueError(f"Invalid amino acid '{aa}'")
    return encoded.reshape(1, -1)

# ------------------------------ Prediction ------------------------------ #

def predict_psi(seq):
    encoded_seq = encode_single_sequence(seq)
    psi = xgb_model.predict(encoded_seq)[0]
    return psi

def predict_psis_for_list(sequences):
    results = []
    for seq in sequences:
        seq = seq.strip().upper()
        try:
            psi = predict_psi(seq)
            results.append((seq, f"{psi:.4f}"))
        except Exception as e:
            results.append((seq, f"Error: {e}"))
    return results

# ------------------------------ Suggestion ------------------------------ #

def random_sequence(constraint_pos=None, constraint_residue=None):
    seq = [random.choice(AA_LIST) for _ in range(NUM_POSITIONS)]
    if constraint_pos and constraint_residue:
        seq[constraint_pos - 1] = constraint_residue.upper()
    return ''.join(seq)

def suggest_sequences(target_psi, constraint_pos=None, constraint_residue=None, num_suggestions=5, num_candidates=5000):
    sequences_set = set()
    sequences = []
    attempts = 0
    max_attempts = num_candidates * 2  # prevent infinite loop

    while len(sequences) < num_candidates and attempts < max_attempts:
        seq = random_sequence(constraint_pos, constraint_residue)
        if seq not in sequences_set:
            sequences_set.add(seq)
            sequences.append(seq)
        attempts += 1

    X_encoded = np.vstack([encode_single_sequence(seq) for seq in sequences])
    preds = xgb_model.predict(X_encoded)
    sorted_seqs = sorted(zip(sequences, preds), key=lambda x: abs(x[1] - target_psi))

    # Ensure uniqueness in the final suggested list
    unique_suggestions = []
    seen = set()
    for seq, psi in sorted_seqs:
        if seq not in seen:
            unique_suggestions.append((seq, psi))
            seen.add(seq)
        if len(unique_suggestions) == num_suggestions:
            break

    return unique_suggestions

# ------------------------------ GUI Callbacks ------------------------------ #

def copy_to_clipboard(widget):
    root.clipboard_clear()
    root.clipboard_append(widget.get("1.0", tk.END).strip())
    messagebox.showinfo("Copied", "Copied to clipboard!")

def validate_numeric(entry):
    try:
        float(entry.get())
        entry.config(background='white')
    except ValueError:
        entry.config(background='salmon')

def predict_callback():
    sequences = predict_text.get("1.0", tk.END).strip().splitlines()
    if not sequences:
        messagebox.showwarning("Input needed", "Enter sequences to predict.")
        return
    results = predict_psis_for_list(sequences)
    predict_output.delete("1.0", tk.END)
    for seq, psi in results:
        predict_output.insert(tk.END, f"{seq}: {psi}\n")

def suggest_callback():
    try:
        target = float(suggest_target_entry.get())
    except ValueError:
        messagebox.showerror("Error", "Enter a valid numeric PSI.")
        return
    suggest_output.delete("1.0", tk.END)
    suggestions = suggest_sequences(target)
    for seq, psi in suggestions:
        suggest_output.insert(tk.END, f"{seq} ({psi:.4f})\n")

def constrained_suggest_callback():
    try:
        target = float(constraint_target_entry.get())
        pos = int(constraint_pos_entry.get())
        residue = constraint_residue_entry.get().strip().upper()
        num_suggestions = int(num_suggestions_entry.get())
        num_candidates = int(num_candidates_entry.get())

        if pos not in range(1, 6) or residue not in AA_LIST:
            raise ValueError("Position or residue invalid.")
        if num_suggestions <= 0 or num_candidates <= 0:
            raise ValueError("Suggestions and candidates must be positive integers.")
    except ValueError as e:
        messagebox.showerror("Error", f"Invalid input: {e}")
        return

    constraint_output.delete("1.0", tk.END)
    suggestions = suggest_sequences(
        target,
        constraint_pos=pos,
        constraint_residue=residue,
        num_suggestions=num_suggestions,
        num_candidates=num_candidates
    )
    for seq, psi in suggestions:
        constraint_output.insert(tk.END, f"{seq} ({psi:.4f})\n")


# ------------------------------ GUI ------------------------------ #

root = tk.Tk()
root.title("PSI Predictor (XGBoost)")

style = ttk.Style(root)
style.theme_use('clam')

# Notebook tabs
tabs = ttk.Notebook(root)
tabs.pack(expand=True, fill="both")

# Predict Tab
tab_predict = ttk.Frame(tabs)
tabs.add(tab_predict, text="Predict PSI")
ttk.Label(tab_predict, text="Enter sequences:").pack(anchor="w", padx=10, pady=5)
predict_text = scrolledtext.ScrolledText(tab_predict, width=50, height=10, font=('Consolas', 10))
predict_text.pack(padx=10, pady=5)
ttk.Button(tab_predict, text="Predict", command=predict_callback).pack(padx=10, pady=5)
predict_output = scrolledtext.ScrolledText(tab_predict, width=50, height=10, font=('Consolas', 10))
predict_output.pack(padx=10, pady=5)
ttk.Button(tab_predict, text="Copy Results", command=lambda: copy_to_clipboard(predict_output)).pack(padx=10, pady=5)

# Suggest Tab
tab_suggest = ttk.Frame(tabs)
tabs.add(tab_suggest, text="Suggest by PSI")
frame_suggest = ttk.Frame(tab_suggest)
frame_suggest.pack(padx=10, pady=10)
ttk.Label(frame_suggest, text="Target PSI:").pack(side="left", padx=5)
suggest_target_entry = ttk.Entry(frame_suggest, width=10)
suggest_target_entry.pack(side="left")
suggest_target_entry.bind("<FocusOut>", lambda _: validate_numeric(suggest_target_entry))
ttk.Button(frame_suggest, text="Suggest", command=suggest_callback).pack(side="left", padx=5)
suggest_output = scrolledtext.ScrolledText(tab_suggest, width=50, height=10, font=('Consolas', 10))
suggest_output.pack(padx=10, pady=5)
ttk.Button(tab_suggest, text="Copy Suggestions", command=lambda: copy_to_clipboard(suggest_output)).pack(padx=10, pady=5)



# Constrained Suggest Tab
# --- GUI Constraint Tab Updated --- #
tab_constraint = ttk.Frame(tabs)
tabs.add(tab_constraint, text="Suggest with Constraint")

frame_constraint = ttk.Frame(tab_constraint)
frame_constraint.pack(padx=10, pady=10)

ttk.Label(frame_constraint, text="Target PSI:").grid(row=0, column=0, padx=5, pady=2)
constraint_target_entry = ttk.Entry(frame_constraint, width=10)
constraint_target_entry.grid(row=0, column=1)
constraint_target_entry.bind("<FocusOut>", lambda _: validate_numeric(constraint_target_entry))

ttk.Label(frame_constraint, text="Position (1-5):").grid(row=1, column=0, padx=5, pady=2)
constraint_pos_entry = ttk.Entry(frame_constraint, width=10)
constraint_pos_entry.grid(row=1, column=1)

ttk.Label(frame_constraint, text="Residue:").grid(row=2, column=0, padx=5, pady=2)
constraint_residue_entry = ttk.Entry(frame_constraint, width=10)
constraint_residue_entry.grid(row=2, column=1)

ttk.Label(frame_constraint, text="Number of suggestions:").grid(row=3, column=0, padx=5, pady=2)
num_suggestions_entry = ttk.Entry(frame_constraint, width=10)
num_suggestions_entry.grid(row=3, column=1)
num_suggestions_entry.insert(0, "5")

ttk.Label(frame_constraint, text="Number of candidates:").grid(row=4, column=0, padx=5, pady=2)
num_candidates_entry = ttk.Entry(frame_constraint, width=10)
num_candidates_entry.grid(row=4, column=1)
num_candidates_entry.insert(0, "5000")

ttk.Button(frame_constraint, text="Suggest", command=constrained_suggest_callback).grid(row=5, column=0, columnspan=2, pady=5)

constraint_output = scrolledtext.ScrolledText(tab_constraint, width=50, height=10, font=('Consolas', 10))
constraint_output.pack(padx=10, pady=5)

ttk.Button(tab_constraint, text="Copy Suggestions", command=lambda: copy_to_clipboard(constraint_output)).pack(padx=10, pady=5)

# Run the GUI
root.mainloop()
