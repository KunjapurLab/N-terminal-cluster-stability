#!/usr/bin/env python3

import pickle
import sqlite3
import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt
import os 

def find_existing_file(candidates):
    seen = set()
    for path in candidates:
        normalized = os.path.normpath(path)
        if normalized in seen:
            continue
        seen.add(normalized)
        if os.path.isfile(normalized):
            return normalized
    return None

############################################
# Data Loading and Preprocessing Functions
############################################
def load_data(db_file):
    """
    Connect to the SQLite database and load Sequence, B1, B2, B3, B4, and PSI from ResultList.
    """
    conn = sqlite3.connect(db_file)
    query = """
        SELECT 
            Sequence, 
            B1, 
            B2, 
            B3, 
            B4, 
            PSI
        FROM ResultList
    """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

def filter_sequences_by_reads(df, min_reads=20):
    """
    Compute total_reads = B1 + B2 + B3 + B4 and filter out rows with fewer than min_reads.
    """
    df['Total_reads'] = df['B1'] + df['B2'] + df['B3'] + df['B4']
    return df[df['Total_reads'] >= min_reads].copy()

def one_hot_encode_sequences(sequences):
    """
    One-hot encode each 5-residue sequence into a 100-dimensional feature vector.
    Returns the feature matrix X and the list of feature names.
    """
    aa_list = [
        'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G',
        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
        'T', 'W', 'Y', 'V'
    ]
    aa_to_idx = {aa: i for i, aa in enumerate(aa_list)}
    num_samples = len(sequences)
    num_positions = 5
    num_aa = 20

    X = np.zeros((num_samples, num_positions * num_aa), dtype=np.float32)
    feature_names = []
    for pos in range(num_positions):
        for aa in aa_list:
            feature_names.append(f"Pos{pos+1}_{aa}")

    for row_idx, seq in enumerate(sequences):
        if len(seq) != num_positions:
            continue
        for pos in range(num_positions):
            aa = seq[pos]
            if aa in aa_to_idx:
                col_idx = pos * num_aa + aa_to_idx[aa]
                X[row_idx, col_idx] = 1.0

    return X, feature_names

############################################
# Main Script
############################################
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    cache_file = os.path.join(script_dir, "shap_interaction_values.pkl")

    # 1. Import the pre-trained model
    model_candidates = [
        os.path.join(script_dir, "20Feb2025 WT BL21 model.pkl"),
        os.path.join(script_dir, "N-FIVE WT model.pkl"),
        os.path.join(repo_root, "N-FIVE", "N-FIVE WT model.pkl"),
    ]
    model_path = find_existing_file(model_candidates)
    if not model_path:
        raise FileNotFoundError(
            "Model file not found. Searched:\n" + "\n".join(f"- {p}" for p in model_candidates)
        )

    with open(model_path, "rb") as file:
        final_model = pickle.load(file)
    print(f"Model loaded successfully from: {model_path}")

    # 2. Load and preprocess data
    db_candidates = [
        os.path.join(script_dir, "Feb2025 NGS for Dec2024 WT resort.db"),
        os.path.join(repo_root, "Databases", "Feb2025 NGS for Dec2024 WT resort.db"),
    ]
    db_file = find_existing_file(db_candidates)
    if not db_file:
        raise FileNotFoundError(
            "Database file not found. Searched:\n" + "\n".join(f"- {p}" for p in db_candidates)
        )
    min_reads = 20

    print(f"Loading data from {db_file}...")
    df = load_data(db_file)
    print(f"Loaded {len(df)} rows.")

    df = filter_sequences_by_reads(df, min_reads)
    print(f"{len(df)} rows remain after filtering for >= {min_reads} reads.")

    sequences = df["Sequence"].values
    X, feature_names = one_hot_encode_sequences(sequences)
    print("Feature matrix shape:", X.shape)

    # 3. Sample the data for SHAP interaction analysis if necessary
    sample_size = 2000
    if X.shape[0] > sample_size:
        np.random.seed(42)
        idx = np.random.choice(X.shape[0], sample_size, replace=False)
        X_sample = X[idx]
    else:
        X_sample = X

    # 4. Create SHAP TreeExplainer and compute or load interaction values
    explainer = shap.TreeExplainer(final_model)
    
    if os.path.exists(cache_file):
        print("Loading cached SHAP interaction values...")
        with open(cache_file, "rb") as f:
            shap_interaction_values = pickle.load(f)
    else:
        print("Calculating SHAP interaction values...")
        shap_interaction_values = explainer.shap_interaction_values(X_sample)
        with open(cache_file, "wb") as f:
            pickle.dump(shap_interaction_values, f)      
    
    # Convert feature_names to a NumPy array for compatibility with SHAP plotting
    feature_names = np.array(feature_names)
    
    # 5. Visualize the interaction effects with a summary plot
    print("Generating SHAP interaction summary plot...")
    shap.summary_plot(shap_interaction_values, X_sample, feature_names=feature_names)

if __name__ == "__main__":
    main()
