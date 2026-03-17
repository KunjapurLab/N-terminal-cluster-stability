#!/usr/bin/env python3
import os
from typing import Any

import joblib
import numpy as np
import shap
from flask import Flask, jsonify, render_template, request

AA_LIST = list("ARNDCEQGHILKMFPSTWYV")
AA_TO_IDX = {aa: idx for idx, aa in enumerate(AA_LIST)}
NUM_POSITIONS = 5
NUM_AA = len(AA_LIST)

app = Flask(__name__)


def find_model_file() -> str:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    n_five_dir = os.path.dirname(script_dir)
    repo_root = os.path.dirname(n_five_dir)
    candidates = [
        os.path.join(n_five_dir, "N-FIVE WT model.pkl"),
        os.path.join(n_five_dir, "20Feb2025 WT BL21 model.pkl"),
        os.path.join(script_dir, "N-FIVE WT model.pkl"),
        os.path.join(script_dir, "20Feb2025 WT BL21 model.pkl"),
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

    raise FileNotFoundError(
        "Model file not found. Searched:\n" + "\n".join(f"- {p}" for p in candidates)
    )


def encode_single_sequence(seq: str) -> np.ndarray:
    seq = seq.strip().upper()
    if len(seq) != NUM_POSITIONS:
        raise ValueError("Sequence must be exactly 5 amino acids long.")

    encoded = np.zeros(NUM_POSITIONS * NUM_AA, dtype=np.float32)
    for pos, aa in enumerate(seq):
        if aa not in AA_TO_IDX:
            raise ValueError(f"Invalid amino acid '{aa}'.")
        encoded[pos * NUM_AA + AA_TO_IDX[aa]] = 1.0
    return encoded.reshape(1, -1)


MODEL_PATH = find_model_file()
MODEL = joblib.load(MODEL_PATH)
EXPLAINER = shap.TreeExplainer(MODEL)


def _as_float(value: Any) -> float:
    if isinstance(value, (np.floating, float, int)):
        return float(value)
    if isinstance(value, np.ndarray):
        return float(value.reshape(-1)[0])
    return float(value)


@app.get("/")
def index():
    return render_template("index.html", model_path=MODEL_PATH.replace("\\", "/"))


@app.get("/api/health")
def health():
    return jsonify({"status": "ok", "model_path": MODEL_PATH})


@app.post("/api/predict")
def predict():
    payload = request.get_json(silent=True) or {}
    sequences = payload.get("sequences", [])
    if not isinstance(sequences, list):
        return jsonify({"error": "Expected 'sequences' as a list."}), 400

    results = []
    for raw_seq in sequences:
        seq = str(raw_seq).strip().upper()
        if not seq:
            continue
        try:
            encoded = encode_single_sequence(seq)
            psi = _as_float(MODEL.predict(encoded)[0])
            results.append({"sequence": seq, "psi": round(psi, 4)})
        except Exception as exc:
            results.append({"sequence": seq, "error": str(exc)})

    return jsonify({"results": results, "count": len(results)})


@app.post("/api/shap")
def shap_values():
    payload = request.get_json(silent=True) or {}
    seq = str(payload.get("sequence", "")).strip().upper()
    if not seq:
        return jsonify({"error": "Please provide a sequence."}), 400

    try:
        encoded = encode_single_sequence(seq)
        psi = _as_float(MODEL.predict(encoded)[0])

        shap_values_raw = EXPLAINER.shap_values(encoded)
        if isinstance(shap_values_raw, list):
            shap_vector = np.asarray(shap_values_raw[0]).reshape(-1)
        else:
            shap_vector = np.asarray(shap_values_raw).reshape(-1)

        expected_value = EXPLAINER.expected_value
        if isinstance(expected_value, (list, np.ndarray)):
            base_value = _as_float(np.asarray(expected_value).reshape(-1)[0])
        else:
            base_value = _as_float(expected_value)

        position_breakdown = []
        for pos in range(NUM_POSITIONS):
            aa = seq[pos]
            aa_idx = AA_TO_IDX[aa]
            selected_idx = pos * NUM_AA + aa_idx
            pos_slice = slice(pos * NUM_AA, (pos + 1) * NUM_AA)
            position_breakdown.append(
                {
                    "position": pos + 1,
                    "residue": aa,
                    "selected_feature": f"Pos{pos + 1}_{aa}",
                    "selected_contribution": round(_as_float(shap_vector[selected_idx]), 6),
                    "position_total": round(_as_float(np.sum(shap_vector[pos_slice])), 6),
                }
            )

        feature_rows = []
        for pos in range(NUM_POSITIONS):
            for aa in AA_LIST:
                idx = pos * NUM_AA + AA_TO_IDX[aa]
                feature_rows.append(
                    {"feature": f"Pos{pos + 1}_{aa}", "contribution": _as_float(shap_vector[idx])}
                )
        top_features = sorted(feature_rows, key=lambda row: abs(row["contribution"]), reverse=True)[:12]
        for row in top_features:
            row["contribution"] = round(row["contribution"], 6)

        return jsonify(
            {
                "sequence": seq,
                "psi": round(psi, 4),
                "base_value": round(base_value, 6),
                "position_breakdown": position_breakdown,
                "top_features": top_features,
            }
        )
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8000, debug=False)
