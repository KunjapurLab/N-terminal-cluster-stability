# N-FIVE Web Studio

Local web app for:
- Batch PSI prediction
- SHAP attribution for a single sequence

## Run (Windows PowerShell)

From the repository root:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python "N-FIVE/webapp/app.py"
```

Open:

`http://127.0.0.1:8000`

## Notes

- The app auto-locates `N-FIVE WT model.pkl` in `N-FIVE/`.
- Input sequences must be exactly 5 amino acids.
