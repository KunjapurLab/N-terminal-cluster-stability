function setStatus(id, message, isError = false) {
  const el = document.getElementById(id);
  el.textContent = message;
  el.style.color = isError ? "#b91c1c" : "";
}

function switchTab(tabName) {
  document.querySelectorAll(".tab").forEach((btn) => {
    const active = btn.dataset.tab === tabName;
    btn.classList.toggle("active", active);
    btn.setAttribute("aria-selected", active ? "true" : "false");
  });
  document.querySelectorAll(".panel").forEach((panel) => {
    panel.classList.toggle("active-panel", panel.id === tabName);
  });
}

async function runPrediction() {
  const input = document.getElementById("predict-input").value;
  const sequences = input
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line.length > 0);

  if (sequences.length === 0) {
    setStatus("predict-status", "Enter at least one sequence.", true);
    return;
  }

  setStatus("predict-status", "Running predictions...");
  const tbody = document.querySelector("#predict-table tbody");
  tbody.innerHTML = "";

  try {
    const res = await fetch("/api/predict", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequences }),
    });
    const data = await res.json();
    if (!res.ok) {
      throw new Error(data.error || "Prediction request failed.");
    }

    for (const row of data.results) {
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>${row.sequence || ""}</td>
        <td>${row.psi ?? "-"}</td>
        <td>${row.error ? row.error : "ok"}</td>
      `;
      tbody.appendChild(tr);
    }
    setStatus("predict-status", `Completed ${data.count} sequence(s).`);
  } catch (err) {
    setStatus("predict-status", err.message, true);
  }
}

function renderPositionBars(rows) {
  const container = document.getElementById("position-bars");
  container.innerHTML = "";
  const maxAbs = Math.max(0.000001, ...rows.map((r) => Math.abs(r.position_total)));

  rows.forEach((row) => {
    const width = Math.max(2, Math.round((Math.abs(row.position_total) / maxAbs) * 100));
    const signClass = row.position_total >= 0 ? "pos" : "neg";
    const entry = document.createElement("div");
    entry.className = "bar-row";
    entry.innerHTML = `
      <div>P${row.position} (${row.residue})</div>
      <div class="bar-track"><div class="bar-fill ${signClass}" style="width:${width}%"></div></div>
      <div>${row.position_total.toFixed(4)}</div>
    `;
    container.appendChild(entry);
  });
}

function renderShapTable(features) {
  const tbody = document.querySelector("#shap-table tbody");
  tbody.innerHTML = "";
  features.forEach((feature) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${feature.feature}</td>
      <td>${feature.contribution.toFixed(6)}</td>
    `;
    tbody.appendChild(tr);
  });
}

async function runShap() {
  const sequence = document.getElementById("shap-input").value.trim().toUpperCase();
  if (!sequence) {
    setStatus("shap-status", "Enter one sequence.", true);
    return;
  }

  setStatus("shap-status", "Computing SHAP values...");
  document.getElementById("kpi-psi").textContent = "-";
  document.getElementById("kpi-base").textContent = "-";

  try {
    const res = await fetch("/api/shap", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequence }),
    });
    const data = await res.json();
    if (!res.ok) {
      throw new Error(data.error || "SHAP request failed.");
    }

    document.getElementById("kpi-psi").textContent = data.psi.toFixed(4);
    document.getElementById("kpi-base").textContent = data.base_value.toFixed(6);
    renderPositionBars(data.position_breakdown || []);
    renderShapTable(data.top_features || []);
    setStatus("shap-status", `SHAP completed for ${data.sequence}.`);
  } catch (err) {
    setStatus("shap-status", err.message, true);
  }
}

document.querySelectorAll(".tab").forEach((tab) => {
  tab.addEventListener("click", () => switchTab(tab.dataset.tab));
});

document.getElementById("predict-btn").addEventListener("click", runPrediction);
document.getElementById("predict-clear").addEventListener("click", () => {
  document.getElementById("predict-input").value = "";
  document.querySelector("#predict-table tbody").innerHTML = "";
  setStatus("predict-status", "");
});
document.getElementById("shap-btn").addEventListener("click", runShap);

document.getElementById("shap-input").addEventListener("keydown", (event) => {
  if (event.key === "Enter") {
    event.preventDefault();
    runShap();
  }
});
