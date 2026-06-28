#!/usr/bin/env python3
"""
CG_Solver vs Eigen_CG comparison for larger grids:
  RB 80x38, FT 200x100, LDC 100x100

Reads sor_log.csv and walltime.txt from Output directories.
Also compares VTK field data (p, u, v, T) to verify near-identical results.

Output: pictures/fig74_eigen_comparison_large.pdf
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

BASE = ("/Users/thomaslutz/Documents/Master CSE/2.Semester/CFD Lab"
        "/CFD_Lab_Own/example_cases")
PICTURES = ("/Users/thomaslutz/Documents/Master CSE/2.Semester/CFD Lab"
            "/Project Propsosal/Output_analysis_PCG_SSOR/pictures")

CASES = [
    {
        "label":    "RB 80×38",
        "cg_dir":   "RayleighBenard/RayleighBenard_cg_80x38_1_1_Output_1_1",
        "eigen_dir":"RayleighBenard/RayleighBenard_eigen_cg_80x38_1_1_Output_1_1",
    },
    {
        "label":    "FT 200×100",
        "cg_dir":   "FluidTrap/FluidTrap_200x100_cg_1_1_Output_1_1",
        "eigen_dir":"FluidTrap/FluidTrap_200x100_eigen_cg_1_1_Output_1_1",
    },
    {
        "label":    "LDC 100×100",
        "cg_dir":   "LidDrivenCavity50/LidDrivenCavity100_cg_1_1_Output_1_1",
        "eigen_dir":"LidDrivenCavity50/LidDrivenCavity100_eigen_cg_1_1_Output_1_1",
    },
]

COLOR_CG    = "#4e79a7"
COLOR_EIGEN = "#e15759"


# ---------------------------------------------------------------------------
# Helpers: sor_log / walltime
# ---------------------------------------------------------------------------

def load_log(dirpath):
    p = os.path.join(BASE, dirpath, "sor_log.csv")
    if not os.path.isfile(p):
        return None
    df = pd.read_csv(p)
    df.columns = [c.strip() for c in df.columns]
    return df

def load_walltime(dirpath):
    p = os.path.join(BASE, dirpath, "walltime.txt")
    if not os.path.isfile(p):
        return None
    return float(open(p).read().strip())

def iter_col(df):
    return next((c for c in df.columns if "iter" in c.lower()), None)

def gather(cases):
    rows = []
    for c in cases:
        for key, dname in [("CG_Solver", c["cg_dir"]), ("Eigen_CG", c["eigen_dir"])]:
            log = load_log(dname)
            wt  = load_walltime(dname)
            if log is None or wt is None:
                rows.append(dict(label=c["label"], solver=key,
                                 avg_iter=np.nan, walltime=np.nan,
                                 total_iter=np.nan, tpi=np.nan))
                continue
            col        = iter_col(log)
            avg_iter   = log[col].mean()
            total_iter = log[col].sum()
            tpi        = (wt / total_iter * 1000) if total_iter > 0 else np.nan
            rows.append(dict(label=c["label"], solver=key,
                             avg_iter=avg_iter, walltime=wt,
                             total_iter=total_iter, tpi=tpi))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# VTK field comparison
# ---------------------------------------------------------------------------

def parse_vtk_field(path, field_name, n_components=1):
    """Read a field from VTK FIELD FieldData format. Returns flat numpy array."""
    data = []
    n_vals = 0
    inside = False
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            # Header line: "fieldname n_comp n_tuples dtype"
            parts = stripped.split()
            if len(parts) == 4 and parts[0] == field_name:
                n_vals = int(parts[1]) * int(parts[2])
                inside = True
                continue
            if inside:
                if stripped and not stripped[0].isdigit() and stripped[0] != '-':
                    break
                if stripped:
                    data.extend(float(x) for x in stripped.split())
                if len(data) >= n_vals:
                    break
    arr = np.array(data[:n_vals]) if data else np.array([])
    if n_components > 1 and len(arr) > 0:
        return arr.reshape(-1, n_components)
    return arr

def last_vtk(dirpath):
    pattern = os.path.join(BASE, dirpath, "*.vtk")
    files = sorted(glob.glob(pattern))
    return files[-1] if files else None

def compare_fields(cg_dir, eigen_dir, has_temp=False):
    """Compare last VTK snapshot between CG and Eigen. Returns dict of max/L2 diffs."""
    cg_vtk    = last_vtk(cg_dir)
    eigen_vtk = last_vtk(eigen_dir)
    if cg_vtk is None or eigen_vtk is None:
        return None

    results = {}
    # pressure
    p_cg    = parse_vtk_field(cg_vtk,    "pressure")
    p_eigen = parse_vtk_field(eigen_vtk, "pressure")
    if len(p_cg) and len(p_eigen) and len(p_cg) == len(p_eigen):
        diff = np.abs(p_cg - p_eigen)
        results["p_max"]  = diff.max()
        results["p_l2"]   = np.linalg.norm(diff) / len(diff)**0.5
        results["p_rel"]  = diff.max() / (np.abs(p_cg).max() + 1e-30)

    # velocity (3 components per point)
    v_cg    = parse_vtk_field(cg_vtk,    "velocity", n_components=3)
    v_eigen = parse_vtk_field(eigen_vtk, "velocity", n_components=3)
    if len(v_cg) and len(v_eigen) and v_cg.shape == v_eigen.shape:
        diff_u = np.abs(v_cg[:,0] - v_eigen[:,0])
        diff_v = np.abs(v_cg[:,1] - v_eigen[:,1])
        results["u_max"] = diff_u.max()
        results["v_max"] = diff_v.max()
        results["u_rel"] = diff_u.max() / (np.abs(v_cg[:,0]).max() + 1e-30)
        results["v_rel"] = diff_v.max() / (np.abs(v_cg[:,1]).max() + 1e-30)

    # temperature (only for energy cases)
    if has_temp:
        t_cg    = parse_vtk_field(cg_vtk,    "temperature")
        t_eigen = parse_vtk_field(eigen_vtk, "temperature")
        if len(t_cg) and len(t_eigen) and len(t_cg) == len(t_eigen):
            diff = np.abs(t_cg - t_eigen)
            results["T_max"] = diff.max()
            results["T_rel"] = diff.max() / (np.abs(t_cg).max() + 1e-30)

    return results


# ---------------------------------------------------------------------------
# Build data
# ---------------------------------------------------------------------------

df = gather(CASES)

print("\n=== Summary ===")
print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print()

print("=== Field Comparison (last VTK snapshot) ===")
has_temp = [True, True, False]
field_results = []
for case, ht in zip(CASES, has_temp):
    res = compare_fields(case["cg_dir"], case["eigen_dir"], has_temp=ht)
    print(f"\n{case['label']}:")
    if res:
        for k, v in res.items():
            print(f"  {k:8s} = {v:.3e}")
        field_results.append((case["label"], res))
    else:
        print("  VTK files not found")
        field_results.append((case["label"], {}))


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def grouped_bars(ax, df, col, ylabel, title, fmt="{:.2f}"):
    labels  = [c["label"] for c in CASES]
    cg_vals = [df.loc[(df.label==l) & (df.solver=="CG_Solver"),  col].values[0]
               if len(df.loc[(df.label==l) & (df.solver=="CG_Solver"), col]) else np.nan
               for l in labels]
    eg_vals = [df.loc[(df.label==l) & (df.solver=="Eigen_CG"), col].values[0]
               if len(df.loc[(df.label==l) & (df.solver=="Eigen_CG"), col]) else np.nan
               for l in labels]
    x = np.arange(len(labels))
    w = 0.35
    b1 = ax.bar(x - w/2, cg_vals, w, label="CG_Solver (ours)", color=COLOR_CG)
    b2 = ax.bar(x + w/2, eg_vals, w, label="Eigen_CG",          color=COLOR_EIGEN)
    all_vals = [v for v in cg_vals + eg_vals if not np.isnan(v)]
    vmax = max(all_vals) if all_vals else 1.0
    for bar, v in list(zip(b1, cg_vals)) + list(zip(b2, eg_vals)):
        if not np.isnan(v):
            ax.text(bar.get_x() + bar.get_width()/2, v + vmax*0.01,
                    fmt.format(v), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel); ax.set_title(title)
    ax.legend(fontsize=8)

# Page 1: bar charts
fig1, axes = plt.subplots(1, 3, figsize=(16, 5))
grouped_bars(axes[0], df, "avg_iter",  "Avg. iterations / timestep",
             "Avg. iterations per timestep")
grouped_bars(axes[1], df, "walltime",  "Wall time [s]",
             "Total wall time", fmt="{:.1f}s")
grouped_bars(axes[2], df, "tpi",       "Time per iteration [ms]",
             "Time per iteration [ms]", fmt="{:.3f}")
fig1.suptitle("CG_Solver vs Eigen_CG — larger grids (1×1, warm start)",
              fontsize=12, fontweight="bold")
plt.tight_layout()

# Page 2: per-timestep iteration traces
fig2, axes2 = plt.subplots(1, len(CASES), figsize=(16, 4), sharey=False)
for ax, case in zip(axes2, CASES):
    cg_log = load_log(case["cg_dir"])
    eg_log = load_log(case["eigen_dir"])
    if cg_log is not None:
        col = iter_col(cg_log)
        ax.plot(cg_log["t"], cg_log[col], color=COLOR_CG,    alpha=0.7, lw=0.8, label="CG_Solver")
    if eg_log is not None:
        col = iter_col(eg_log)
        ax.plot(eg_log["t"], eg_log[col], color=COLOR_EIGEN, alpha=0.7, lw=0.8, label="Eigen_CG")
    ax.set_xlabel("t"); ax.set_ylabel("Iterations")
    ax.set_title(case["label"]); ax.legend(fontsize=8)
fig2.suptitle("Iterations per timestep — CG_Solver vs Eigen_CG (large grids)", fontsize=11)
plt.tight_layout()

# Page 3: field differences table
fig3, ax3 = plt.subplots(figsize=(12, 5))
ax3.axis("off")
rows_table = []
for label, res in field_results:
    if not res:
        continue
    rows_table.append([
        label,
        f"{res.get('p_max',np.nan):.2e}",
        f"{res.get('p_rel',np.nan):.2e}",
        f"{res.get('u_max',np.nan):.2e}",
        f"{res.get('u_rel',np.nan):.2e}",
        f"{res.get('v_max',np.nan):.2e}",
        f"{res.get('v_rel',np.nan):.2e}",
        f"{res.get('T_max',np.nan):.2e}" if "T_max" in res else "—",
    ])
cols = ["Case", "p max Δ", "p rel Δ", "u max Δ", "u rel Δ", "v max Δ", "v rel Δ", "T max Δ"]
if not rows_table:
    ax3.text(0.5, 0.5, "No VTK data found", ha="center", va="center")
    rows_table = [["—"]*len(cols)]
tbl = ax3.table(cellText=rows_table, colLabels=cols, loc="center", cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(9)
tbl.auto_set_column_width(list(range(len(cols))))
ax3.set_title("Field differences: CG_Solver vs Eigen_CG (last VTK snapshot)\n"
              "near-identical results confirm both solvers solve the same system",
              fontsize=10, pad=20)
plt.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

os.makedirs(PICTURES, exist_ok=True)
out = os.path.join(PICTURES, "fig74_eigen_comparison_large.pdf")
with PdfPages(out) as pdf:
    pdf.savefig(fig1, bbox_inches="tight")
    pdf.savefig(fig2, bbox_inches="tight")
    pdf.savefig(fig3, bbox_inches="tight")
plt.close("all")

print(f"\nSaved: {out}")
print("Done.")
