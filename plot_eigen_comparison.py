#!/usr/bin/env python3
"""
CG_Solver vs Eigen_CG comparison: avg iterations/step, walltime, time/iter.

Usage:
    python3 plot_eigen_comparison.py

Reads sor_log.csv and walltime.txt from the Output directories of
CG_STANDARD and EIGEN_CG runs (1x1 only, Eigen is serial).

Output: pictures/fig73_eigen_comparison.pdf
"""

import os
import pandas as pd
import numpy as np
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
        "label":    "RB 40×18",
        "cg_dir":   "RayleighBenard/RayleighBenard_cg_1_1_Output_1_1",
        "eigen_dir":"RayleighBenard/RayleighBenard_eigen_cg_1_1_Output_1_1",
    },
    {
        "label":    "FT 100×50",
        "cg_dir":   "FluidTrap/FluidTrap_s100_cg_1_1_Output_1_1",
        "eigen_dir":"FluidTrap/FluidTrap_s100_eigen_cg_1_1_Output_1_1",
    },
    {
        "label":    "LDC 50×50",
        "cg_dir":   "LidDrivenCavity50/LidDrivenCavity50_cg_1_1_Output_1_1",
        "eigen_dir":"LidDrivenCavity50/LidDrivenCavity50_eigen_cg_1_1_Output_1_1",
    },
]

COLOR_CG    = "#4e79a7"
COLOR_EIGEN = "#e15759"


# ---------------------------------------------------------------------------
# Helpers
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
    """Return name of the iteration column."""
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
# Plot helpers
# ---------------------------------------------------------------------------

def grouped_bars(ax, df, col, ylabel, title, fmt="{:.2f}"):
    labels  = [c["label"] for c in CASES]
    cg_vals = [df.loc[(df.label==l) & (df.solver=="CG_Solver"),  col].values[0]
               if len(df.loc[(df.label==l) & (df.solver=="CG_Solver"),  col]) else np.nan
               for l in labels]
    eg_vals = [df.loc[(df.label==l) & (df.solver=="Eigen_CG"), col].values[0]
               if len(df.loc[(df.label==l) & (df.solver=="Eigen_CG"), col]) else np.nan
               for l in labels]

    x = np.arange(len(labels))
    w = 0.35
    b1 = ax.bar(x - w/2, cg_vals,  w, label="CG_Solver (ours)", color=COLOR_CG)
    b2 = ax.bar(x + w/2, eg_vals,  w, label="Eigen_CG",          color=COLOR_EIGEN)

    # Value labels on bars
    all_vals = [v for v in cg_vals + eg_vals if not np.isnan(v)]
    vmax = max(all_vals) if all_vals else 1.0
    for bar, v in list(zip(b1, cg_vals)) + list(zip(b2, eg_vals)):
        if not np.isnan(v):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    v + vmax * 0.01,
                    fmt.format(v), ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)


# ---------------------------------------------------------------------------
# Build data and plots
# ---------------------------------------------------------------------------

df = gather(CASES)

# Print table to console for quick sanity check
print("\n=== Summary ===")
print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print()

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

grouped_bars(axes[0], df, "avg_iter",
             "Avg. iterations / timestep",
             "Avg. iterations per timestep\n(lower = faster convergence)")

grouped_bars(axes[1], df, "walltime",
             "Wall time [s]",
             "Total wall time",
             fmt="{:.1f}s")

grouped_bars(axes[2], df, "tpi",
             "Time per iteration [ms]",
             "Time per iteration [ms]\n(lower = cheaper matvec)",
             fmt="{:.3f}")

fig.suptitle("CG_Solver (hand-written) vs Eigen_CG (library) — 1×1, warm start",
             fontsize=12, fontweight="bold")

# Note about iteration count caveat
fig.text(0.5, -0.02,
         "Note: iteration counts may differ by ±1 due to different internal convergence\n"
         "check placement (CG_Solver: pre-check; Eigen: post-check). "
         "Time/iter comparison is the most meaningful metric.",
         ha="center", fontsize=8, color="gray")

plt.tight_layout()

# ---------------------------------------------------------------------------
# Per-timestep iteration traces (one page per case)
# ---------------------------------------------------------------------------

fig2, axes2 = plt.subplots(1, len(CASES), figsize=(16, 4), sharey=False)

for ax, case in zip(axes2, CASES):
    cg_log  = load_log(case["cg_dir"])
    eg_log  = load_log(case["eigen_dir"])

    if cg_log is not None:
        col = iter_col(cg_log)
        ax.plot(cg_log["t"], cg_log[col],
                color=COLOR_CG,    alpha=0.7, lw=0.8, label="CG_Solver")
    if eg_log is not None:
        col = iter_col(eg_log)
        ax.plot(eg_log["t"], eg_log[col],
                color=COLOR_EIGEN, alpha=0.7, lw=0.8, label="Eigen_CG")

    ax.set_xlabel("t")
    ax.set_ylabel("Iterations")
    ax.set_title(case["label"])
    ax.legend(fontsize=8)

fig2.suptitle("Iterations per timestep over time — CG_Solver vs Eigen_CG",
              fontsize=11)
plt.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

os.makedirs(PICTURES, exist_ok=True)
out = os.path.join(PICTURES, "fig73_eigen_comparison.pdf")
with PdfPages(out) as pdf:
    pdf.savefig(fig,  bbox_inches="tight")
    pdf.savefig(fig2, bbox_inches="tight")
plt.close("all")

print(f"Saved: {out}")
print("Done.")
