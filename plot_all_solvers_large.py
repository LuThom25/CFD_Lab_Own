#!/usr/bin/env python3
"""
SOR vs CG vs PCG_SSOR vs Eigen_CG — full 4-solver comparison, large grids (1x1).

Cases:
  RB  80x36:   all 4 solvers
  FT  200x100: all 4 solvers
  LDC 100x100: all 4 solvers

Output: pictures/fig76_all_solvers_large.pdf
"""

import os
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
        "label": "RB 80×36",
        "solvers": {
            "SOR":      "RayleighBenard/RayleighBenard_sor_80x36_1_1_Output_1_1",
            "CG":       "RayleighBenard/RayleighBenard_cg_80x36_1_1_Output_1_1",
            "PCG_SSOR": "RayleighBenard/RayleighBenard_pcg_ssor_80x36_1_1_Output_1_1",
            "Eigen_CG": "RayleighBenard/RayleighBenard_eigen_cg_80x36_1_1_Output_1_1",
        },
    },
    {
        "label": "FT 200×100",
        "solvers": {
            "SOR":      "FluidTrap/FluidTrap_200x100_sor_1_1_Output_1_1",
            "CG":       "FluidTrap/FluidTrap_200x100_cg_1_1_Output_1_1",
            "PCG_SSOR": "FluidTrap/FluidTrap_200x100_pcg_1_1_Output_1_1",
            "Eigen_CG": "FluidTrap/FluidTrap_200x100_eigen_cg_1_1_Output_1_1",
        },
    },
    {
        "label": "LDC 100×100",
        "solvers": {
            "SOR":      "LidDrivenCavity50/LidDrivenCavity100_sor_1_1_Output_1_1",
            "CG":       "LidDrivenCavity50/LidDrivenCavity100_cg_1_1_Output_1_1",
            "PCG_SSOR": "LidDrivenCavity50/LidDrivenCavity100_pcg_ssor_1_1_Output_1_1",
            "Eigen_CG": "LidDrivenCavity50/LidDrivenCavity100_eigen_cg_1_1_Output_1_1",
        },
    },
]

SOLVERS = ["SOR", "CG", "PCG_SSOR", "Eigen_CG"]
COLORS  = {"SOR": "#f28e2b", "CG": "#4e79a7", "PCG_SSOR": "#59a14f", "Eigen_CG": "#e15759"}
LABELS  = {"SOR": "SOR", "CG": "CG (ours)", "PCG_SSOR": "PCG_SSOR (ours)", "Eigen_CG": "Eigen CG"}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_log(dirpath):
    p = os.path.join(BASE, dirpath, "sor_log.csv")
    if not os.path.isfile(p):
        return None
    try:
        df = pd.read_csv(p)
        if df.empty:
            return None
        df.columns = [c.strip() for c in df.columns]
        return df
    except Exception:
        return None

def load_walltime(dirpath):
    p = os.path.join(BASE, dirpath, "walltime.txt")
    if not os.path.isfile(p):
        return None
    return float(open(p).read().strip())

def iter_col(df):
    return next((c for c in df.columns if "iter" in c.lower()), None)

def gather():
    rows = []
    for case in CASES:
        for solver, dirpath in case["solvers"].items():
            log = load_log(dirpath)
            wt  = load_walltime(dirpath)
            if log is None or wt is None:
                rows.append(dict(case=case["label"], solver=solver,
                                 avg_iter=np.nan, walltime=np.nan,
                                 total_iter=np.nan, tpi=np.nan))
                continue
            col        = iter_col(log)
            avg_iter   = log[col].mean()
            total_iter = log[col].sum()
            tpi        = (wt / total_iter * 1000) if total_iter > 0 else np.nan
            rows.append(dict(case=case["label"], solver=solver,
                             avg_iter=avg_iter, walltime=wt,
                             total_iter=total_iter, tpi=tpi))
    return pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def grouped_bars(ax, df, metric, ylabel, title, fmt="{:.2f}", logy=False):
    case_labels = [c["label"] for c in CASES]
    x = np.arange(len(case_labels))
    n = len(SOLVERS)
    w = 0.18
    offsets = np.linspace(-(n-1)/2, (n-1)/2, n) * w

    for offset, solver in zip(offsets, SOLVERS):
        vals = [df.loc[(df.case == cl) & (df.solver == solver), metric].values
                for cl in case_labels]
        vals = [v[0] if len(v) else np.nan for v in vals]
        bars = ax.bar(x + offset, vals, w, label=LABELS[solver],
                      color=COLORS[solver])
        vmax = max((v for v in vals if not np.isnan(v)), default=1.0)
        for bar, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        v + vmax * 0.01,
                        fmt.format(v), ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(case_labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

df = gather()

print("\n=== Summary (large grids) ===")
print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

missing = df[df["walltime"].isna()]
if not missing.empty:
    print("\n⚠ Missing data (run not finished yet):")
    for _, r in missing.iterrows():
        print(f"  {r['case']} — {r['solver']}")

# Page 1: bar charts
fig1, axes = plt.subplots(1, 3, figsize=(17, 5))

grouped_bars(axes[0], df, "avg_iter",
             "Avg. iterations / timestep",
             "Avg. iterations per timestep\n(lower = better convergence)",
             fmt="{:.2f}")

grouped_bars(axes[1], df, "walltime",
             "Wall time [s]",
             "Total wall time",
             fmt="{:.0f}")

grouped_bars(axes[2], df, "tpi",
             "Time per iteration [ms]",
             "Time per iteration [ms]\n(lower = cheaper matvec)",
             fmt="{:.2f}")

fig1.suptitle("SOR vs CG vs PCG_SSOR vs Eigen_CG — large grids (2× small), 1×1, warm start",
              fontsize=12, fontweight="bold")
fig1.text(0.5, -0.02,
          "SOR: stationary (no warm start). CG/PCG/Eigen_CG: warm start from previous timestep.\n"
          "Domains: RB 80×36, FT 200×100, LDC 100×100 (each 2× the small-grid case).",
          ha="center", fontsize=8, color="gray")
plt.tight_layout()

# Page 2: per-timestep iteration traces
fig2, axes2 = plt.subplots(1, len(CASES), figsize=(17, 4))
for ax, case in zip(axes2, CASES):
    for solver, dirpath in case["solvers"].items():
        log = load_log(dirpath)
        if log is None:
            continue
        col = iter_col(log)
        ax.plot(log["t"], log[col],
                color=COLORS[solver], alpha=0.75, lw=0.8, label=LABELS[solver])
    ax.set_xlabel("Simulation time t")
    ax.set_ylabel("Iterations")
    ax.set_title(case["label"])
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
fig2.suptitle("Iterations per timestep over time — large grids (2× small), 1×1",
              fontsize=11, fontweight="bold")
plt.tight_layout()

# Page 3: speedup vs SOR
fig3, axes3 = plt.subplots(1, 2, figsize=(13, 5))
case_labels = [c["label"] for c in CASES]
x = np.arange(len(case_labels))
n_other = len(SOLVERS) - 1
w = 0.22
offsets = np.linspace(-(n_other-1)/2, (n_other-1)/2, n_other) * w

for ax, metric, ylabel, title in [
    (axes3[0], "walltime", "Speedup (SOR / solver walltime)", "Walltime speedup vs SOR"),
    (axes3[1], "avg_iter", "Iter. reduction (SOR / solver avg)", "Iteration reduction vs SOR"),
]:
    for offset, solver in zip(offsets, [s for s in SOLVERS if s != "SOR"]):
        speedups = []
        for cl in case_labels:
            sor_val = df.loc[(df.case == cl) & (df.solver == "SOR"), metric].values
            sol_val = df.loc[(df.case == cl) & (df.solver == solver), metric].values
            if (len(sor_val) and len(sol_val)
                    and not np.isnan(sol_val[0]) and sol_val[0] > 0
                    and not np.isnan(sor_val[0])):
                speedups.append(sor_val[0] / sol_val[0])
            else:
                speedups.append(np.nan)
        bars = ax.bar(x + offset, speedups, w, label=LABELS[solver],
                      color=COLORS[solver])
        for bar, v in zip(bars, speedups):
            if not np.isnan(v):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        v + 0.02, f"{v:.2f}×",
                        ha="center", va="bottom", fontsize=8)

    ax.axhline(1.0, color="gray", lw=1, ls="--", label="SOR baseline")
    ax.set_xticks(x)
    ax.set_xticklabels(case_labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)

fig3.suptitle("Speedup relative to SOR — large grids (2× small), 1×1", fontsize=12, fontweight="bold")
plt.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
os.makedirs(PICTURES, exist_ok=True)
out = os.path.join(PICTURES, "fig76_all_solvers_large.pdf")
with PdfPages(out) as pdf:
    pdf.savefig(fig1, bbox_inches="tight")
    pdf.savefig(fig2, bbox_inches="tight")
    pdf.savefig(fig3, bbox_inches="tight")
plt.close("all")

print(f"\nSaved: {out}")
print("Done.")
