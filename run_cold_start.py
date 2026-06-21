#!/usr/bin/env python3
"""
Run all cold-start serial (1x1) cases and generate comparison plots.

Cases:
  - LidDrivenCavity50  50x50   (t_end=50)
  - RayleighBenard     40x18   (t_end=10000)
  - FluidTrap s100     100x50  (t_end=2000)

Solvers: SOR, CG, PCG-SSOR  (each 1x1 only)

Usage:
  python3 run_cold_start.py          # run simulations + plots
  python3 run_cold_start.py --plot   # plots only (skip simulations)
"""

import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "figure.dpi": 150})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BIN         = os.path.join(SCRIPT_DIR, "build", "fluidchen")
CASES_ROOT  = os.path.join(SCRIPT_DIR, "example_cases")
PICTURES    = ("/Users/thomaslutz/Documents/Master CSE/2.Semester/CFD Lab"
               "/Project Propsosal/Output_analysis_PCG_SSOR/pictures")

CASES = [
    {"key": "LidDrivenCavity50", "label": "LDC 50×50",          "prefix": "LidDrivenCavity50"},
    {"key": "RayleighBenard",    "label": "RayleighBenard 40×18","prefix": "RayleighBenard"},
    {"key": "FluidTrap",         "label": "FluidTrap 100×50",    "prefix": "FluidTrap_s100"},
]
SOLVERS = [
    {"key": "sor", "label": "SOR",      "color": "#e05252"},
    {"key": "cg",  "label": "CG",       "color": "#5276e0"},
    {"key": "pcg", "label": "PCG-SSOR", "color": "#52a85e"},
]

FIG_START = 63

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def dat_path(case, solver):
    prefix = case["prefix"]
    d = os.path.join(CASES_ROOT, case["key"])
    return os.path.join(d, f"{prefix}_{solver['key']}_cold_1_1.dat")


def output_dir(case, solver):
    prefix = case["prefix"]
    pattern = os.path.join(CASES_ROOT, case["key"],
                           f"{prefix}_{solver['key']}_cold_1_1_Output_1_1")
    hits = glob.glob(pattern)
    return hits[0] if hits else pattern  # return expected path even if missing


def read_walltime(odir):
    path = os.path.join(odir, "walltime.txt")
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        m = re.search(r"([\d.]+)", f.read())
        return float(m.group(1)) if m else None


def read_sor_log(odir):
    path = os.path.join(odir, "sor_log.csv")
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def fig_path(num, name):
    return os.path.join(PICTURES, f"fig{num:02d}_{name}.pdf")


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------

def run_simulations():
    if not os.path.isfile(BIN):
        sys.exit(f"ERROR: binary not found at {BIN} — run cmake --build build first")

    for case in CASES:
        print(f"\n=== {case['label']} ===")
        for solver in SOLVERS:
            dat = dat_path(case, solver)
            odir = output_dir(case, solver)
            wt = read_walltime(odir)
            if wt is not None:
                print(f"  [skip] {solver['label']} — already done ({wt:.1f}s)")
                continue
            print(f"  [run]  {solver['label']} ...", flush=True)
            result = subprocess.run(
                ["mpirun", "-n", "1", BIN, dat],
                capture_output=True, text=True
            )
            wt = read_walltime(odir)
            status = f"{wt:.1f}s" if wt else "ERROR"
            print(f"         → {status}")
            if result.returncode != 0:
                print(f"         STDERR: {result.stderr[:300]}")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_walltime(fig_num):
    """Grouped bar: SOR / CG / PCG side by side per case."""
    fig, axes = plt.subplots(1, len(CASES), figsize=(13, 4.5))
    for ax, case in zip(axes, CASES):
        walltimes = []
        labels = []
        colors = []
        for solver in SOLVERS:
            wt = read_walltime(output_dir(case, solver))
            if wt is not None:
                walltimes.append(wt)
                labels.append(solver["label"])
                colors.append(solver["color"])

        x = np.arange(len(walltimes))
        bars = ax.bar(x, walltimes, color=colors, width=0.5)
        for bar, wt in zip(bars, walltimes):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    wt + max(walltimes) * 0.01,
                    f"{wt:.1f}s", ha="center", va="bottom", fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Wall time [s]")
        ax.set_title(case["label"])

    fig.suptitle("Cold-start wall time — serial (1×1)", fontsize=12)
    plt.tight_layout()
    out = fig_path(fig_num, "cold_walltime")
    plt.savefig(out)
    plt.close()
    print(f"  fig{fig_num:02d}: {os.path.basename(out)}")
    return fig_num + 1


def plot_iter_timeseries(fig_num):
    """Iteration count over physical time, one figure per case."""
    for case in CASES:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        has_data = False
        for solver in SOLVERS:
            df = read_sor_log(output_dir(case, solver))
            if df is None or df.empty:
                continue
            time_col = next((c for c in df.columns if "t" == c.lower()), None)
            iter_col = next((c for c in df.columns if "iter" in c.lower()), None)
            if time_col is None or iter_col is None:
                continue
            ax.plot(df[time_col], df[iter_col],
                    label=solver["label"], color=solver["color"],
                    linewidth=0.9, alpha=0.85)
            has_data = True

        if not has_data:
            plt.close()
            continue

        ax.set_xlabel("Physical time [s]")
        ax.set_ylabel("Pressure solver iterations")
        ax.set_title(f"Cold-start iterations over time — {case['label']}")
        ax.legend()
        plt.tight_layout()
        out = fig_path(fig_num, f"cold_iters_ts_{case['key'].lower()}")
        plt.savefig(out)
        plt.close()
        print(f"  fig{fig_num:02d}: {os.path.basename(out)}")
        fig_num += 1

    return fig_num


def plot_avg_iters(fig_num):
    """Average iterations bar chart, all cases side by side."""
    fig, axes = plt.subplots(1, len(CASES), figsize=(13, 4.5))
    for ax, case in zip(axes, CASES):
        avgs = []
        labels = []
        colors = []
        for solver in SOLVERS:
            df = read_sor_log(output_dir(case, solver))
            if df is None or df.empty:
                continue
            iter_col = next((c for c in df.columns if "iter" in c.lower()), None)
            if iter_col is None:
                continue
            avgs.append(df[iter_col].mean())
            labels.append(solver["label"])
            colors.append(solver["color"])

        x = np.arange(len(avgs))
        bars = ax.bar(x, avgs, color=colors, width=0.5)
        for bar, a in zip(bars, avgs):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    a + max(avgs) * 0.01,
                    f"{a:.1f}", ha="center", va="bottom", fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Avg. iterations / timestep")
        ax.set_title(case["label"])

    fig.suptitle("Cold-start avg. iterations — serial (1×1)", fontsize=12)
    plt.tight_layout()
    out = fig_path(fig_num, "cold_avg_iters")
    plt.savefig(out)
    plt.close()
    print(f"  fig{fig_num:02d}: {os.path.basename(out)}")
    return fig_num + 1


def generate_plots():
    os.makedirs(PICTURES, exist_ok=True)
    fig_num = FIG_START
    print("\n=== Generating plots ===")
    fig_num = plot_walltime(fig_num)
    fig_num = plot_iter_timeseries(fig_num)
    fig_num = plot_avg_iters(fig_num)
    print(f"\nDone. Next available figure: fig{fig_num:02d}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot", action="store_true",
                        help="Only generate plots, skip simulations")
    args = parser.parse_args()

    if not args.plot:
        run_simulations()

    generate_plots()
