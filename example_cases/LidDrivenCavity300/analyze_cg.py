"""
Strong Scaling Analysis — 150×150 Lid-Driven Cavity
Compares SOR_STANDARD, CG_STANDARD and PCG_SSOR
Decompositions: serial, 1×1, 2×2, 1×4

Run after all simulations are complete:
    python3 analyze_cg.py
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

BASE      = os.path.dirname(os.path.abspath(__file__))
DOWNLOADS = os.path.expanduser("~/Downloads")

# (stem, label, n_procs, linestyle, color)
SOR_CONFIGS = [
    ("LidDrivenCavity300_sor_serial", "SOR Serial",    1, "-",  "#90CAF9"),
    ("LidDrivenCavity300_sor_1_1",    "SOR MPI (1×1)", 1, "--", "#42A5F5"),
    ("LidDrivenCavity300_sor_2_2",    "SOR MPI (2×2)", 4, "-",  "#1E88E5"),
    ("LidDrivenCavity300_sor_1_4",    "SOR MPI (1×4)", 4, "--", "#1565C0"),
]

CG_CONFIGS = [
    ("LidDrivenCavity300_cg_serial", "CG Serial",    1, "-",  "#A5D6A7"),
    ("LidDrivenCavity300_cg_1_1",    "CG MPI (1×1)", 1, "--", "#66BB6A"),
    ("LidDrivenCavity300_cg_2_2",    "CG MPI (2×2)", 4, "-",  "#43A047"),
    ("LidDrivenCavity300_cg_1_4",    "CG MPI (1×4)", 4, "--", "#1B5E20"),
]

PCG_CONFIGS = [
    ("LidDrivenCavity300_pcg_serial", "PCG Serial",    1, "-",  "#CE93D8"),
    ("LidDrivenCavity300_pcg_1_1",    "PCG MPI (1×1)", 1, "--", "#AB47BC"),
    ("LidDrivenCavity300_pcg_2_2",    "PCG MPI (2×2)", 4, "-",  "#8E24AA"),
    ("LidDrivenCavity300_pcg_1_4",    "PCG MPI (1×4)", 4, "--", "#4A148C"),
]

ALL_CONFIGS = SOR_CONFIGS + CG_CONFIGS + PCG_CONFIGS

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def iproc_jproc(stem):
    parts = stem.rsplit("_", 2)
    if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
        return int(parts[1]), int(parts[2])
    return 1, 1

def output_dir(stem):
    ip, jp = iproc_jproc(stem)
    return os.path.join(BASE, f"{stem}_Output_{ip}_{jp}")

def load_csv(stem):
    csv = os.path.join(output_dir(stem), "sor_log.csv")
    return pd.read_csv(csv) if os.path.exists(csv) else None

def load_walltime(stem):
    wt = os.path.join(output_dir(stem), "walltime.txt")
    return float(open(wt).read().strip()) if os.path.exists(wt) else None

def save(fig, name):
    for dest in [BASE, DOWNLOADS]:
        fig.savefig(os.path.join(dest, name), dpi=150)
    print(f"  Saved: {name}")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

print("Loading data...")
dfs   = {}
times = {}
for stem, label, *_ in ALL_CONFIGS:
    dfs[stem]   = load_csv(stem)
    times[stem] = load_walltime(stem)
    df, wt = dfs[stem], times[stem]
    status = f"{len(df)} steps, {wt:.1f}s" if df is not None and wt else \
             f"{len(df)} steps"             if df is not None       else "not found"
    print(f"  {label:25s}: {status}")
print()

# ---------------------------------------------------------------------------
# Plot helper: iterations per step
# ---------------------------------------------------------------------------

def iters_plot(configs, title, outfile):
    if not any(dfs.get(s) is not None for s, *_ in configs):
        print(f"  [skip] {title}")
        return
    fig, ax = plt.subplots(figsize=(9, 4))
    for stem, label, _, ls, col in configs:
        if dfs.get(stem) is None:
            continue
        ax.plot(dfs[stem]["timestep"], dfs[stem]["sor_iters"],
                label=label, linestyle=ls, color=col, linewidth=1.4, alpha=0.9)
    ax.set_xlabel("Timestep", fontsize=13)
    ax.set_ylabel("Solver Iterations", fontsize=13)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(fontsize=10, framealpha=0.8)
    ax.tick_params(labelsize=11)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    save(fig, outfile)
    plt.close()

# ---------------------------------------------------------------------------
# Plot 1: Per-solver iteration plots + comparisons
# ---------------------------------------------------------------------------

print("Generating iteration plots...")
iters_plot(SOR_CONFIGS, "SOR iterations/step — 150×150 Lid-Driven Cavity",      "plot_ldc_sor_iters.pdf")
iters_plot(CG_CONFIGS,  "CG iterations/step — 150×150 Lid-Driven Cavity",        "plot_ldc_cg_iters.pdf")
iters_plot(PCG_CONFIGS, "PCG_SSOR iterations/step — 150×150 Lid-Driven Cavity",  "plot_ldc_pcg_iters.pdf")

iters_plot(
    [("LidDrivenCavity300_sor_serial", "SOR Serial",  1, "-", "#1E88E5"),
     ("LidDrivenCavity300_cg_serial",  "CG Serial",   1, "-", "#43A047"),
     ("LidDrivenCavity300_pcg_serial", "PCG Serial",  1, "-", "#8E24AA")],
    "SOR vs CG vs PCG — iterations/step (Serial, 150×150)",
    "plot_ldc_all_solvers_serial.pdf")

iters_plot(
    [("LidDrivenCavity300_sor_1_1", "SOR (1×1)", 1, "-", "#1E88E5"),
     ("LidDrivenCavity300_cg_1_1",  "CG (1×1)",  1, "-", "#43A047"),
     ("LidDrivenCavity300_pcg_1_1", "PCG (1×1)", 1, "-", "#8E24AA")],
    "SOR vs CG vs PCG — iterations/step (1×1, 150×150)",
    "plot_ldc_all_solvers_1_1.pdf")

iters_plot(
    [("LidDrivenCavity300_sor_2_2", "SOR (2×2)", 4, "-", "#1E88E5"),
     ("LidDrivenCavity300_cg_2_2",  "CG (2×2)",  4, "-", "#43A047"),
     ("LidDrivenCavity300_pcg_2_2", "PCG (2×2)", 4, "-", "#8E24AA")],
    "SOR vs CG vs PCG — iterations/step (2×2, 150×150)",
    "plot_ldc_all_solvers_2_2.pdf")

iters_plot(
    [("LidDrivenCavity300_sor_1_4", "SOR (1×4)", 4, "-", "#1E88E5"),
     ("LidDrivenCavity300_cg_1_4",  "CG (1×4)",  4, "-", "#43A047"),
     ("LidDrivenCavity300_pcg_1_4", "PCG (1×4)", 4, "-", "#8E24AA")],
    "SOR vs CG vs PCG — iterations/step (1×4, 150×150)",
    "plot_ldc_all_solvers_1_4.pdf")

iters_plot(
    [("LidDrivenCavity300_sor_2_2", "SOR (2×2)", 4, "-",  "#1E88E5"),
     ("LidDrivenCavity300_sor_1_4", "SOR (1×4)", 4, "--", "#90CAF9"),
     ("LidDrivenCavity300_cg_2_2",  "CG (2×2)",  4, "-",  "#43A047"),
     ("LidDrivenCavity300_cg_1_4",  "CG (1×4)",  4, "--", "#A5D6A7"),
     ("LidDrivenCavity300_pcg_2_2", "PCG (2×2)", 4, "-",  "#8E24AA"),
     ("LidDrivenCavity300_pcg_1_4", "PCG (1×4)", 4, "--", "#CE93D8")],
    "SOR vs CG vs PCG — iterations/step (4 Ranks, 150×150)",
    "plot_ldc_all_solvers_4ranks.pdf")

# ---------------------------------------------------------------------------
# Plot 2: Wall time — all 3 solvers grouped
# ---------------------------------------------------------------------------

print("Generating wall time plot...")

decomp_labels = ["Serial", "1×1", "2×2", "1×4"]
sor_vals  = [times.get(s) for s, *_ in SOR_CONFIGS]
cg_vals   = [times.get(s) for s, *_ in CG_CONFIGS]
pcg_vals  = [times.get(s) for s, *_ in PCG_CONFIGS]
all_vals  = [v for v in sor_vals + cg_vals + pcg_vals if v]

if all_vals:
    fig, ax = plt.subplots(figsize=(12, 5))
    w, n, mx = 0.25, len(SOR_CONFIGS), max(all_vals)
    xs = list(range(n))

    def bar_group(vals, offset, color, label):
        bars = ax.bar([x + offset for x in xs],
                      [v if v else 0 for v in vals],
                      width=w, color=color, label=label, edgecolor="white", alpha=0.85)
        for bar, val in zip(bars, vals):
            if val:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + mx*0.01,
                        f"{val:.0f}s", ha="center", va="bottom", fontsize=9)

    bar_group(sor_vals, -w,  "#1E88E5", "SOR_STANDARD")
    bar_group(cg_vals,   0,  "#43A047", "CG_STANDARD")
    bar_group(pcg_vals,  w,  "#8E24AA", "PCG_SSOR")

    ax.set_xticks(xs)
    ax.set_xticklabels(decomp_labels, fontsize=12)
    ax.set_ylabel("Wall time (s)", fontsize=13)
    ax.set_title("Wall time: SOR vs CG vs PCG_SSOR — 150×150 Lid-Driven Cavity",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)
    ax.set_ylim(top=mx * 1.18)
    plt.tight_layout()
    save(fig, "plot_ldc_walltime.pdf")
    plt.close()

# ---------------------------------------------------------------------------
# Plot 3: Speedup
# ---------------------------------------------------------------------------

print("Generating speedup plot...")

baseline_sor = times.get("LidDrivenCavity300_sor_serial") or times.get("LidDrivenCavity300_sor_1_1")
baseline_cg  = times.get("LidDrivenCavity300_cg_serial")  or times.get("LidDrivenCavity300_cg_1_1")
baseline_pcg = times.get("LidDrivenCavity300_pcg_serial") or times.get("LidDrivenCavity300_pcg_1_1")

if any([baseline_sor, baseline_cg, baseline_pcg]):
    fig, ax = plt.subplots(figsize=(9, 4))

    def speedup_line(vals, baseline, color, marker, label):
        if not baseline:
            return
        su = [baseline / v if v else None for v in vals]
        xs_ = [i for i, v in enumerate(su) if v is not None]
        ys_ = [v for v in su if v is not None]
        if not xs_:
            return
        ax.plot(xs_, ys_, marker=marker, color=color, label=label, linewidth=2)
        for x, y in zip(xs_, ys_):
            ax.annotate(f"{y:.2f}×", (x, y), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=10, color=color)

    speedup_line(sor_vals,  baseline_sor, "#1E88E5", "s", "SOR Speedup")
    speedup_line(cg_vals,   baseline_cg,  "#43A047", "o", "CG Speedup")
    speedup_line(pcg_vals,  baseline_pcg, "#8E24AA", "^", "PCG Speedup")

    ax.set_xticks(range(n))
    ax.set_xticklabels(decomp_labels, fontsize=12)
    ax.set_ylabel("Speedup (vs. Serial baseline)", fontsize=13)
    ax.set_title("Strong Scaling Speedup — 150×150 Lid-Driven Cavity",
                 fontsize=14, fontweight="bold")
    ax.axhline(y=1, color="gray", linestyle=":", alpha=0.5)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    save(fig, "plot_ldc_speedup.pdf")
    plt.close()

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

print("\n=== Summary Table — 150×150 Lid-Driven Cavity ===")
print(f"{'Config':<28} {'Wall [s]':>10} {'Speedup':>10} {'Avg iters':>12} {'Max iters':>12}")
print("-" * 75)
for configs, baseline in [(SOR_CONFIGS, baseline_sor), (CG_CONFIGS, baseline_cg), (PCG_CONFIGS, baseline_pcg)]:
    for stem, label, *_ in configs:
        wt = times.get(stem)
        df = dfs.get(stem)
        su     = f"{baseline/wt:.2f}×" if wt and baseline else "—"
        wt_str = f"{wt:.1f}"           if wt              else "—"
        avg    = f"{df['sor_iters'].mean():.1f}" if df is not None else "—"
        mx_str = f"{df['sor_iters'].max()}"      if df is not None else "—"
        print(f"  {label:<26} {wt_str:>10} {su:>10} {avg:>12} {mx_str:>12}")
    print()

print("Done.")
