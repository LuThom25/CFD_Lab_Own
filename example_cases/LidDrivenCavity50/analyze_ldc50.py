"""
Strong Scaling Analysis — 50×50 Lid-Driven Cavity (small-grid comparison)
Compares SOR_STANDARD, CG_STANDARD and PCG_SSOR
Decompositions: 1×1, 2×2

Run after all simulations are complete:
    python3 analyze_ldc50.py
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

BASE      = os.path.dirname(os.path.abspath(__file__))
PICTURES  = os.path.join(
    os.path.expanduser("~"),
    "Documents/Master CSE/2.Semester/CFD Lab/Project Propsosal/Output_analysis_PCG_SSOR/pictures"
)
os.makedirs(PICTURES, exist_ok=True)

SOR_CONFIGS = [
    ("LidDrivenCavity50_sor_1_1", "SOR (1×1)", 1, "-",  "#1E88E5"),
    ("LidDrivenCavity50_sor_2_2", "SOR (2×2)", 4, "--", "#1565C0"),
]
CG_CONFIGS = [
    ("LidDrivenCavity50_cg_1_1",  "CG (1×1)",  1, "-",  "#43A047"),
    ("LidDrivenCavity50_cg_2_2",  "CG (2×2)",  4, "--", "#1B5E20"),
]
PCG_CONFIGS = [
    ("LidDrivenCavity50_pcg_1_1", "PCG (1×1)", 1, "-",  "#8E24AA"),
    ("LidDrivenCavity50_pcg_2_2", "PCG (2×2)", 4, "--", "#4A148C"),
]
ALL_CONFIGS = SOR_CONFIGS + CG_CONFIGS + PCG_CONFIGS

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def output_dir(stem):
    parts = stem.rsplit("_", 2)
    ip, jp = parts[1], parts[2]
    return os.path.join(BASE, f"{stem}_Output_{ip}_{jp}")

def load_csv(stem):
    path = os.path.join(output_dir(stem), "sor_log.csv")
    return pd.read_csv(path) if os.path.exists(path) else None

def load_walltime(stem):
    path = os.path.join(output_dir(stem), "walltime.txt")
    return float(open(path).read().strip()) if os.path.exists(path) else None

def save(fig, name):
    for dest in [BASE, PICTURES]:
        fig.savefig(os.path.join(dest, name), dpi=150, bbox_inches="tight")
    print(f"  Saved: {name}")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

print("Loading data...")
dfs, times = {}, {}
for stem, label, *_ in ALL_CONFIGS:
    dfs[stem]   = load_csv(stem)
    times[stem] = load_walltime(stem)
    df, wt = dfs[stem], times[stem]
    status = f"{len(df)} steps, {wt:.2f}s" if df is not None and wt else \
             f"{len(df)} steps"             if df is not None       else "not found"
    print(f"  {label:20s}: {status}")
print()

# ---------------------------------------------------------------------------
# Plot 1: Iterations per step — all solvers 1×1
# ---------------------------------------------------------------------------

print("Generating iteration plots...")

def iters_plot(configs, title, outfile):
    if not any(dfs.get(s) is not None for s, *_ in configs):
        print(f"  [skip] {title}"); return
    fig, ax = plt.subplots(figsize=(9, 4))
    for stem, label, _, ls, col in configs:
        if dfs.get(stem) is None: continue
        df = dfs[stem]
        ax.plot(df["timestep"], df["sor_iters"],
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

iters_plot(SOR_CONFIGS, "SOR iterations/step — 50×50 Lid-Driven Cavity",     "plot_ldc50_sor_iters.pdf")
iters_plot(CG_CONFIGS,  "CG iterations/step — 50×50 Lid-Driven Cavity",      "plot_ldc50_cg_iters.pdf")
iters_plot(PCG_CONFIGS, "PCG_SSOR iterations/step — 50×50 Lid-Driven Cavity","plot_ldc50_pcg_iters.pdf")

iters_plot(
    [("LidDrivenCavity50_sor_1_1", "SOR (1×1)", 1, "-", "#1E88E5"),
     ("LidDrivenCavity50_cg_1_1",  "CG (1×1)",  1, "-", "#43A047"),
     ("LidDrivenCavity50_pcg_1_1", "PCG (1×1)", 1, "-", "#8E24AA")],
    "SOR vs CG vs PCG — iterations/step (1×1, 50×50 Lid-Driven Cavity)",
    "plot_ldc50_all_solvers_1_1.pdf")

iters_plot(
    [("LidDrivenCavity50_sor_2_2", "SOR (2×2)", 4, "-", "#1E88E5"),
     ("LidDrivenCavity50_cg_2_2",  "CG (2×2)",  4, "-", "#43A047"),
     ("LidDrivenCavity50_pcg_2_2", "PCG (2×2)", 4, "-", "#8E24AA")],
    "SOR vs CG vs PCG — iterations/step (2×2, 50×50 Lid-Driven Cavity)",
    "plot_ldc50_all_solvers_2_2.pdf")

# ---------------------------------------------------------------------------
# Plot 2: Wall time bar chart
# ---------------------------------------------------------------------------

print("Generating wall time plot...")

decomp_labels = ["1×1", "2×2"]
sor_vals  = [times.get(s) for s, *_ in SOR_CONFIGS]
cg_vals   = [times.get(s) for s, *_ in CG_CONFIGS]
pcg_vals  = [times.get(s) for s, *_ in PCG_CONFIGS]
all_vals  = [v for v in sor_vals + cg_vals + pcg_vals if v]

if all_vals:
    fig, ax = plt.subplots(figsize=(8, 5))
    w, n, mx = 0.25, len(SOR_CONFIGS), max(all_vals)
    xs = list(range(n))

    def bar_group(vals, offset, color, label):
        bars = ax.bar([x + offset for x in xs],
                      [v if v else 0 for v in vals],
                      width=w, color=color, label=label, edgecolor="white", alpha=0.85)
        for bar, val in zip(bars, vals):
            if val:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + mx*0.01,
                        f"{val:.2f}s", ha="center", va="bottom", fontsize=9)

    bar_group(sor_vals, -w,  "#1E88E5", "SOR_STANDARD")
    bar_group(cg_vals,   0,  "#43A047", "CG_STANDARD")
    bar_group(pcg_vals,  w,  "#8E24AA", "PCG_SSOR")

    ax.set_xticks(xs)
    ax.set_xticklabels(decomp_labels, fontsize=12)
    ax.set_ylabel("Wall time (s)", fontsize=13)
    ax.set_title("Wall time: SOR vs CG vs PCG_SSOR — 50×50 Lid-Driven Cavity",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)
    ax.set_ylim(top=mx * 1.18)
    plt.tight_layout()
    save(fig, "plot_ldc50_walltime.pdf")
    plt.close()

# ---------------------------------------------------------------------------
# Plot 3: Speedup
# ---------------------------------------------------------------------------

print("Generating speedup plot...")

baseline_sor = times.get("LidDrivenCavity50_sor_1_1")
baseline_cg  = times.get("LidDrivenCavity50_cg_1_1")
baseline_pcg = times.get("LidDrivenCavity50_pcg_1_1")

if any([baseline_sor, baseline_cg, baseline_pcg]):
    fig, ax = plt.subplots(figsize=(7, 4))

    def speedup_line(vals, baseline, color, marker, label):
        if not baseline: return
        su  = [baseline / v if v else None for v in vals]
        xs_ = [i for i, v in enumerate(su) if v is not None]
        ys_ = [v for v in su if v is not None]
        if not xs_: return
        ax.plot(xs_, ys_, marker=marker, color=color, label=label, linewidth=2)
        for x, y in zip(xs_, ys_):
            ax.annotate(f"{y:.2f}×", (x, y), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=10, color=color)

    speedup_line(sor_vals,  baseline_sor, "#1E88E5", "s", "SOR Speedup")
    speedup_line(cg_vals,   baseline_cg,  "#43A047", "o", "CG Speedup")
    speedup_line(pcg_vals,  baseline_pcg, "#8E24AA", "^", "PCG Speedup")

    ax.set_xticks(range(n))
    ax.set_xticklabels(decomp_labels, fontsize=12)
    ax.set_ylabel("Speedup (vs. 1×1 baseline)", fontsize=13)
    ax.set_title("Strong Scaling Speedup — 50×50 Lid-Driven Cavity",
                 fontsize=14, fontweight="bold")
    ax.axhline(y=1, color="gray", linestyle=":", alpha=0.5)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    save(fig, "plot_ldc50_speedup.pdf")
    plt.close()

# ---------------------------------------------------------------------------
# Plot 4: Grid size comparison — LDC 50×50 vs 150×150 (1×1 only)
# ---------------------------------------------------------------------------

print("Generating grid size comparison plot...")

ldc300_base = os.path.join(os.path.dirname(BASE), "LidDrivenCavity300")

def load_wt_external(folder, stem):
    ip, jp = stem.rsplit("_", 2)[1], stem.rsplit("_", 2)[2]
    path = os.path.join(folder, f"{stem}_Output_{ip}_{jp}", "walltime.txt")
    return float(open(path).read().strip()) if os.path.exists(path) else None

ldc300_walltimes = {
    "SOR 150×150": load_wt_external(ldc300_base, "LidDrivenCavity300_sor_1_1"),
    "CG 150×150":  load_wt_external(ldc300_base, "LidDrivenCavity300_cg_1_1"),
    "PCG 150×150": load_wt_external(ldc300_base, "LidDrivenCavity300_pcg_1_1"),
}
ldc50_walltimes = {
    "SOR 50×50":  times.get("LidDrivenCavity50_sor_1_1"),
    "CG 50×50":   times.get("LidDrivenCavity50_cg_1_1"),
    "PCG 50×50":  times.get("LidDrivenCavity50_pcg_1_1"),
}

all_wt = [v for v in list(ldc300_walltimes.values()) + list(ldc50_walltimes.values()) if v]
if all_wt:
    fig, ax = plt.subplots(figsize=(9, 5))
    solvers = ["SOR", "CG", "PCG"]
    colors  = ["#1E88E5", "#43A047", "#8E24AA"]
    xs = list(range(len(solvers)))
    w  = 0.35
    mx = max(all_wt)

    vals_300 = [ldc300_walltimes.get(f"{s} 150×150") for s in solvers]
    vals_50  = [ldc50_walltimes.get(f"{s} 50×50")    for s in solvers]

    bars1 = ax.bar([x - w/2 for x in xs], [v or 0 for v in vals_300],
                   width=w, color=colors, alpha=0.6, edgecolor="white", label="150×150")
    bars2 = ax.bar([x + w/2 for x in xs], [v or 0 for v in vals_50],
                   width=w, color=colors, alpha=1.0, edgecolor="white", label="50×50",
                   hatch="//")

    for bar, val in zip(list(bars1) + list(bars2),
                        [v or 0 for v in vals_300] + [v or 0 for v in vals_50]):
        if val:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + mx*0.01,
                    f"{val:.1f}s", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(xs)
    ax.set_xticklabels(solvers, fontsize=13)
    ax.set_ylabel("Wall time (s)", fontsize=13)
    ax.set_title("Wall time comparison: 50×50 vs 150×150 Lid-Driven Cavity (1×1)",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)
    ax.set_ylim(top=mx * 1.2)
    plt.tight_layout()
    save(fig, "plot_ldc50_vs_ldc300.pdf")
    plt.close()

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

print("\n=== Summary Table — 50×50 Lid-Driven Cavity ===")
print(f"{'Config':<24} {'Wall [s]':>10} {'Speedup':>10} {'Avg iters':>12} {'Max iters':>12}")
print("-" * 70)
for configs, baseline in [(SOR_CONFIGS, baseline_sor),
                           (CG_CONFIGS,  baseline_cg),
                           (PCG_CONFIGS, baseline_pcg)]:
    for stem, label, *_ in configs:
        wt = times.get(stem)
        df = dfs.get(stem)
        su     = f"{baseline/wt:.2f}×" if wt and baseline else "—"
        wt_str = f"{wt:.2f}"           if wt              else "—"
        avg    = f"{df['sor_iters'].mean():.2f}" if df is not None else "—"
        mx_str = f"{df['sor_iters'].max()}"      if df is not None else "—"
        print(f"  {label:<22} {wt_str:>10} {su:>10} {avg:>12} {mx_str:>12}")
    print()

print("Done.")
