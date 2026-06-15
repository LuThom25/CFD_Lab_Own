"""
Strong Scaling Analysis — 150×150 Lid-Driven Cavity
Vergleich SOR_STANDARD vs CG_STANDARD
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
    ("LidDrivenCavity300_sor_1_1",    "SOR MPI (1×1)", 1, "--", "#FFCC02"),
    ("LidDrivenCavity300_sor_2_2",    "SOR MPI (2×2)", 4, "-",  "#A5D6A7"),
    ("LidDrivenCavity300_sor_1_4",    "SOR MPI (1×4)", 4, "--", "#EF9A9A"),
]

CG_CONFIGS = [
    ("LidDrivenCavity300_cg_serial", "CG Serial",    1, "-",  "#2196F3"),
    ("LidDrivenCavity300_cg_1_1",    "CG MPI (1×1)", 1, "--", "#FF9800"),
    ("LidDrivenCavity300_cg_2_2",    "CG MPI (2×2)", 4, "-",  "#4CAF50"),
    ("LidDrivenCavity300_cg_1_4",    "CG MPI (1×4)", 4, "--", "#F44336"),
]

ALL_CONFIGS = SOR_CONFIGS + CG_CONFIGS

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
    print(f"  {label:22s}: {status}")
print()

# ---------------------------------------------------------------------------
# Plot 1: Solver iterations/step — SOR vs CG (serial only, clear comparison)
# ---------------------------------------------------------------------------

def iters_plot(configs, title, outfile):
    avail = [s for s, *_ in configs if dfs.get(s) is not None]
    if not avail:
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

print("Generating iteration plots...")
iters_plot(SOR_CONFIGS,
           "SOR iterations/step — 150×150 Lid-Driven Cavity",
           "plot_sor_iters.pdf")
iters_plot(CG_CONFIGS,
           "CG iterations/step — 150×150 Lid-Driven Cavity",
           "plot_cg_iters.pdf")

# Serial SOR vs CG direct comparison
serial_both = [
    ("LidDrivenCavity300_sor_serial", "SOR Serial", 1, "-",  "#2196F3"),
    ("LidDrivenCavity300_cg_serial",  "CG Serial",  1, "-",  "#4CAF50"),
]
iters_plot(serial_both,
           "SOR vs CG iterations/step — Serial (150×150)",
           "plot_sor_vs_cg_iters_serial.pdf")

# 4-rank comparison
ranks4_both = [
    ("LidDrivenCavity300_sor_2_2", "SOR (2×2)", 4, "-",  "#2196F3"),
    ("LidDrivenCavity300_sor_1_4", "SOR (1×4)", 4, "--", "#FF9800"),
    ("LidDrivenCavity300_cg_2_2",  "CG (2×2)",  4, "-",  "#4CAF50"),
    ("LidDrivenCavity300_cg_1_4",  "CG (1×4)",  4, "--", "#F44336"),
]
iters_plot(ranks4_both,
           "SOR vs CG iterations/step — 4 Ranks (150×150)",
           "plot_sor_vs_cg_iters_4ranks.pdf")

# ---------------------------------------------------------------------------
# Plot 2: Wall time SOR vs CG
# ---------------------------------------------------------------------------

print("Generating wall time plot...")

sor_avail = [(s, l, c) for s, l, _, _, c in SOR_CONFIGS if times.get(s)]
cg_avail  = [(s, l, c) for s, l, _, _, c in CG_CONFIGS  if times.get(s)]

if sor_avail or cg_avail:
    fig, ax = plt.subplots(figsize=(10, 5))
    w = 0.35
    n = max(len(sor_avail), len(cg_avail))

    if sor_avail:
        vals = [times[s] for s, *_ in sor_avail]
        bars = ax.bar([i - w/2 for i in range(len(vals))], vals,
                      width=w, color="#2196F3", label="SOR_STANDARD", edgecolor="white", alpha=0.85)
        mx = max(times[s] for s, *_ in (sor_avail + cg_avail) if times.get(s))
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + mx*0.01,
                    f"{val:.0f}s", ha="center", va="bottom", fontsize=10)

    if cg_avail:
        vals = [times[s] for s, *_ in cg_avail]
        bars = ax.bar([i + w/2 for i in range(len(vals))], vals,
                      width=w, color="#4CAF50", label="CG_STANDARD", edgecolor="white", alpha=0.85)
        mx = max(times[s] for s, *_ in (sor_avail + cg_avail) if times.get(s))
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + mx*0.01,
                    f"{val:.0f}s", ha="center", va="bottom", fontsize=10)

    tick_lbls = [l.replace("SOR ","").replace("CG ","") for _, l, _ in (sor_avail or cg_avail)]
    ax.set_xticks(range(n))
    ax.set_xticklabels(tick_lbls, fontsize=12)
    ax.set_ylabel("Wall time (s)", fontsize=13)
    ax.set_title("Wall time: SOR vs CG — 150×150 Lid-Driven Cavity", fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)
    all_vals = [times[s] for s, *_ in sor_avail + cg_avail]
    ax.set_ylim(top=max(all_vals) * 1.15)
    plt.tight_layout()
    save(fig, "plot_walltime.pdf")
    plt.close()

# ---------------------------------------------------------------------------
# Plot 3: Speedup
# ---------------------------------------------------------------------------

print("Generating speedup plot...")

baseline_sor = times.get("LidDrivenCavity300_sor_serial") or times.get("LidDrivenCavity300_sor_1_1")
baseline_cg  = times.get("LidDrivenCavity300_cg_serial")  or times.get("LidDrivenCavity300_cg_1_1")

if baseline_sor or baseline_cg:
    fig, ax = plt.subplots(figsize=(9, 4))

    if baseline_sor and sor_avail:
        su = [baseline_sor / times[s] for s, *_ in sor_avail]
        lbls = [l for _, l, _ in sor_avail]
        ax.plot(range(len(su)), su, marker="s", color="#2196F3",
                label="SOR Speedup", linewidth=2, linestyle="--")
        for i, v in enumerate(su):
            ax.annotate(f"{v:.2f}×", (i, v), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=10, color="#2196F3")

    if baseline_cg and cg_avail:
        su = [baseline_cg / times[s] for s, *_ in cg_avail]
        lbls = [l for _, l, _ in cg_avail]
        ax.plot(range(len(su)), su, marker="o", color="#4CAF50",
                label="CG Speedup", linewidth=2)
        for i, v in enumerate(su):
            ax.annotate(f"{v:.2f}×", (i, v), textcoords="offset points",
                        xytext=(0, -15), ha="center", fontsize=10, color="#4CAF50")

    tick_lbls = [l.replace("SOR ","").replace("CG ","")
                 for _, l, _ in (sor_avail or cg_avail)]
    ax.set_xticks(range(len(tick_lbls)))
    ax.set_xticklabels(tick_lbls, fontsize=12)
    ax.set_ylabel("Speedup (vs. Serial)", fontsize=13)
    ax.set_title("Strong Scaling Speedup — 150×150 Lid-Driven Cavity", fontsize=14, fontweight="bold")
    ax.axhline(y=1, color="gray", linestyle=":", alpha=0.5)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    save(fig, "plot_speedup.pdf")
    plt.close()

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

print("\n=== Summary Table ===")
print(f"{'Config':<25} {'Wall [s]':>10} {'Speedup':>10} {'Avg iters':>12} {'Max iters':>12}")
print("-" * 72)
for configs, baseline in [(SOR_CONFIGS, baseline_sor), (CG_CONFIGS, baseline_cg)]:
    for stem, label, *_ in configs:
        wt = times.get(stem)
        df = dfs.get(stem)
        su     = f"{baseline/wt:.2f}×" if wt and baseline else "—"
        wt_str = f"{wt:.1f}"           if wt              else "—"
        avg    = f"{df['sor_iters'].mean():.1f}" if df is not None else "—"
        mx     = f"{df['sor_iters'].max()}"      if df is not None else "—"
        print(f"  {label:<23} {wt_str:>10} {su:>10} {avg:>12} {mx:>12}")
    print()

print("Done.")
