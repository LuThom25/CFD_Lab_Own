"""
WS3 Convergence Study — 300x300 Lid-Driven Cavity
Produces 3 SOR-iteration plots (grouped by process count),
a runtime bar chart, and a speedup plot.

Run after all simulations are complete:
    python3 analyze.py
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE = os.path.dirname(os.path.abspath(__file__))
DOWNLOADS = os.path.expanduser("~/Downloads")

# Each entry: (dat_stem, display_label, n_procs, line_style, color)
# dat_stem determines the output directory via:
#   {dat_stem}_Output_{iproc}_{jproc}/
ALL_CONFIGS = [
    # stem                        label        np  ls      color
    ("LidDrivenCavity300_serial", "Serial",    1,  "-",    "#2196F3"),
    ("LidDrivenCavity300_1_1",    "MPI (1×1)", 1,  "--",   "#FF9800"),
    ("LidDrivenCavity300_2_2",    "MPI (2×2)", 4,  "-",    "#4CAF50"),
    ("LidDrivenCavity300_1_4",    "MPI (1×4)", 4,  "--",   "#F44336"),
    ("LidDrivenCavity300_4_1",    "MPI (4×1)", 4,  ":",    "#9C27B0"),
    ("LidDrivenCavity300_3_2",    "MPI (3×2)", 6,  "-",    "#009688"),
    ("LidDrivenCavity300_2_3",    "MPI (2×3)", 6,  "--",   "#FF5722"),
    ("LidDrivenCavity300_3_3",    "MPI (3×3)", 9,  "-.",   "#795548"),
]

# The 3 groups for the SOR-iteration plots
GROUPS = [
    {
        "title":   "SOR iterations — Serial vs MPI (1×1)",
        "stems":   ["LidDrivenCavity300_serial", "LidDrivenCavity300_1_1"],
        "outfile": "plot_sor_serial_vs_11.pdf",
    },
    {
        "title":   "SOR iterations — 4-rank decompositions",
        "stems":   ["LidDrivenCavity300_4_1",
                    "LidDrivenCavity300_2_2",
                    "LidDrivenCavity300_1_4"],
        "outfile": "plot_sor_4ranks.pdf",
    },
    {
        "title":   "SOR iterations — 6- and 9-rank decompositions",
        "stems":   ["LidDrivenCavity300_3_2",
                    "LidDrivenCavity300_2_3",
                    "LidDrivenCavity300_3_3"],
        "outfile": "plot_sor_6_9ranks.pdf",
    },
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def output_dir(stem, iproc, jproc):
    return os.path.join(BASE, f"{stem}_Output_{iproc}_{jproc}")

def iproc_jproc_from_stem(stem):
    """Extract iproc/jproc from the dat stem or from the config table."""
    for s, _, np, ls, col in ALL_CONFIGS:
        if s == stem:
            # parse from stem suffix e.g. "..._1_4"
            parts = stem.rsplit("_", 2)
            if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
                return int(parts[1]), int(parts[2])
            return 1, 1  # serial fallback
    return 1, 1

def load_csv(stem):
    ip, jp = iproc_jproc_from_stem(stem)
    d = output_dir(stem, ip, jp)
    csv = os.path.join(d, "sor_log.csv")
    if not os.path.exists(csv):
        return None
    return pd.read_csv(csv)

def load_walltime(stem):
    ip, jp = iproc_jproc_from_stem(stem)
    d = output_dir(stem, ip, jp)
    wt = os.path.join(d, "walltime.txt")
    if not os.path.exists(wt):
        return None
    with open(wt) as f:
        return float(f.read().strip())

def cfg(stem):
    """Return (label, n_procs, linestyle, color) for a stem."""
    for s, lbl, np, ls, col in ALL_CONFIGS:
        if s == stem:
            return lbl, np, ls, col
    return stem, 1, "-", "black"

# ---------------------------------------------------------------------------
# Load all data
# ---------------------------------------------------------------------------

print("Loading data...")
dfs   = {}   # stem -> DataFrame or None
times = {}   # stem -> float or None

for stem, label, *_ in ALL_CONFIGS:
    df = load_csv(stem)
    wt = load_walltime(stem)
    dfs[stem]   = df
    times[stem] = wt
    status = f"{len(df)} steps, {wt:.1f}s" if df is not None and wt else \
             f"{len(df)} steps" if df is not None else "not found"
    print(f"  {label:15s}: {status}")

print()

# ---------------------------------------------------------------------------
# Plot helper
# ---------------------------------------------------------------------------

def sor_plot(group):
    stems = group["stems"]
    available = [s for s in stems if dfs[s] is not None]
    if not available:
        print(f"  [skip] {group['title']} — no data")
        return

    fig, ax = plt.subplots(figsize=(9, 4))

    for stem in stems:
        if dfs[stem] is None:
            continue
        df = dfs[stem]
        lbl, np_, ls, col = cfg(stem)
        ax.plot(df["timestep"], df["sor_iters"],
                label=lbl, linestyle=ls, color=col,
                linewidth=1.4, alpha=0.9)

    ax.set_xlabel("Timestep (step)", fontsize=13)
    ax.set_ylabel("SOR Iterations (iters)", fontsize=13)
    ax.set_title(group["title"], fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=11, framealpha=0.8)
    ax.tick_params(axis="both", labelsize=11)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    out = os.path.join(BASE, group["outfile"])
    plt.savefig(out, dpi=150)
    out_dl = os.path.join(DOWNLOADS, group["outfile"])
    plt.savefig(out_dl, dpi=150)
    print(f"  Saved: {group['outfile']}")
    plt.close()

# ---------------------------------------------------------------------------
# Produce SOR plots
# ---------------------------------------------------------------------------

print("Generating SOR plots...")
for group in GROUPS:
    sor_plot(group)

# ---------------------------------------------------------------------------
# Combined runtime plot: walltime bar (left) + speedup curve (right)
# ---------------------------------------------------------------------------

print("Generating runtime plot...")

# Collect available data
avail = []
for stem, label, np_, ls, col in ALL_CONFIGS:
    if times[stem] is not None:
        avail.append((stem, label, np_, col, times[stem]))

# Baseline for speedup (prefer serial, fallback to 1x1)
baseline_stem = "LidDrivenCavity300_serial"
if times.get(baseline_stem) is None:
    baseline_stem = "LidDrivenCavity300_1_1"

if avail:
    lbls  = [a[1] for a in avail]
    vals  = [a[4] for a in avail]
    cols  = [a[3] for a in avail]

    fig, ax_bar = plt.subplots(figsize=(9, 4))

    # --- Wall time bar chart ---
    bars = ax_bar.bar(lbls, vals, color=cols, edgecolor="white", width=0.6)
    for bar, val in zip(bars, vals):
        ax_bar.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(vals) * 0.01,
                    f"{val:.0f}s", ha="center", va="bottom", fontsize=11)
    ax_bar.set_ylabel("Wall time (s)", fontsize=13)
    ax_bar.set_title("Wall time per setup", fontsize=14, fontweight="bold")
    ax_bar.grid(True, axis="y", linestyle="--", alpha=0.35)
    ax_bar.tick_params(axis="x", rotation=25, labelsize=11)
    ax_bar.tick_params(axis="y", labelsize=11)
    ax_bar.set_ylim(top=max(vals) * 1.12)

    plt.tight_layout()
    out = os.path.join(BASE, "plot_runtime.pdf")
    plt.savefig(out, dpi=150)
    out_dl = os.path.join(DOWNLOADS, "plot_runtime.pdf")
    plt.savefig(out_dl, dpi=150)
    print(f"  Saved: plot_runtime.pdf")
    plt.close()

print("\nDone.")
