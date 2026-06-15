"""
Strong Scaling Analysis — 40×18 Rayleigh-Bénard Convection
Decompositions: (1×1), (2×2), (4×1)  — max 2 subdomains in vertical direction.
Compares SOR_STANDARD and CG_STANDARD:
  - Solver iterations per timestep
  - Wall time and speedup

Run after all simulations are complete:
    python3 analyze_rb.py
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

# (dat_stem, display_label, n_procs, decomp_label, line_style, color)
SOR_CONFIGS = [
    ("RayleighBenard_sor_1_1", "SOR (1×1)", 1,  "1×1", "-",  "#2196F3"),
    ("RayleighBenard_sor_2_2", "SOR (2×2)", 4,  "2×2", "--", "#FF9800"),
    ("RayleighBenard_sor_4_1", "SOR (4×1)", 4,  "4×1", ":",  "#9C27B0"),
]

CG_CONFIGS = [
    ("RayleighBenard_cg_1_1", "CG (1×1)", 1, "1×1", "-",  "#4CAF50"),
    ("RayleighBenard_cg_2_2", "CG (2×2)", 4, "2×2", "--", "#F44336"),
    ("RayleighBenard_cg_4_1", "CG (4×1)", 4, "4×1", ":",  "#009688"),
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
    if not os.path.exists(csv):
        return None
    return pd.read_csv(csv)

def load_walltime(stem):
    wt = os.path.join(output_dir(stem), "walltime.txt")
    if not os.path.exists(wt):
        return None
    with open(wt) as f:
        return float(f.read().strip())

# ---------------------------------------------------------------------------
# Load all data
# ---------------------------------------------------------------------------

print("Loading simulation data...")
dfs   = {}
times = {}

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
# Plot 1: Solver iterations per timestep — SOR vs CG for each decomposition
# ---------------------------------------------------------------------------

print("Generating solver-iteration plots...")

def iters_plot(configs, title, outfile):
    available = [(s, l, np_, _, ls, col) for s, l, np_, _, ls, col in configs
                 if dfs.get(s) is not None]
    if not available:
        print(f"  [skip] {title} — no data")
        return
    fig, ax = plt.subplots(figsize=(9, 4))
    for stem, label, _, _, ls, col in configs:
        if dfs.get(stem) is None:
            continue
        df = dfs[stem]
        ax.plot(df["timestep"], df["sor_iters"],
                label=label, linestyle=ls, color=col, linewidth=1.4, alpha=0.9)
    ax.set_xlabel("Timestep", fontsize=13)
    ax.set_ylabel("Solver Iterations", fontsize=13)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=10, framealpha=0.8)
    ax.tick_params(axis="both", labelsize=11)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    for dest in [BASE, DOWNLOADS]:
        plt.savefig(os.path.join(dest, outfile), dpi=150)
    print(f"  Saved: {outfile}")
    plt.close()

iters_plot(SOR_CONFIGS,
           "SOR iterations/step — Rayleigh-Bénard (40×18)",
           "plot_rb_sor_iters.pdf")
iters_plot(CG_CONFIGS,
           "CG iterations/step — Rayleigh-Bénard (40×18)",
           "plot_rb_cg_iters.pdf")

# Combined: SOR vs CG for 1×1 decomposition
combined_1_1 = [
    ("RayleighBenard_sor_1_1", "SOR (1×1)", 1, "1×1", "-",  "#2196F3"),
    ("RayleighBenard_cg_1_1",  "CG (1×1)",  1, "1×1", "-",  "#4CAF50"),
]
iters_plot(combined_1_1,
           "SOR vs CG iterations/step — Rayleigh-Bénard 1×1",
           "plot_rb_sor_vs_cg_1_1.pdf")

# ---------------------------------------------------------------------------
# Plot 2: Wall time comparison SOR vs CG
# ---------------------------------------------------------------------------

print("Generating wall time comparison...")

fig, ax = plt.subplots(figsize=(10, 5))
w = 0.35

sor_lbls = [l for s, l, *_ in SOR_CONFIGS if times.get(s) is not None]
sor_vals = [times[s] for s, *_ in SOR_CONFIGS if times.get(s) is not None]
sor_cols = [col for s, l, np_, _, ls, col in SOR_CONFIGS if times.get(s) is not None]

cg_lbls = [l for s, l, *_ in CG_CONFIGS if times.get(s) is not None]
cg_vals = [times[s] for s, *_ in CG_CONFIGS if times.get(s) is not None]
cg_cols = [col for s, l, np_, _, ls, col in CG_CONFIGS if times.get(s) is not None]

n = max(len(sor_vals), len(cg_vals))

if sor_vals:
    bars = ax.bar([i - w/2 for i in range(len(sor_vals))], sor_vals,
                  width=w, color="#2196F3", label="SOR_STANDARD", edgecolor="white", alpha=0.85)
    for bar, val in zip(bars, sor_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(sor_vals + cg_vals) * 0.01,
                f"{val:.1f}s", ha="center", va="bottom", fontsize=10)

if cg_vals:
    bars = ax.bar([i + w/2 for i in range(len(cg_vals))], cg_vals,
                  width=w, color="#4CAF50", label="CG_STANDARD", edgecolor="white", alpha=0.85)
    for bar, val in zip(bars, cg_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(sor_vals + cg_vals) * 0.01,
                f"{val:.1f}s", ha="center", va="bottom", fontsize=10)

tick_lbls = [l.replace("SOR ", "").replace("CG ", "") for l in sor_lbls or cg_lbls]
ax.set_xticks(range(n))
ax.set_xticklabels(tick_lbls, fontsize=12)
ax.set_ylabel("Wall time (s)", fontsize=13)
ax.set_title("Wall time: SOR vs CG — Rayleigh-Bénard (40×18)", fontsize=14, fontweight="bold")
ax.legend(fontsize=11)
ax.grid(True, axis="y", linestyle="--", alpha=0.35)
all_vals = sor_vals + cg_vals
if all_vals:
    ax.set_ylim(top=max(all_vals) * 1.15)

plt.tight_layout()
for dest in [BASE, DOWNLOADS]:
    plt.savefig(os.path.join(dest, "plot_rb_walltime.pdf"), dpi=150)
print("  Saved: plot_rb_walltime.pdf")
plt.close()

# ---------------------------------------------------------------------------
# Plot 3: Speedup — SOR and CG, baseline = serial (1×1)
# ---------------------------------------------------------------------------

print("Generating speedup plot...")

baseline_sor = times.get("RayleighBenard_sor_1_1")
baseline_cg  = times.get("RayleighBenard_cg_1_1")

fig, ax = plt.subplots(figsize=(9, 4))
decomp_labels = ["1×1", "2×2", "4×1"]

if baseline_sor and sor_vals:
    su_sor = [baseline_sor / t for t in sor_vals]
    ax.plot(range(len(su_sor)), su_sor, marker="s", color="#2196F3",
            label="SOR Speedup", linewidth=2, linestyle="--")
    for i, su in enumerate(su_sor):
        ax.annotate(f"{su:.2f}×", (i, su), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=10, color="#2196F3")

if baseline_cg and cg_vals:
    su_cg = [baseline_cg / t for t in cg_vals]
    ax.plot(range(len(su_cg)), su_cg, marker="o", color="#4CAF50",
            label="CG Speedup", linewidth=2)
    for i, su in enumerate(su_cg):
        ax.annotate(f"{su:.2f}×", (i, su), textcoords="offset points",
                    xytext=(0, -15), ha="center", fontsize=10, color="#4CAF50")

ax.set_xticks(range(len(decomp_labels[:n])))
ax.set_xticklabels(decomp_labels[:n], fontsize=12)
ax.set_ylabel("Speedup (vs. 1×1 baseline)", fontsize=13)
ax.set_title("Strong Scaling Speedup — Rayleigh-Bénard (40×18)", fontsize=14, fontweight="bold")
ax.axhline(y=1, color="gray", linestyle=":", alpha=0.5)
ax.legend(fontsize=11)
ax.grid(True, linestyle="--", alpha=0.35)
ax.set_ylim(bottom=0)

plt.tight_layout()
for dest in [BASE, DOWNLOADS]:
    plt.savefig(os.path.join(dest, "plot_rb_speedup.pdf"), dpi=150)
print("  Saved: plot_rb_speedup.pdf")
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
        su = f"{baseline/wt:.2f}×" if wt and baseline else "—"
        wt_str = f"{wt:.1f}" if wt else "—"
        avg = f"{df['sor_iters'].mean():.1f}" if df is not None else "—"
        mx  = f"{df['sor_iters'].max()}" if df is not None else "—"
        print(f"  {label:<23} {wt_str:>10} {su:>10} {avg:>12} {mx:>12}")
    print()

print("Done.")
