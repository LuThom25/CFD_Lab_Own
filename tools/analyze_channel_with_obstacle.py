"""
Analysis: Kármán Vortex Street (Channel With Obstacle) — Worksheet 1.3.2
=========================================================================
This script reads the VTK output files of the ChannelWithObstacle simulations
(Re=2000, Re=200, Re=40) and produces for each case a velocity time series
at probe point (x=3.05, y=1.05) in the wake — showing whether the flow
is unsteady (vortex shedding) or reaches a steady state.

Physical background:
  For a tilted flat plate the onset of periodic vortex shedding occurs at
  Re_D ~ 50–100 (based on effective plate width D ≈ 0.3).
  With D = 0.3:
    Re_D = U * D / nu = 0.3 / nu
  Critical nu ≈ 0.3/50 ≈ 0.006  →  Re_channel = 2/nu ≈ 333

Usage (run from the repository root):
  python3 tools/analyze_channel_with_obstacle.py
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Simulation parameters ─────────────────────────────────────────────────────
IMAX, JMAX = 100, 20
XLENGTH    = 10.0
YLENGTH    = 2.0
DX, DY     = XLENGTH / IMAX, YLENGTH / JMAX
DT_VALUE   = 0.2    # output interval [s]

# Probe point — same as circular obstacle case for direct comparison
# i=30 → x = (30+0.5)*0.1 = 3.05
# j=10 → y = (10+0.5)*0.1 = 1.05
i_probe = 30
j_probe = 10
x_probe = (i_probe + 0.5) * DX
y_probe = (j_probe + 0.5) * DY

# Effective plate width (approximate from PGM geometry)
D_PLATE = 0.3
UIN     = 1.0

# Base directory (relative to this script)
BASE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "example_cases", "ChannelWithObstacle"
)

CASES = [
    ("Re2000", 0.001, "tab:red"),
    ("Re200",  0.01,  "tab:orange"),
    ("Re40",   0.05,  "tab:green"),
]


# ── VTK parser ────────────────────────────────────────────────────────────────
def read_field(lines, keyword, n):
    """
    Finds the line starting with 'keyword' and reads n float values after it.
    Used to extract velocity data from the ASCII VTK FIELD block.
    """
    for idx, line in enumerate(lines):
        if line.strip().startswith(keyword):
            data = []; i = idx + 1
            while len(data) < n:
                data.extend(lines[i].split()); i += 1
            return np.array(data[:n], dtype=float)
    raise ValueError(f"Keyword '{keyword}' not found.")

def parse_vtk(filepath):
    """
    Reads ux and uy cell-centred velocity fields from one VTK file.
    VTK cell ordering: i (x) varies fastest → reshape(JMAX, IMAX).T
    """
    with open(filepath) as f:
        lines = f.readlines()
    n   = IMAX * JMAX
    vel = read_field(lines, "velocity 3 2000", n * 3)
    ux  = vel[0::3].reshape(JMAX, IMAX).T
    uy  = vel[1::3].reshape(JMAX, IMAX).T
    return ux, uy


# ── Read time series for each case ────────────────────────────────────────────
results = {}

for case, nu, color in CASES:
    # Re2000 uses the original (default) output directory name
    if case == "Re2000":
        outdir = os.path.join(BASE, "ChannelWithObstacle_Output")
        prefix = "ChannelWithObstacle"
    else:
        outdir = os.path.join(BASE, f"ChannelWithObstacle_{case}_Output")
        prefix = f"ChannelWithObstacle_{case}"

    files = sorted(
        glob.glob(os.path.join(outdir, f"{prefix}_0.*.vtk")),
        key=lambda f: int(os.path.basename(f).split("_0.")[-1].replace(".vtk", ""))
    )
    if not files:
        print(f"WARNING: No VTK files found for {case} in {outdir}. Skipping.")
        continue

    times = np.arange(len(files)) * DT_VALUE
    ux_p  = np.zeros(len(files))
    uy_p  = np.zeros(len(files))

    print(f"Reading {len(files)} files for {case} ...", flush=True)
    for k, fp in enumerate(files):
        ux, uy   = parse_vtk(fp)
        ux_p[k]  = ux[i_probe, j_probe]
        uy_p[k]  = uy[i_probe, j_probe]

    # Steady-state check over last 50 snapshots (t = 20–30 s)
    last50    = np.sqrt(ux_p[-50:]**2 + uy_p[-50:]**2)
    mean_val  = float(np.mean(last50))
    std_val   = float(np.std(last50))
    variation = std_val / mean_val * 100 if mean_val > 0 else 0.0
    Re_D      = UIN * D_PLATE / nu
    Re_ch     = UIN * YLENGTH / nu

    # Estimate shedding period from uy zero-crossings
    T_shed, St_D = None, None
    if variation > 2.0:
        crossings = np.where(np.diff(np.sign(uy_p[-50:])))[0]
        if len(crossings) >= 2:
            T_shed = 2.0 * float(np.mean(np.diff(times[-50:][crossings])))
            St_D   = D_PLATE / (T_shed * UIN)

    results[case] = dict(
        nu=nu, Re_ch=Re_ch, Re_D=Re_D, color=color, outdir=outdir, prefix=prefix,
        times=times, ux_p=ux_p, uy_p=uy_p,
        mean=mean_val, std=std_val, variation=variation,
        T_shed=T_shed, St_D=St_D
    )

    print(f"  Re_channel = {Re_ch:.0f},  Re_D = {Re_D:.1f}")
    print(f"  Mean |u| = {mean_val:.4f},  Std = {std_val:.5f},  "
          f"Rel. variation = {variation:.2f}%")
    if variation > 2.0:
        shed_str = f"T ≈ {T_shed:.2f} s" if T_shed else "T = n/a"
        st_str   = f"St_D = {St_D:.3f}" if St_D else ""
        print(f"  → UNSTEADY — {shed_str},  {st_str}")
    else:
        print(f"  → STEADY STATE reached")
    print()


# ── Helper: build annotation text ─────────────────────────────────────────────
def make_annotation(r):
    if r["variation"] > 2.0:
        lines = [f"UNSTEADY", f"Rel. variation = {r['variation']:.1f}%"]
        if r["T_shed"] is not None:
            lines.append(f"T ≈ {r['T_shed']:.2f} s")
        if r["St_D"] is not None:
            lines.append(f"St_D = {r['St_D']:.3f}")
        return "\n".join(lines)
    else:
        return f"STEADY STATE\nRel. variation = {r['variation']:.2f}%"


# ── Plot: one subplot per Re case (combined figure) ───────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
fig.suptitle(
    f"Tilted Plate Obstacle — Steady-State Assessment at different Re\n"
    f"Probe point: x = {x_probe:.2f},  y = {y_probe:.2f}   (D_plate ≈ {D_PLATE})",
    fontsize=13, fontweight="bold"
)

for ax, (case, nu, color) in zip(axes, CASES):
    if case not in results:
        ax.set_visible(False)
        continue
    r = results[case]

    # Determine y-axis limits with a 20% margin, minimum ±0.1 headroom
    all_vals  = np.concatenate([r["ux_p"], r["uy_p"]])
    vmin, vmax = float(np.min(all_vals)), float(np.max(all_vals))
    margin    = max(0.15 * (vmax - vmin), 0.1)
    ax.set_ylim(vmin - margin, vmax + margin)

    ax.plot(r["times"], r["ux_p"], color=color, linewidth=1.2, label="ux")
    ax.plot(r["times"], r["uy_p"], color=color, linewidth=1.2,
            linestyle="--", label="uy")
    ax.axhline(0, color="k", linewidth=0.5, alpha=0.4)
    ax.set_ylabel("velocity  [–]", fontsize=10)
    ax.set_title(
        f"Re_channel = {r['Re_ch']:.0f}  |  Re_D = {r['Re_D']:.0f}  "
        f"(ν = {nu})", fontsize=11
    )
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)

    # Annotation box — lower right
    ax.text(0.98, 0.04, make_annotation(r),
            transform=ax.transAxes, fontsize=9, va="bottom", ha="right",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))

axes[-1].set_xlabel("t  [s]", fontsize=11)
plt.tight_layout()

# Save combined comparison plot in the ChannelWithObstacle root folder
combined = os.path.join(BASE, "ChannelWithObstacle_ReComparison.png")
fig.savefig(combined, dpi=150, bbox_inches="tight")
print(f"Combined plot saved: {combined}")

# Save individual plot in each output directory
for case, nu, color in CASES:
    if case not in results:
        continue
    r = results[case]

    fig_i, ax_i = plt.subplots(figsize=(13, 4.5))
    fig_i.suptitle(
        f"Tilted Plate Obstacle  —  {case}  "
        f"[Re_channel = {r['Re_ch']:.0f},  Re_D = {r['Re_D']:.0f},  grid {IMAX}×{JMAX}]",
        fontsize=13, fontweight="bold"
    )

    all_vals  = np.concatenate([r["ux_p"], r["uy_p"]])
    vmin, vmax = float(np.min(all_vals)), float(np.max(all_vals))
    margin    = max(0.15 * (vmax - vmin), 0.1)
    ax_i.set_ylim(vmin - margin, vmax + margin)

    ax_i.plot(r["times"], r["ux_p"], color=r["color"], linewidth=1.2,
              label=f"ux  at  ({x_probe:.2f}, {y_probe:.2f})")
    ax_i.plot(r["times"], r["uy_p"], color=r["color"], linewidth=1.2,
              linestyle="--", label=f"uy  at  ({x_probe:.2f}, {y_probe:.2f})")
    ax_i.axhline(0, color="k", linewidth=0.5, alpha=0.5)
    ax_i.set_xlabel("t  [s]", fontsize=11)
    ax_i.set_ylabel("velocity  [–]", fontsize=11)
    ax_i.set_title(
        f"Velocity time series at probe point in the wake\n"
        f"Re_channel = {r['Re_ch']:.0f}  |  Re_D = {r['Re_D']:.0f}  (D ≈ {D_PLATE})",
        fontsize=11
    )
    ax_i.legend(fontsize=10, loc="upper left")
    ax_i.grid(True, alpha=0.35)

    ax_i.text(0.98, 0.04, make_annotation(r),
              transform=ax_i.transAxes, fontsize=9, va="bottom", ha="right",
              bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))

    plt.tight_layout()
    out_path = os.path.join(r["outdir"], f"{r['prefix']}_analysis.png")
    fig_i.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig_i)
    print(f"Individual plot saved: {out_path}")
