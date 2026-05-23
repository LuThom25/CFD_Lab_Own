"""
Analysis: Channel With Circular Obstacle — Re comparison
=========================================================
This script reads the VTK output files of the CircleObstacle simulations
(Re=2000, Re=200, Re=40) and produces for each case a velocity time series
at probe point (x=3.05, y=1.05) in the wake — showing whether the flow
is unsteady (vortex shedding) or reaches a steady state.

All three cases use the same probe point as the tilted-plate
ChannelWithObstacle case, allowing direct comparison.

Physical background:
  For a circular cylinder the onset of periodic vortex shedding occurs at
  Re_D ≈ 47 (based on cylinder diameter D). Below this threshold the wake
  is steady; above it alternating vortices are shed periodically.
  With D ≈ 0.7 in our setup:
    Re_D = U * D / nu = 0.7 / nu
  Critical nu ≈ 0.7/47 ≈ 0.015  →  Re_channel = 2/nu ≈ 133

Usage (run from the repository root):
  python3 tools/analyze_circle_obstacle.py
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

# Probe point — same as tilted-plate case for direct comparison
# i=30 → x = (30+0.5)*0.1 = 3.05
# j=10 → y = (10+0.5)*0.1 = 1.05
i_probe = 30
j_probe = 10
x_probe = (i_probe + 0.5) * DX
y_probe = (j_probe + 0.5) * DY

# Effective cylinder diameter (radius=3.5 cells * dy=0.1 * 2)
D = 0.7

# Base directory (relative to this script)
BASE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "example_cases", "ChannelWithCircleObstacle"
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
    outdir = os.path.join(BASE, f"CircleObstacle_{case}_Output")
    files  = sorted(
        glob.glob(os.path.join(outdir, f"CircleObstacle_{case}_0.*.vtk")),
        key=lambda f: int(os.path.basename(f).split("_0.")[-1].replace(".vtk", ""))
    )
    times = np.arange(len(files)) * DT_VALUE
    ux_p  = np.zeros(len(files))
    uy_p  = np.zeros(len(files))

    print(f"Reading {len(files)} files for {case} ...", flush=True)
    for k, fp in enumerate(files):
        ux, uy   = parse_vtk(fp)
        ux_p[k]  = ux[i_probe, j_probe]
        uy_p[k]  = uy[i_probe, j_probe]

    # Steady-state check over last 50 snapshots (t = 20-30 s)
    last50    = np.sqrt(ux_p[-50:]**2 + uy_p[-50:]**2)
    mean_val  = float(np.mean(last50))
    std_val   = float(np.std(last50))
    variation = std_val / mean_val * 100 if mean_val > 0 else 0.0
    Re_D      = D / nu
    Re_ch     = YLENGTH / nu

    # Estimate shedding period from uy zero-crossings
    T_shed, St_D = None, None
    if variation > 2.0:
        crossings = np.where(np.diff(np.sign(uy_p[-50:])))[0]
        if len(crossings) >= 2:
            T_shed = 2.0 * float(np.mean(np.diff(times[-50:][crossings])))
            St_D   = D / (T_shed * 1.0)

    results[case] = dict(
        nu=nu, Re_ch=Re_ch, Re_D=Re_D, color=color,
        times=times, ux_p=ux_p, uy_p=uy_p,
        mean=mean_val, std=std_val, variation=variation,
        T_shed=T_shed, St_D=St_D
    )

    print(f"  Re_channel = {Re_ch:.0f},  Re_D = {Re_D:.1f}")
    print(f"  Mean |u| = {mean_val:.4f},  Std = {std_val:.5f},  "
          f"Rel. variation = {variation:.2f}%")
    if variation > 2.0:
        print(f"  → UNSTEADY — T ≈ {T_shed:.2f} s,  St_D = {St_D:.3f}")
    else:
        print(f"  → STEADY STATE reached")
    print()


# ── Plot: one subplot per Re case ─────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
fig.suptitle(
    f"Circular Obstacle — Steady-State Assessment at different Re\n"
    f"Probe point: x = {x_probe:.2f},  y = {y_probe:.2f}   (D ≈ {D})",
    fontsize=13, fontweight="bold"
)

for ax, (case, nu, color) in zip(axes, CASES):
    r = results[case]
    ax.plot(r["times"], r["ux_p"], color=color, linewidth=1.2,
            label="ux")
    ax.plot(r["times"], r["uy_p"], color=color, linewidth=1.2,
            linestyle="--", label="uy")
    ax.axhline(0, color="k", linewidth=0.5, alpha=0.4)
    ax.set_ylabel("velocity  [–]", fontsize=10)
    ax.set_title(
        f"Re_channel = {r['Re_ch']:.0f}  |  Re_D = {r['Re_D']:.0f}  "
        f"(ν = {nu})", fontsize=11
    )
    # Legend in upper left — away from the annotation box
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)

    # Status annotation box in lower right
    if r["variation"] > 2.0:
        status = (f"UNSTEADY\nRel. variation = {r['variation']:.1f}%\n"
                  f"T ≈ {r['T_shed']:.2f} s,  St_D = {r['St_D']:.3f}")
    else:
        status = f"STEADY STATE\nRel. variation = {r['variation']:.2f}%"

    ax.text(0.98, 0.04, status,
            transform=ax.transAxes, fontsize=9, va="bottom", ha="right",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))

axes[-1].set_xlabel("t  [s]", fontsize=11)
plt.tight_layout()

# Save individual plot in each output directory
for case, nu, color in CASES:
    outdir   = os.path.join(BASE, f"CircleObstacle_{case}_Output")
    out_path = os.path.join(outdir, f"CircleObstacle_{case}_analysis.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {out_path}")

# Save combined comparison plot in the CircleObstacle root folder
combined = os.path.join(BASE, "CircleObstacle_ReComparison.png")
fig.savefig(combined, dpi=150, bbox_inches="tight")
print(f"Combined plot saved: {combined}")
