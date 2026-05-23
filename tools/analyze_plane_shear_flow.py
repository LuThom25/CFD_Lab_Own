"""
Analysis: Plane Shear Flow (Poiseuille Flow) — Worksheet 1.3.1
==============================================================
This script reads the last VTK output file of the PlaneShearFlow simulation
and produces two plots:

  1) Velocity profile u(y) at the centre cross-section x = xlength/2
     → simulated values compared against the analytical Poiseuille solution

  2) Pressure profile along the channel centre line
     → verifies the linear pressure drop and computes dp/dx

Analytical solution (derived from the incompressible NS momentum equation):
  The x-momentum equation for fully-developed, steady, 2D channel flow reduces to:
      0 = -dp/dx  +  nu * d²u/dy²
  =>  d²u/dy² = (1/nu) * dp/dx = Re * dp/dx
  Integrating twice with no-slip BCs u(0)=0, u(h)=0:
      u(y) = (Re/2) * (dp/dx) * y * (y - h)
  The pressure gradient dp/dx is obtained from the mean-velocity condition:
      u_mean = UIN = 1.0  =>  dp/dx = -12*nu*u_mean / h^2 = -0.3
  Substituting gives the parabolic profile:
      u(y) = (10/2) * (-0.3) * y * (y-2) = 1.5 * y * (2-y)

Usage (run from the repository root):
  python3 tools/analyze_plane_shear_flow.py
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend: saves directly to file, no display needed
import matplotlib.pyplot as plt


# ── Simulation parameters (must match the .dat file) ─────────────────────────
# These define the grid size and physical dimensions of the domain.
IMAX    = 100     # number of fluid cells in x-direction
JMAX    = 20      # number of fluid cells in y-direction
XLENGTH = 10.0    # domain length in x [m or dimensionless units]
YLENGTH = 2.0     # domain height in y (= channel width h)
NU      = 0.1     # kinematic viscosity (nu = 1/Re in dimensionless form)
UIN     = 1.0     # prescribed inlet velocity = mean channel velocity

# Derived quantities
DX = XLENGTH / IMAX    # cell size in x = 0.10
DY = YLENGTH / JMAX    # cell size in y = 0.10
RE = 1.0 / NU          # Reynolds number = 10
H  = YLENGTH           # channel height = 2.0

# Analytical pressure gradient derived from u_mean = UIN:
#   Integrating u(y) over [0,h] and dividing by h gives:
#   u_mean = -(Re * dp/dx * h^2) / 12
#   => dp/dx = -12 * nu * u_mean / h^2 = -12 * 0.1 * 1.0 / 4 = -0.3
# Negative sign: pressure decreases in the flow direction (+x), which drives the flow.
DPDX = -12.0 * NU * UIN / H**2   # = -0.3

# Path to the VTK output directory (relative to this script's location)
OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "example_cases", "PlaneShearFlow", "PlaneShearFlow_Output"
)


# ── VTK file parser ───────────────────────────────────────────────────────────
def read_field(lines, keyword, n_total):
    """
    Searches the list of file lines for the first line that starts with
    'keyword', then reads exactly n_total floating-point numbers from the
    lines that follow it.

    The ASCII VTK FIELD format stores each named array as:
        <name> <n_components> <n_tuples> <datatype>
        value0 value1 value2 ...

    Parameters
    ----------
    lines   : list of str — all lines of the VTK file
    keyword : str         — beginning of the header line to search for
                            (e.g. "pressure 1 2000" or "velocity 3 2000")
    n_total : int         — total number of scalar values to read
                            (= n_components * n_tuples)

    Returns
    -------
    numpy array of shape (n_total,) with dtype float64
    """
    for idx, line in enumerate(lines):
        if line.strip().startswith(keyword):
            data = []
            i = idx + 1
            # Keep reading lines until we have collected enough numbers.
            # Lines may contain varying numbers of values (VTK writes ~9 per line).
            while len(data) < n_total:
                data.extend(lines[i].split())
                i += 1
            return np.array(data[:n_total], dtype=float)
    raise ValueError(f"VTK keyword '{keyword}' not found in file.")


def parse_vtk(filepath):
    """
    Parses one ASCII VTK StructuredGrid file and extracts the cell-centred
    pressure and x-velocity fields.

    VTK cell ordering for a structured grid with IMAX × JMAX cells:
      - Index varies fastest in x (i-direction), slowest in y (j-direction).
      - Cell (i, j)  <->  flat index  k = j * IMAX + i
      - Flat array reshaped as (JMAX, IMAX) gives [j, i] indexing.
      - Transposing gives (IMAX, JMAX) with [i, j] indexing (more natural).

    Parameters
    ----------
    filepath : str — full path to the .vtk file

    Returns
    -------
    p_2d  : ndarray (IMAX, JMAX) — pressure at each cell centre
    ux_2d : ndarray (IMAX, JMAX) — x-velocity at each cell centre
    """
    with open(filepath) as f:
        lines = f.readlines()

    n_cells = IMAX * JMAX    # 100 * 20 = 2000 cells total

    # Read pressure: 1 scalar per cell → n_total = 2000
    p_flat = read_field(lines, "pressure 1 2000", n_cells)

    # Read cell-centred velocity: 3 components (ux, uy, uz) per cell → n_total = 6000
    # We use "velocity 3 2000" to match the CELL_DATA block (not "velocity 3 2121"
    # which belongs to the POINT_DATA block written separately by the solver).
    vel_flat = read_field(lines, "velocity 3 2000", n_cells * 3)

    # Extract only the x-component (ux): every 3rd value starting at index 0
    ux_flat = vel_flat[0::3]

    # Reshape from flat (2000,) to 2D (IMAX, JMAX):
    #   Step 1: reshape(JMAX, IMAX)  →  array[j, i]  (j is the outer/slow index)
    #   Step 2: .T                    →  array[i, j]  (more convenient for column access)
    p_2d  = p_flat.reshape(JMAX, IMAX).T
    ux_2d = ux_flat.reshape(JMAX, IMAX).T

    return p_2d, ux_2d


# ── Load the last (= steady-state) VTK file ───────────────────────────────────
# The VTK files are named PlaneShearFlow_0.<timestep>.vtk, where <timestep>
# is the integer time-step counter. The last file corresponds to t = t_end = 30 s.
# We sort by the integer timestep (not lexicographically) to get the correct order.
vtk_files = sorted(
    glob.glob(os.path.join(OUTPUT_DIR, "PlaneShearFlow_0.*.vtk")),
    key=lambda f: int(f.split("_0.")[-1].replace(".vtk", ""))
)
if not vtk_files:
    raise FileNotFoundError(
        f"No VTK files found in {OUTPUT_DIR}.\n"
        "Run the simulation first:\n"
        "  ./build/fluidchen example_cases/PlaneShearFlow/PlaneShearFlow.dat"
    )

latest = vtk_files[-1]
step   = int(latest.split("_0.")[-1].replace(".vtk", ""))
# Physical time: each time step uses dt_adaptive ≈ 0.0125 s
# (limited by the viscous CFL: tau * Re * dx^2*dy^2 / (2*(dx^2+dy^2)) = 0.5*10*0.01*0.01/0.02 = 0.0125)
t_phys = step * 0.0125
print(f"Reading: {os.path.basename(latest)}  (step {step}, t ≈ {t_phys:.1f} s)")

p_2d, ux_2d = parse_vtk(latest)


# ── Cell-centre coordinates ───────────────────────────────────────────────────
# Cell (i, j) has its centre at x = (i + 0.5)*dx, y = (j + 0.5)*dy.
# These are the physical positions at which the VTK cell data is stored.
x_centers = (np.arange(IMAX) + 0.5) * DX    # shape (100,): 0.05, 0.15, ..., 9.95
y_centers = (np.arange(JMAX) + 0.5) * DY    # shape ( 20,): 0.05, 0.15, ..., 1.95


# ── 1) Velocity profile u(y) at the centre cross-section ─────────────────────
# We extract the column at i = IMAX//2 = 50, which corresponds to x ≈ 5.05.
# At this position the flow is already fully developed (entry length ≈ 0.06*Re*h ≈ 1.2).
i_mid = IMAX // 2                    # column index i = 50
x_mid = x_centers[i_mid]            # physical x ≈ 5.05
u_sim = ux_2d[i_mid, :]             # u at all 20 cell centres for this column

# Analytical solution on a fine y-grid for a smooth curve,
# and at the same y-positions as the simulation cells for error computation.
y_fine       = np.linspace(0.0, H, 500)
u_ana_fine   = (RE / 2.0) * DPDX * y_fine   * (y_fine   - H)  # smooth reference curve
u_ana_coarse = (RE / 2.0) * DPDX * y_centers * (y_centers - H) # at cell centres

# Point-wise absolute error between simulation and analytical solution
abs_error = np.abs(u_sim - u_ana_coarse)
max_err   = float(np.max(abs_error))
rel_err   = max_err / float(np.max(u_ana_fine)) * 100.0   # relative to u_max = 1.5

print(f"\nVelocity profile u(y) at x ≈ {x_mid:.2f}:")
print(f"  Max absolute error : {max_err:.5f}  m/s")
print(f"  Relative error     : {rel_err:.2f}%  (w.r.t. u_max = {np.max(u_ana_fine):.2f} m/s)")


# ── 2) Pressure drop along the channel centre line ────────────────────────────
# We take the row at j = JMAX//2 = 10, corresponding to y ≈ 1.05 (mid-height).
# For fully-developed Poiseuille flow the pressure is uniform across the height
# at any given x, so the exact j-choice does not matter much.
j_mid = JMAX // 2
y_mid = y_centers[j_mid]            # physical y ≈ 1.05
p_x   = p_2d[:, j_mid]             # pressure at all x-positions along this line

# Total pressure drop Δp = p_inlet − p_outlet.
# The outflow BC sets p = 0 at the right boundary, so the outlet cell pressure
# is close to (but not exactly) zero.
delta_p_sim = float(p_x[0] - p_x[-1])

# Analytical pressure drop: Δp = |dp/dx| * L
delta_p_ana = abs(DPDX) * XLENGTH   # 0.3 * 10 = 3.0

# Linear regression of p(x) to extract the slope dp/dx from the simulation.
# For fully-developed flow dp/dx should be constant along x.
coeffs   = np.polyfit(x_centers, p_x, 1)
dpdx_fit = float(coeffs[0])

print(f"\nPressure drop (Δp = p_inlet − p_outlet):")
print(f"  Simulation : Δp = {delta_p_sim:.4f}   dp/dx = {dpdx_fit:.4f}")
print(f"  Analytical : Δp = {delta_p_ana:.4f}   dp/dx = {DPDX:.4f}")
print(f"  Deviation  : {abs(delta_p_sim - delta_p_ana) / delta_p_ana * 100:.2f}%")
print(f"  Note: ~4% deviation is expected because the inlet block profile")
print(f"        causes extra pressure loss in the entry region (x < ~1.2).")


# ── Plots ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle(
    f"Plane Shear Flow  —  t = {t_phys:.1f} s   "
    f"[Re = {RE:.0f},  grid {IMAX}×{JMAX}]",
    fontsize=13, fontweight="bold"
)

# ---- Left plot: velocity profile ----
ax = axes[0]

# Analytical solution: smooth black line as the ground-truth reference
ax.plot(u_ana_fine, y_fine, "k-", linewidth=2.5,
        label=r"Analytical:  $u(y) = 1.5\,y\,(2-y)$")

# Simulated profile: open red circles at cell-centre y-positions
ax.plot(u_sim, y_centers, "ro", markersize=7, markerfacecolor="none", linewidth=0,
        label=f"Simulation  (x = {x_mid:.2f})")

ax.set_xlabel("x-velocity  u  [–]", fontsize=12)
ax.set_ylabel("y  [–]", fontsize=12)
ax.set_title("Velocity profile u(y)\nat cross-section x = xlength / 2", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.35)
ax.set_xlim(left=0.0)
ax.set_ylim(0.0, H)

# Error annotation box — placed in the lower-right corner to avoid overlapping the legend
ax.text(0.97, 0.08,
        f"Max. error:  {max_err:.4f}\nRel. error:  {rel_err:.2f}%",
        transform=ax.transAxes, fontsize=9, va="bottom", ha="right",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))

# ---- Right plot: pressure profile ----
ax = axes[1]

# Simulated pressure along the centre line
ax.plot(x_centers, p_x, "b.-", markersize=5, linewidth=1.5,
        label=f"Simulation  (y = {y_mid:.2f})")

# Linear fit through the simulated pressure (should be a straight line for developed flow)
x_ends = np.array([x_centers[0], x_centers[-1]])
ax.plot(x_ends, np.polyval(coeffs, x_ends), "k--", linewidth=2.0,
        label=f"Linear fit:  dp/dx = {dpdx_fit:.4f}")

# Analytical pressure line: p(x) = |dp/dx| * (L - x)  with p(L) = 0 (outflow BC)
ax.plot(x_ends, abs(DPDX) * (XLENGTH - x_ends), "g-", linewidth=2.0, alpha=0.75,
        label=f"Analytical:  dp/dx = {DPDX:.4f}")

ax.set_xlabel("x  [–]", fontsize=12)
ax.set_ylabel("Pressure  p  [–]", fontsize=12)
ax.set_title(f"Pressure profile along channel centre line\n(y = {y_mid:.2f})", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.35)

# Summary annotation box in the lower-left corner
ax.text(0.04, 0.06,
        f"Δp  Simulation :  {delta_p_sim:.3f}\n"
        f"Δp  Analytical  :  {delta_p_ana:.3f}\n"
        f"Deviation        :  {abs(delta_p_sim - delta_p_ana) / delta_p_ana * 100:.2f}%",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))

plt.tight_layout()

out_path = os.path.join(OUTPUT_DIR, "PlaneShearFlow_analysis.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nPlot saved to: {out_path}")
