#!/usr/bin/env python3
"""
Evaluate drag and lift coefficients for the Airfoil Wind Tunnel case.

Usage:
    python3 tools/eval_airfoil.py <output_dir>

    <output_dir>: path to the case output directory, e.g.
        example_cases/AirfoilWindTunnel/AirfoilWindTunnel_Re200_1_1_Output

The script reads the last VTK file in that directory and computes CD and CL
via a surface integral over all obstacle cells adjacent to the fluid.

Physical parameters (must match the .dat file):
    U_inf = 1.0, c = 1.0 (chord), rho = 1.0, nu = 0.005 (Re=200)
"""

import os
import sys
import glob
import re

import numpy as np

# ---------------------------------------------------------------------------
# Physical parameters — keep in sync with the .dat file
# ---------------------------------------------------------------------------
U_INF  = 1.0    # free-stream velocity
C      = 1.0    # chord length (physical)
RHO    = 1.0    # fluid density (incompressible, normalised)
NU     = 0.005  # kinematic viscosity  (Re = U*C/NU = 200)

# Reference drag coefficient at Re=200, NACA0012, alpha=0
CD_REF = 0.28   # Liu et al. 1998

# Dynamic pressure
Q_INF  = 0.5 * RHO * U_INF**2 * C  # = 0.5


# ---------------------------------------------------------------------------
# VTK reader (structured grid, ASCII)
# ---------------------------------------------------------------------------
def read_last_vtk(output_dir: str):
    """Return (imax, jmax, pressure_2d, u_2d, v_2d, obstacle_2d) from the
    last VTK file in output_dir.  Arrays are indexed [i, j] (0-based,
    excluding halo), matching the internal grid layout."""

    pattern = os.path.join(output_dir, "*.vtk")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No VTK files found in {output_dir}")

    vtk_file = files[-1]
    print(f"Reading: {os.path.basename(vtk_file)}")

    with open(vtk_file) as f:
        lines = f.readlines()

    # Parse DIMENSIONS line → imax+1, jmax+1, 1
    dims = None
    for i, line in enumerate(lines):
        if line.startswith("DIMENSIONS"):
            parts = line.split()
            dims = int(parts[1]) - 1, int(parts[2]) - 1  # (imax, jmax)
            break
    if dims is None:
        raise ValueError("DIMENSIONS not found in VTK file")
    imax, jmax = dims

    def read_scalar(name: str) -> np.ndarray:
        for i, line in enumerate(lines):
            if name in line and "SCALARS" in line:
                # skip LOOKUP_TABLE line
                data_start = i + 2
                vals = []
                j = data_start
                while len(vals) < imax * jmax:
                    vals.extend(float(x) for x in lines[j].split())
                    j += 1
                arr = np.array(vals).reshape((jmax, imax))
                return arr  # [j, i]
        raise ValueError(f"Scalar '{name}' not found in VTK file")

    def read_vector(name: str) -> tuple:
        for i, line in enumerate(lines):
            if name in line and "VECTORS" in line:
                data_start = i + 1
                vals = []
                j = data_start
                while len(vals) < imax * jmax * 3:
                    vals.extend(float(x) for x in lines[j].split())
                    j += 1
                arr = np.array(vals).reshape((jmax, imax, 3))
                return arr[:, :, 0], arr[:, :, 1]  # u[j,i], v[j,i]
        raise ValueError(f"Vector '{name}' not found in VTK file")

    p_ji = read_scalar("pressure")
    obs_ji = read_scalar("obstacle")
    try:
        u_ji, v_ji = read_vector("velocity")
    except ValueError:
        u_ji = np.zeros((jmax, imax))
        v_ji = np.zeros((jmax, imax))

    return imax, jmax, p_ji, u_ji, v_ji, obs_ji


# ---------------------------------------------------------------------------
# Force integration
# ---------------------------------------------------------------------------
def compute_forces(imax, jmax, p_ji, u_ji, v_ji, obs_ji, dx, dy):
    """
    Integrate pressure drag and lift over airfoil surface cells.

    For each obstacle cell adjacent to a fluid cell, the outward normal points
    from the obstacle into the fluid.  The pressure force on the fluid side
    equals p * n * dA (per unit depth).  We accumulate the reaction force on
    the airfoil: F = -p * n * dA.

    Returns (F_drag, F_lift).
    """
    F_drag = 0.0
    F_lift = 0.0

    # Obstacle id=1 in VTK (FIXED_WALL), fluid id=0
    OBSTACLE_ID = 1.0

    for j in range(jmax):
        for i in range(imax):
            if obs_ji[j, i] != OBSTACLE_ID:
                continue  # not an obstacle cell

            # Check each face for fluid neighbour
            # Right face (east): normal = (+1, 0), dA = dy
            if i + 1 < imax and obs_ji[j, i + 1] == 0.0:
                p_face = 0.5 * (p_ji[j, i] + p_ji[j, i + 1])
                F_drag += p_face * dy   # +x → drag
            # Left face (west): normal = (-1, 0), dA = dy
            if i - 1 >= 0 and obs_ji[j, i - 1] == 0.0:
                p_face = 0.5 * (p_ji[j, i] + p_ji[j, i - 1])
                F_drag -= p_face * dy   # -x → negative drag
            # Top face (north): normal = (0, +1), dA = dx
            if j + 1 < jmax and obs_ji[j + 1, i] == 0.0:
                p_face = 0.5 * (p_ji[j, i] + p_ji[j + 1, i])
                F_lift += p_face * dx   # +y → lift
            # Bottom face (south): normal = (0, -1), dA = dx
            if j - 1 >= 0 and obs_ji[j - 1, i] == 0.0:
                p_face = 0.5 * (p_ji[j, i] + p_ji[j - 1, i])
                F_lift -= p_face * dx   # -y → negative lift

    return F_drag, F_lift


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    output_dir = sys.argv[1]
    if not os.path.isdir(output_dir):
        print(f"Error: directory not found: {output_dir}")
        sys.exit(1)

    imax, jmax, p_ji, u_ji, v_ji, obs_ji = read_last_vtk(output_dir)

    # Grid spacing (from .dat: xlength=4.0, ylength=2.0)
    xlength, ylength = 4.0, 2.0
    dx = xlength / imax
    dy = ylength / jmax

    print(f"Grid: {imax} x {jmax},  dx={dx:.4f},  dy={dy:.4f}")
    print(f"Obstacle cells: {int(np.sum(obs_ji == 1.0))}")

    F_drag, F_lift = compute_forces(imax, jmax, p_ji, u_ji, v_ji, obs_ji, dx, dy)

    CD = F_drag / Q_INF
    CL = F_lift / Q_INF

    print()
    print(f"  F_drag = {F_drag:.6f} N/m")
    print(f"  F_lift = {F_lift:.6f} N/m")
    print()
    print(f"  CD = {CD:.4f}   (reference Re=200: {CD_REF:.2f})")
    print(f"  CL = {CL:.4f}   (expected ~0 for alpha=0)")
    print()
    if abs(CD) > 0:
        err = abs(CD - CD_REF) / CD_REF * 100
        print(f"  CD error vs literature: {err:.1f}%")


if __name__ == "__main__":
    main()
