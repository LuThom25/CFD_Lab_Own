#!/usr/bin/env python3
"""
Generate a PGM geometry file for the Airfoil Wind Tunnel case study.

Domain layout (DOMAIN_W x DOMAIN_H cells):
  - Left column  (col=0):          Inflow  (pixel=1)
  - Right column (col=DOMAIN_W-1): Outflow (pixel=2)
  - Top/bottom rows:               Wall    (pixel=3)
  - NACA 0012 profile:             Wall    (pixel=3), centered in domain
  - Everything else:               Fluid   (pixel=0)

Forbidden cells (obstacle cells with fluid on opposite sides, i.e. B_NS or
B_EW) are iteratively converted to fluid, producing a blunt trailing edge.
"""

import math
import os
import numpy as np

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
IMAX     = 200     # interior cells in x (matches imax in .dat file)
JMAX     = 80      # interior cells in y (matches jmax in .dat file)

# PGM dimensions include ghost/halo cells on all sides: (imax+2) x (jmax+2)
DOMAIN_W = IMAX + 2   # = 202
DOMAIN_H = JMAX + 2   # = 82

# CHORD=80: max half-thickness = 0.12*80/2 ≈ 4.8 → ~10 cells thick → clearly visible
CHORD    = 80      # airfoil chord length in interior grid cells
OFFSET_X = 40      # x-offset: 0.5c upstream of leading edge
OFFSET_Y = (JMAX // 2) + 1  # y-center in PGM row coords (mid of interior)
ALPHA_DEG = 0.0    # angle of attack in degrees

OUTPUT_DIR = os.path.join(os.path.dirname(__file__),
                          "..", "example_cases", "AirfoilWindTunnel")
OUTPUT_PGM = os.path.join(OUTPUT_DIR, "AirfoilWindTunnel.pgm")

# PGM pixel values
FLUID   = 0
INFLOW  = 1
OUTFLOW = 2
WALL    = 3


# ---------------------------------------------------------------------------
# NACA 0012 thickness distribution
# ---------------------------------------------------------------------------
def naca0012_half_thickness(x_norm: float) -> float:
    """Return half-thickness y_t at normalised chord position x in [0, 1]."""
    if x_norm < 0.0 or x_norm > 1.0:
        return 0.0
    x = x_norm
    return 5.0 * 0.12 * (
        0.2969 * math.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1015 * x**4
    )


def rotate(dx: float, dy: float, alpha_rad: float):
    """Rotate vector (dx, dy) by -alpha_rad (move from wind to body frame)."""
    ca, sa = math.cos(alpha_rad), math.sin(alpha_rad)
    return dx * ca + dy * sa, -dx * sa + dy * ca


# ---------------------------------------------------------------------------
# Build PGM array  (row 0 = physical bottom, row DOMAIN_H-1 = physical top)
# ---------------------------------------------------------------------------
def build_pgm() -> np.ndarray:
    pgm = np.zeros((DOMAIN_H, DOMAIN_W), dtype=np.int32)

    # Domain boundaries — matching ChannelWithObstacle.pgm convention:
    # Ghost rows (row 0, row DOMAIN_H-1): physical top/bottom walls.
    # Ghost columns (col 0, col DOMAIN_W-1): Inflow and Outflow.
    pgm[0, :]    = WALL     # bottom ghost row (physical bottom wall)
    pgm[-1, :]   = WALL     # top ghost row    (physical top wall)
    pgm[:, 0]    = INFLOW   # left ghost column
    pgm[:, -1]   = OUTFLOW  # right ghost column

    alpha_rad = math.radians(ALPHA_DEG)

    # Rasterise NACA 0012 in interior cells: rows 1..DOMAIN_H-2, cols 1..DOMAIN_W-2
    for row in range(1, DOMAIN_H - 1):
        for col in range(1, DOMAIN_W - 1):
            # Cell centre relative to leading edge.
            # Interior cols: col=1..DOMAIN_W-2 (cols 0 and DOMAIN_W-1 are ghost).
            dx = col - OFFSET_X + 0.5
            dy = row - OFFSET_Y + 0.5

            # Rotate into body frame
            bx, by = rotate(dx, dy, alpha_rad)

            x_norm = bx / CHORD
            y_norm = by / CHORD

            if 0.0 <= x_norm <= 1.0:
                yt = naca0012_half_thickness(x_norm)
                if -yt <= y_norm <= yt:
                    pgm[row, col] = WALL

    return pgm


# ---------------------------------------------------------------------------
# Forbidden-cell elimination
# ---------------------------------------------------------------------------
def eliminate_forbidden_cells(pgm: np.ndarray) -> int:
    """
    Convert forbidden obstacle cells to fluid in-place.
    A forbidden cell is a wall cell that has fluid on both opposite sides
    (top+bottom or left+right) or has 3+ fluid neighbours.
    Returns the number of cells converted.
    """
    H, W = pgm.shape
    total_converted = 0

    while True:
        converted = 0
        for row in range(1, H - 1):
            for col in range(1, W - 1):
                if pgm[row, col] != WALL:
                    continue

                n_fluid = sum([
                    pgm[row - 1, col] == FLUID,  # below  (physical south)
                    pgm[row + 1, col] == FLUID,  # above  (physical north)
                    pgm[row, col - 1] == FLUID,  # left
                    pgm[row, col + 1] == FLUID,  # right
                ])

                fluid_ns = (pgm[row - 1, col] == FLUID and
                            pgm[row + 1, col] == FLUID)
                fluid_ew = (pgm[row, col - 1] == FLUID and
                            pgm[row, col + 1] == FLUID)

                if fluid_ns or fluid_ew or n_fluid >= 3:
                    pgm[row, col] = FLUID
                    converted += 1

        total_converted += converted
        if converted == 0:
            break

    return total_converted


# ---------------------------------------------------------------------------
# Write PGM (P2 ASCII)
# ---------------------------------------------------------------------------
def write_pgm(pgm: np.ndarray, path: str) -> None:
    H, W = pgm.shape
    maxval = max(INFLOW, OUTFLOW, WALL)
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w") as f:
        f.write("P2\n")
        f.write(f"# NACA0012 airfoil, chord={CHORD}, alpha={ALPHA_DEG:.1f}deg\n")
        f.write(f"{W} {H}\n")
        f.write(f"{maxval}\n")
        # PGM row 0 = top of image; solver reads bottom-to-top, so we write
        # physical top row first (row index DOMAIN_H-1).
        for row in range(H - 1, -1, -1):
            f.write(" ".join(str(pgm[row, col]) for col in range(W)))
            f.write("\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"Building domain {DOMAIN_W}x{DOMAIN_H}, chord={CHORD} cells, "
          f"alpha={ALPHA_DEG}deg")

    pgm = build_pgm()

    wall_before = int(np.sum(pgm == WALL))
    n_fixed = eliminate_forbidden_cells(pgm)
    wall_after = int(np.sum(pgm == WALL))
    print(f"Forbidden-cell elimination: {n_fixed} cells converted to fluid "
          f"({wall_before} → {wall_after} wall cells)")

    write_pgm(pgm, OUTPUT_PGM)
    print(f"Written: {OUTPUT_PGM}")

    # Quick ASCII preview (80-char wide terminal)
    print("\nGeometry preview (. = fluid, # = wall, > = inflow, < = outflow):")
    scale = max(1, DOMAIN_H // 20)
    for row in range(DOMAIN_H - 1, -1, -scale):
        line = ""
        for col in range(0, DOMAIN_W, max(1, DOMAIN_W // 80)):
            v = pgm[row, col]
            line += {FLUID: ".", INFLOW: ">", OUTFLOW: "<", WALL: "#"}.get(v, "?")
        print(line)


if __name__ == "__main__":
    main()
