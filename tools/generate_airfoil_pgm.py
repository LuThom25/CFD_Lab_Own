#!/usr/bin/env python3
"""
Generate a PGM geometry file for the Airfoil Wind Tunnel case study.

Usage:
    python3 tools/generate_airfoil_pgm.py [--alpha DEGREES]

    --alpha: angle of attack in degrees (default 0).
             Output filename: AirfoilWindTunnel.pgm        for alpha = 0
                              AirfoilWindTunnel_alpha5.pgm  for alpha = 5, etc.

The script rasterises a NACA 0012 airfoil onto a rectangular grid, removes
forbidden cells at the thin trailing edge, and writes the result as a P2 PGM
file that the CFD solver reads as its geometry input.

Physical setup
--------------
  Free-stream:  U_inf = 1.0  (left to right, +x direction)
  Chord:        c = 1.6  (physical units, = 80 cells * dx)
  Reynolds:     Re = U_inf * c / nu = 1.0 * 1.6 / 0.008 = 200
  Domain:       4.0 x 1.6 physical  (200 x 80 interior cells)
"""

import argparse
import math
import os

import numpy as np

# ---------------------------------------------------------------------------
# Domain parameters
# ---------------------------------------------------------------------------

# imax / jmax: number of interior cells in x and y.
# These values are chosen so that:
#   - the chord spans exactly CHORD = 80 cells
#   - grid spacing dx = dy = xlength / imax = 4.0 / 200 = 0.02  (square cells)
#   - physical chord: 80 * 0.02 = 1.6  -> Re = 1.0 * 1.6 / 0.008 = 200
#   - ~10 cells across the airfoil maximum thickness (12% of chord = 9.6 cells)
#     which is the minimum to resolve the boundary-layer profile at Re = 200
IMAX = 200
JMAX = 80

# PGM dimensions include one ghost/halo cell on every side:
#   columns 0 and DOMAIN_W-1 are Inflow / Outflow ghost cells
#   rows    0 and DOMAIN_H-1 are top and bottom wall ghost cells
DOMAIN_W = IMAX + 2   # = 202
DOMAIN_H = JMAX + 2   # = 82

# CHORD: airfoil chord length measured in interior grid cells.
# 80 cells chosen so the profile has enough detail:
#   max half-thickness = 6% * 80 = 4.8 cells  -> ~10 cells full thickness
#   too few cells would staircase the curved surface too coarsely
CHORD = 80

# OFFSET_X: x-position (interior column index) of the leading edge.
# 40 cells = 0.5 * chord upstream of the LE gives the flow enough space to
# develop a uniform profile before reaching the airfoil.
OFFSET_X = 40

# OFFSET_Y: y-position (interior row index) of the airfoil centreline.
# JMAX // 2 + 1 = 41 places the LE exactly at mid-height of the domain.
# The +1 shifts the centreline one cell above the true midpoint so that the
# symmetric profile (which spans +/-yt cells around the chord) is centred
# rather than offset by half a cell.
OFFSET_Y = (JMAX // 2) + 1   # = 41

OUTPUT_DIR = os.path.join(os.path.dirname(__file__),
                          '..', 'example_cases', 'AirfoilWindTunnel')

# PGM pixel values — must match the solver's geometry reader
FLUID   = 0   # interior flow cell
INFLOW  = 1   # left boundary: Dirichlet velocity UIN, VIN
OUTFLOW = 2   # right boundary: Neumann outflow (dU/dx = 0, p = 0)
WALL    = 3   # solid no-slip wall (airfoil surface and channel walls)


# ---------------------------------------------------------------------------
# NACA 0012 thickness distribution
# ---------------------------------------------------------------------------

def naca0012_half_thickness(x_norm: float) -> float:
    """
    Return the local half-thickness y_t at normalised chord position x in [0, 1].

    The NACA 4-digit formula for a symmetric profile with 12% thickness ratio:

        y_t(x) = 5 * t * ( 0.2969*sqrt(x)
                          - 0.1260*x
                          - 0.3516*x^2
                          + 0.2843*x^3
                          - 0.1015*x^4 )

    where t = 0.12 (12% thickness) and x is the fractional chord position.
    The factor 5 * t = 0.6 scales the polynomial to exactly the desired
    maximum thickness.

    Upper surface: +y_t(x),  lower surface: -y_t(x).
    At x = 0 (leading edge):  y_t -> 0  (sharp leading edge in formula,
                                          but the sqrt term gives a finite
                                          leading-edge radius when evaluated
                                          over a finite cell size).
    At x = 1 (trailing edge): y_t ≈ 0.00126 * chord  (not exactly zero —
                                          the standard NACA formula has a
                                          small open trailing edge).
    """
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


# ---------------------------------------------------------------------------
# Rotation helper
# ---------------------------------------------------------------------------

def rotate_to_body_frame(dx: float, dy: float, alpha_rad: float):
    """
    Rotate a displacement vector (dx, dy) from the wind/grid frame into the
    airfoil body frame by applying a clockwise rotation of alpha_rad.

    Wind frame:  x points downstream (+right), y points up.
    Body frame:  x points along the chord from LE to TE,
                 y points perpendicular to the chord (upward from upper surface).

    Rotation matrix for clockwise rotation by alpha:
        [ bx ]   [  cos(alpha)   sin(alpha) ] [ dx ]
        [ by ] = [ -sin(alpha)   cos(alpha) ] [ dy ]

    With alpha > 0 (nose up), the chord tilts upward and the body-frame x-axis
    rotates clockwise relative to the wind-frame x-axis.
    """
    ca = math.cos(alpha_rad)
    sa = math.sin(alpha_rad)
    return dx * ca + dy * sa, -dx * sa + dy * ca


# ---------------------------------------------------------------------------
# Build PGM array
# ---------------------------------------------------------------------------

def build_pgm(alpha_deg: float = 0.0) -> np.ndarray:
    """
    Rasterise the NACA 0012 airfoil at angle of attack alpha_deg onto the grid.

    Coordinate convention (build_pgm / solver convention):
        row 0              = physical BOTTOM (south wall)
        row DOMAIN_H - 1   = physical TOP    (north wall)
        col 0              = Inflow ghost column  (west)
        col DOMAIN_W - 1   = Outflow ghost column (east)

    A cell at (row, col) is marked WALL if its centre lies inside the airfoil.
    The cell centre in grid coordinates is (col + 0.5, row + 0.5), shifted by
    the leading-edge offset so the LE is at (OFFSET_X, OFFSET_Y).

    Returns a 2-D integer array of shape (DOMAIN_H, DOMAIN_W).
    """
    pgm = np.zeros((DOMAIN_H, DOMAIN_W), dtype=np.int32)

    # Channel walls (top and bottom ghost rows)
    pgm[0, :]  = WALL
    pgm[-1, :] = WALL

    # Inflow / outflow ghost columns
    pgm[:, 0]  = INFLOW
    pgm[:, -1] = OUTFLOW

    alpha_rad = math.radians(alpha_deg)

    for row in range(1, DOMAIN_H - 1):      # loop over interior rows
        for col in range(1, DOMAIN_W - 1):  # loop over interior columns

            # Displacement of cell centre from the leading edge (in grid units)
            dx = col - OFFSET_X + 0.5
            dy = row - OFFSET_Y + 0.5

            # Transform to airfoil body frame
            # bx: distance along chord from LE (0 at LE, CHORD at TE)
            # by: distance perpendicular to chord (positive = upper surface)
            bx, by = rotate_to_body_frame(dx, dy, alpha_rad)

            # Normalised chord position
            x_norm = bx / CHORD

            if 0.0 <= x_norm <= 1.0:
                # Half-thickness at this chord station
                yt = naca0012_half_thickness(x_norm)

                # Cell is inside the airfoil if its perpendicular distance
                # from the chord is smaller than the local half-thickness
                y_norm = by / CHORD
                if -yt <= y_norm <= yt:
                    pgm[row, col] = WALL

    return pgm


# ---------------------------------------------------------------------------
# Forbidden-cell elimination
# ---------------------------------------------------------------------------

def eliminate_forbidden_cells(pgm: np.ndarray) -> int:
    """
    Iteratively convert forbidden wall cells to fluid and return the count.

    Why are forbidden cells a problem?
    -----------------------------------
    The staggered-grid solver assigns boundary conditions to wall cells based
    on which of their four neighbours are fluid.  It can handle:

        B_N  (fluid only above)     B_S  (fluid only below)
        B_E  (fluid only right)     B_W  (fluid only left)
        B_NE, B_NW, B_SE, B_SW     (two adjacent fluid neighbours)

    It CANNOT handle:
        B_NS  — fluid above AND below  ("tunnel" through the wall)
        B_EW  — fluid left  AND right  ("tunnel" through the wall)

    These configurations appear where the airfoil is thinner than one cell
    (typically the last few cells near the trailing edge).  A wall cell with
    3 or more fluid neighbours is equally ill-defined.

    The fix is to convert such cells to fluid, effectively creating a blunt
    trailing edge.  This is physically acceptable because:
      - The trailing edge of NACA 0012 has near-zero thickness anyway
      - Blunt trailing edges are common in practical CFD meshes
      - The affected region is limited to 3–5 cells near the very tip

    The elimination is repeated until no forbidden cell remains (convergence
    typically in 1–2 iterations).
    """
    H, W = pgm.shape
    total_converted = 0

    while True:
        converted = 0
        for row in range(1, H - 1):
            for col in range(1, W - 1):
                if pgm[row, col] != WALL:
                    continue

                # Count fluid neighbours
                fluid_n = pgm[row + 1, col] == FLUID   # north (above)
                fluid_s = pgm[row - 1, col] == FLUID   # south (below)
                fluid_e = pgm[row, col + 1] == FLUID   # east  (right)
                fluid_w = pgm[row, col - 1] == FLUID   # west  (left)
                n_fluid = fluid_n + fluid_s + fluid_e + fluid_w

                # Forbidden: fluid on both sides of a single axis, or ≥ 3 fluid
                if (fluid_n and fluid_s) or (fluid_e and fluid_w) or n_fluid >= 3:
                    pgm[row, col] = FLUID
                    converted += 1

        total_converted += converted
        if converted == 0:
            break   # no more forbidden cells -> converged

    return total_converted


# ---------------------------------------------------------------------------
# Write PGM file
# ---------------------------------------------------------------------------

def write_pgm(pgm: np.ndarray, path: str, alpha_deg: float = 0.0) -> None:
    """
    Write the geometry array as a P2 (ASCII greyscale) PGM file.

    PGM format:
        Line 1:  "P2"
        Line 2:  comment (optional)
        Line 3:  "<width> <height>"  (columns x rows of the image)
        Line 4:  "<maxval>"          (maximum pixel value)
        Lines 5+: pixel values row by row, left to right

    Coordinate flip:
        In the PGM file, row 0 of the image is the TOP of the picture.
        In the solver / build_pgm convention, row 0 is the BOTTOM.
        Therefore we write rows from index DOMAIN_H-1 down to 0 so that
        the solver reads the file bottom-to-top and reconstructs the
        original physical layout.
    """
    H, W = pgm.shape
    maxval = max(INFLOW, OUTFLOW, WALL)   # = 3
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, 'w') as f:
        f.write('P2\n')
        f.write(f'# NACA 0012 airfoil  chord={CHORD} cells  alpha={alpha_deg:.1f} deg\n')
        f.write(f'{W} {H}\n')
        f.write(f'{maxval}\n')
        # Write physical top row first (row DOMAIN_H-1) down to bottom row (0)
        for row in range(H - 1, -1, -1):
            f.write(' '.join(str(pgm[row, col]) for col in range(W)))
            f.write('\n')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Generate a NACA 0012 PGM geometry file for the CFD solver.')
    parser.add_argument('--alpha', type=float, default=0.0,
                        help='Angle of attack in degrees (default: 0)')
    args = parser.parse_args()
    alpha_deg = args.alpha

    # Output filename
    if alpha_deg == 0.0:
        output_pgm = os.path.join(OUTPUT_DIR, 'AirfoilWindTunnel.pgm')
    else:
        tag = f'alpha{int(round(alpha_deg))}'
        output_pgm = os.path.join(OUTPUT_DIR, f'AirfoilWindTunnel_{tag}.pgm')

    print(f'Domain:  {DOMAIN_W} x {DOMAIN_H} cells  '
          f'(interior: {IMAX} x {JMAX})')
    print(f'Chord:   {CHORD} cells  '
          f'(physical: {CHORD * 4.0 / IMAX:.3f})')
    print(f'Alpha:   {alpha_deg} deg')
    print(f'Re:      {1.0 * CHORD * 4.0 / IMAX / 0.008:.0f}  '
          f'(U_inf=1.0, c={CHORD * 4.0 / IMAX:.3f}, nu=0.008)')

    # Build geometry
    pgm = build_pgm(alpha_deg)
    wall_before = int(np.sum(pgm == WALL))

    # Remove forbidden cells
    n_fixed = eliminate_forbidden_cells(pgm)
    wall_after = int(np.sum(pgm == WALL))
    print(f'Forbidden cells eliminated: {n_fixed}  '
          f'({wall_before} -> {wall_after} wall cells)')

    # Write file
    write_pgm(pgm, output_pgm, alpha_deg)
    print(f'Written: {output_pgm}')

    # ASCII preview in terminal
    print('\nGeometry preview  (. fluid  # wall  > inflow  < outflow):')
    scale = max(1, DOMAIN_H // 20)
    for row in range(DOMAIN_H - 1, -1, -scale):
        line = ''
        for col in range(0, DOMAIN_W, max(1, DOMAIN_W // 80)):
            v = pgm[row, col]
            line += {FLUID: '.', INFLOW: '>', OUTFLOW: '<', WALL: '#'}.get(v, '?')
        print(line)


if __name__ == '__main__':
    main()
