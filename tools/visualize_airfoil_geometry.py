#!/usr/bin/env python3
"""
Visualize NACA 0012 airfoil geometry for the Airfoil Wind Tunnel case study.

Generates three PNG files in example_cases/AirfoilWindTunnel/:
  geometry_1_full.png  – full domain overview for all 4 angle-of-attack cases
  geometry_2_zoom.png  – zoomed detail view with analytical NACA contour overlay
  geometry_3_all.png   – combined figure (full domain top, detail view bottom)

Each plot shows:
  - Rasterized wall cells (purple) — the actual discrete geometry used by the solver
  - Eliminated forbidden cells (yellow) — cells removed because they would cause
    ill-defined boundary conditions (fluid on both opposite sides)
  - Analytical NACA 0012 formula contour (red) — the smooth mathematical profile
    overlaid on the rasterized grid to show the discretization error
  - Leading Edge (LE, green circle) and Trailing Edge (TE, orange triangle) markers
  - Chord line (dashed) and angle-of-attack arc

Usage:
    python3 tools/visualize_airfoil_geometry.py

Output is written to:
    example_cases/AirfoilWindTunnel/geometry_1_full.png
    example_cases/AirfoilWindTunnel/geometry_2_zoom.png
    example_cases/AirfoilWindTunnel/geometry_3_all.png
"""

import math
import os

import matplotlib
matplotlib.use('Agg')  # non-interactive backend, no display needed
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Arc
from matplotlib.colors import ListedColormap
import numpy as np

# ---------------------------------------------------------------------------
# Domain parameters  (must match generate_airfoil_pgm.py and .dat files)
# ---------------------------------------------------------------------------
IMAX     = 200    # interior cells in x
JMAX     = 80     # interior cells in y
CHORD    = 80     # airfoil chord in grid cells  (= 80 * dx = 80 * 0.02 = 1.6 physical)
OFFSET_X = 40     # LE x-position in grid cells  (= 0.5 * chord upstream)
OFFSET_Y = 41     # LE y-position in grid cells  (vertical domain center)
DOMAIN_W = IMAX + 2
DOMAIN_H = JMAX + 2

# PGM pixel values
FLUID   = 0
INFLOW  = 1
OUTFLOW = 2
WALL    = 3

# Marker colors
COL_LE   = '#27AE60'   # green  — leading edge
COL_TE   = '#E67E22'   # orange — trailing edge
COL_NACA = '#C0392B'   # dark red — analytical NACA contour

# Angles of attack to visualize
ALPHAS = [0, 5, 10, 15]

OUTPUT_DIR = os.path.join(os.path.dirname(__file__),
                          '..', 'example_cases', 'AirfoilWindTunnel')


# ---------------------------------------------------------------------------
# NACA 0012 thickness distribution
# ---------------------------------------------------------------------------
def naca0012_ht_scalar(x: float) -> float:
    """Half-thickness at normalised chord position x in [0, 1]."""
    if x < 0 or x > 1:
        return 0.0
    return 5.0 * 0.12 * (
        0.2969 * math.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1015 * x**4
    )


def naca0012_ht_array(x: np.ndarray) -> np.ndarray:
    """Vectorised version for smooth curve plotting."""
    x = np.maximum(x, 0.0)
    return 5.0 * 0.12 * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1015 * x**4
    )


# ---------------------------------------------------------------------------
# Rebuild PGM in memory (mirrors generate_airfoil_pgm.py logic)
# ---------------------------------------------------------------------------
def build_pgm(alpha_deg: float) -> np.ndarray:
    """
    Rasterise NACA 0012 at angle of attack alpha_deg into a 2-D integer array.
    Row 0 = physical bottom (y = 0), row DOMAIN_H-1 = physical top.
    This is the native build_pgm coordinate system used throughout this script.
    """
    pgm = np.zeros((DOMAIN_H, DOMAIN_W), dtype=np.int32)
    pgm[0, :]  = WALL     # bottom wall
    pgm[-1, :] = WALL     # top wall
    pgm[:, 0]  = INFLOW   # left ghost column
    pgm[:, -1] = OUTFLOW  # right ghost column

    a  = math.radians(alpha_deg)
    ca = math.cos(a)
    sa = math.sin(a)

    for row in range(1, DOMAIN_H - 1):
        for col in range(1, DOMAIN_W - 1):
            # Cell centre displacement from leading edge
            dx = col - OFFSET_X + 0.5
            dy = row - OFFSET_Y + 0.5
            # Rotate into body (airfoil) frame: clockwise by alpha
            bx =  dx * ca + dy * sa
            by = -dx * sa + dy * ca
            xn = bx / CHORD
            if 0.0 <= xn <= 1.0:
                yt = naca0012_ht_scalar(xn)
                if abs(by / CHORD) <= yt:
                    pgm[row, col] = WALL

    return pgm


def eliminate_forbidden_cells(pgm: np.ndarray):
    """
    Remove forbidden obstacle cells iteratively.

    A cell is forbidden when:
      - it has fluid directly above AND below  (B_NS tunnel)
      - it has fluid directly left AND right   (B_EW tunnel)
      - it has 3 or more fluid neighbours

    Such cells cannot have well-defined boundary conditions in the staggered-
    grid solver, so they are converted to fluid.  This produces a blunt
    trailing edge near the thin airfoil tip, which is physically acceptable.

    Returns:
        pgm_clean    cleaned copy of the input array
        eliminated   set of (row, col) tuples that were converted to fluid
    """
    pgm = pgm.copy()
    eliminated = set()
    H, W = pgm.shape
    while True:
        converted = 0
        for row in range(1, H - 1):
            for col in range(1, W - 1):
                if pgm[row, col] != WALL:
                    continue
                fluid_n = pgm[row + 1, col] == FLUID
                fluid_s = pgm[row - 1, col] == FLUID
                fluid_e = pgm[row, col + 1] == FLUID
                fluid_w = pgm[row, col - 1] == FLUID
                n_fluid = fluid_n + fluid_s + fluid_e + fluid_w
                if (fluid_n and fluid_s) or (fluid_e and fluid_w) or n_fluid >= 3:
                    pgm[row, col] = FLUID
                    eliminated.add((row, col))
                    converted += 1
        if converted == 0:
            break
    return pgm, eliminated


# ---------------------------------------------------------------------------
# Analytical NACA contour in display coordinates
# ---------------------------------------------------------------------------
def smooth_contour(alpha_deg: float):
    """
    Return (upper, lower) contour arrays in grid display coordinates.

    The analytical profile is computed in body frame, then rotated back to
    the grid frame (inverse of the rasterisation rotation).

    Coordinate convention: display x = grid column, display y = grid row
    with row 0 at the physical bottom (same as build_pgm).
    """
    a  = math.radians(alpha_deg)
    ca = math.cos(a)
    sa = math.sin(a)
    xs = np.linspace(0.0, 1.0, 600)
    yt = naca0012_ht_array(xs)
    contours = []
    for sign in [+1.0, -1.0]:   # upper (+), lower (-)
        bx = xs * CHORD
        by = sign * yt * CHORD
        # Inverse rotation: counter-clockwise by alpha
        dx = bx * ca - by * sa
        dy = bx * sa + by * ca
        contours.append((OFFSET_X + dx, OFFSET_Y + dy))
    return contours[0], contours[1]


# ---------------------------------------------------------------------------
# LE / TE detection from rasterised cells
# ---------------------------------------------------------------------------
def find_le_te(pgm: np.ndarray, min_thickness: int = 2):
    """
    Locate leading and trailing edge from rasterised wall cells.

    LE = leftmost column with at least min_thickness wall cells (median y).
    TE = rightmost column with at least min_thickness wall cells (median y).

    Returns (le, te, x_min, x_max, y_min, y_max) all in display coordinates
    (column = display x, row = display y, row 0 = bottom).
    """
    interior = pgm[1:-1, 1:-1]
    rows, cols = np.where(interior == WALL)
    disp_x = cols + 1   # interior col → pgm col = display x
    disp_y = rows + 1   # interior row → pgm row = display y

    col_min, col_max = int(cols.min()), int(cols.max())

    le = te = None
    for c in range(col_min, col_max + 1):
        mask = cols == c
        if mask.sum() >= min_thickness:
            le = (float(c + 1), float(np.median(disp_y[mask])))
            break
    for c in range(col_max, col_min - 1, -1):
        mask = cols == c
        if mask.sum() >= min_thickness:
            te = (float(c + 1), float(np.median(disp_y[mask])))
            break

    return (le, te,
            int(disp_x.min()), int(disp_x.max()),
            int(disp_y.min()), int(disp_y.max()))


# ---------------------------------------------------------------------------
# Pre-compute all data
# ---------------------------------------------------------------------------
def precompute():
    results = {}
    for alpha in ALPHAS:
        raw            = build_pgm(alpha)
        clean, elim    = eliminate_forbidden_cells(raw)
        pgm_disp       = clean.copy()
        for (r, c) in elim:
            pgm_disp[r, c] = 4          # value 4 = eliminated (shown in yellow)
        le, te, xmin, xmax, ymin, ymax = find_le_te(clean)
        upper, lower   = smooth_contour(alpha)
        results[alpha] = dict(
            pgm_disp=pgm_disp,
            le=le, te=te,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            upper=upper, lower=lower,
            n_elim=len(elim),
        )
        print(f'  α={alpha:2d}°: LE=({le[0]:.0f},{le[1]:.1f})  '
              f'TE=({te[0]:.0f},{te[1]:.1f})  '
              f'cells x={xmin}–{xmax} y={ymin}–{ymax}  '
              f'forbidden={len(elim)}')
    return results


# ---------------------------------------------------------------------------
# Shared colour map  (values 0–4)
# ---------------------------------------------------------------------------
CMAP = ListedColormap([
    '#D6EAF8',   # 0 = fluid        (light blue)
    '#D5F5E3',   # 1 = inflow       (light green)
    '#FDEBD0',   # 2 = outflow      (light orange)
    '#7D3C98',   # 3 = wall         (purple)
    '#FFD700',   # 4 = eliminated   (yellow)
])


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _add_le_te(ax, le, te, labels=True, ms=8):
    ax.plot(le[0], le[1], 'o', color=COL_LE, ms=ms, zorder=9, mec='white', mew=0.6)
    ax.plot(te[0], te[1], '^', color=COL_TE, ms=ms, zorder=9, mec='white', mew=0.6)
    if labels:
        ax.text(le[0] - 1.5, le[1] + 1.5, 'LE', fontsize=8, color=COL_LE,
                ha='right', va='bottom', fontweight='bold')
        ax.text(te[0] + 1.5, te[1] + 1.5, 'TE', fontsize=8, color=COL_TE,
                ha='left',  va='bottom', fontweight='bold')


def draw_full(ax, alpha, data):
    """Draw full-domain panel for one angle of attack."""
    d = data[alpha]
    ax.imshow(d['pgm_disp'], cmap=CMAP, vmin=0, vmax=4,
              interpolation='nearest', origin='lower', aspect='equal')
    _add_le_te(ax, d['le'], d['te'], labels=False, ms=7)
    ax.set_title(f'Angle of attack α = {alpha}°', fontsize=11, fontweight='bold')
    ax.set_xlabel('Grid cells in x-direction (i)', fontsize=8)
    ax.set_ylabel('Grid cells in y-direction (j)', fontsize=8)
    ax.tick_params(labelsize=7)


def draw_zoom(ax, alpha, data):
    """Draw zoomed detail panel with grid, analytical contour, LE/TE, angle arc."""
    d = data[alpha]
    ax.imshow(d['pgm_disp'], cmap=CMAP, vmin=0, vmax=4,
              interpolation='nearest', origin='lower', aspect='equal')

    xmin, xmax = d['xmin'], d['xmax']
    ymin, ymax = d['ymin'], d['ymax']
    pad_x = max(8, (xmax - xmin) // 5)
    pad_y = max(10, (ymax - ymin) // 2 + 6)
    ax.set_xlim(xmin - pad_x, xmax + pad_x)
    ax.set_ylim(ymin - pad_y, ymax + pad_y)

    # Cell-level grid
    x0, x1 = int(xmin - pad_x), int(xmax + pad_x)
    y0, y1 = int(ymin - pad_y), int(ymax + pad_y)
    ax.set_xticks(np.arange(x0, x1 + 1, 5))
    ax.set_xticks(np.arange(x0, x1 + 1, 1), minor=True)
    ax.set_yticks(np.arange(y0, y1 + 1, 5))
    ax.set_yticks(np.arange(y0, y1 + 1, 1), minor=True)
    ax.grid(which='minor', color='#bbb', lw=0.2, alpha=0.5)
    ax.grid(which='major', color='#888', lw=0.4, alpha=0.6)

    # Analytical NACA 0012 contour
    ax.plot(d['upper'][0], d['upper'][1], '-', color=COL_NACA, lw=1.7, zorder=7)
    ax.plot(d['lower'][0], d['lower'][1], '-', color=COL_NACA, lw=1.7, zorder=7)

    # Chord line and horizontal reference
    le, te = d['le'], d['te']
    ax.plot([le[0], te[0]], [le[1], te[1]], '--', color='#555', lw=1.0,
            alpha=0.8, zorder=5)
    ax.axhline(le[1], color='#aaa', lw=0.8, ls=':', alpha=0.7)

    # Angle-of-attack arc at leading edge
    if alpha > 0:
        r = max(6, (xmax - xmin) // 8)
        ax.add_patch(Arc((le[0], le[1]), 2 * r, 2 * r, angle=0,
                         theta1=0, theta2=alpha, color='#2980B9', lw=2.0, zorder=8))
        mid = math.radians(alpha / 2)
        ax.annotate(f'α = {alpha}°',
                    xy=(le[0] + (r + 2) * math.cos(mid),
                        le[1] + (r + 2) * math.sin(mid)),
                    fontsize=9, color='#1A5276', fontweight='bold')

    _add_le_te(ax, le, te, labels=True, ms=9)

    n = d['n_elim']
    ax.set_title(
        f'α = {alpha}° — Detail view  '
        f'({n} forbidden cell{"s" if n != 1 else ""} eliminated)',
        fontsize=10, fontweight='bold')
    ax.set_xlabel('Grid cells in x-direction (i)', fontsize=8)
    ax.set_ylabel('Grid cells in y-direction (j)', fontsize=8)
    ax.tick_params(labelsize=7)


# ---------------------------------------------------------------------------
# Shared legend entries
# ---------------------------------------------------------------------------
LEG_FULL = [
    mpatches.Patch(color='#D6EAF8', label='Fluid'),
    mpatches.Patch(color='#7D3C98', label='Wall cells (rasterized)'),
    mpatches.Patch(color='#FFD700', label='Eliminated forbidden cells'),
    plt.Line2D([0],[0], marker='o', color=COL_LE, ls='None', ms=9,
               label='Leading Edge (LE)'),
    plt.Line2D([0],[0], marker='^', color=COL_TE, ls='None', ms=9,
               label='Trailing Edge (TE)'),
]
LEG_ZOOM = LEG_FULL + [
    plt.Line2D([0],[0], color=COL_NACA, lw=1.8,
               label='NACA 0012 formula (analytical)'),
    plt.Line2D([0],[0], color='#555', ls='--', lw=1.0, label='Chord line'),
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print('Pre-computing geometry data...')
    data = precompute()

    # ── PNG 1: Full domain overview ──────────────────────────────────────
    print('Saving geometry_1_full.png ...')
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    fig.patch.set_facecolor('white')
    for ci, alpha in enumerate(ALPHAS):
        draw_full(axes[ci], alpha, data)
    fig.legend(handles=LEG_FULL, loc='lower center', ncol=5, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle(
        'NACA 0012 Airfoil Wind Tunnel — Full Domain Overview\n'
        r'Grid: $i_{max}=200$, $j_{max}=80$,  Reynolds number $Re = 200$',
        fontsize=13, fontweight='bold', y=1.04)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'geometry_1_full.png'),
                bbox_inches='tight', dpi=150)
    plt.close(fig)

    # ── PNG 2: Detail / zoomed view ──────────────────────────────────────
    print('Saving geometry_2_zoom.png ...')
    fig, axes = plt.subplots(1, 4, figsize=(22, 6))
    fig.patch.set_facecolor('white')
    for ci, alpha in enumerate(ALPHAS):
        draw_zoom(axes[ci], alpha, data)
    fig.legend(handles=LEG_ZOOM, loc='lower center', ncol=7, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle(
        'NACA 0012 Airfoil Wind Tunnel — Detail View with Analytical Profile Contour\n'
        'Purple = rasterized wall cells · Yellow = eliminated forbidden cells '
        '· Red = NACA formula',
        fontsize=12, fontweight='bold', y=1.04)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'geometry_2_zoom.png'),
                bbox_inches='tight', dpi=150)
    plt.close(fig)

    # ── PNG 3: Combined (full top, detail bottom) ────────────────────────
    print('Saving geometry_3_all.png ...')
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    fig.patch.set_facecolor('white')
    for ci, alpha in enumerate(ALPHAS):
        draw_full(axes[0, ci], alpha, data)
        draw_zoom(axes[1, ci], alpha, data)
    fig.legend(handles=LEG_ZOOM, loc='lower center', ncol=7, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle(
        'NACA 0012 Airfoil Wind Tunnel — Full Domain (top) & Detail View (bottom)\n'
        r'$i_{max}=200$, $j_{max}=80$,  chord length = 80 cells,  $Re = 200$',
        fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'geometry_3_all.png'),
                bbox_inches='tight', dpi=150)
    plt.close(fig)

    print('Done. Files written to:', os.path.abspath(OUTPUT_DIR))


if __name__ == '__main__':
    main()
