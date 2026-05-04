#!/usr/bin/env python3
"""
Worksheet 1 — Simulation Tasks 4–8
====================================
Runs the fluidchen solver with different parameter configurations and
presents formatted comparison tables for each task.

Usage (from this directory):
  python3 run_studies.py          # all tasks (4-8)
  python3 run_studies.py 4        # Task 4 only
  python3 run_studies.py 5 6 7 8  # explicit task list

Tasks:
  4 — Base LDC simulation: 50×50, nu=0.01 (Re=100), adaptive dt, t_end=50
  5 — SOR solver: effect of relaxation factor omega and itermax
  6 — Fixed time step: which dt values lead to a stable simulation?
  7 — Grid refinement: imax = 16, 32, 64, 128, 256 with fixed dt = 0.05
  8 — Viscosity study: nu = 0.01 → 0.0001 (Re = 100 → 10000), adaptive dt
"""

import subprocess, re, time, shutil, sys, tempfile, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")          # headless – no GUI needed
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PIL import Image, ImageDraw, ImageFont

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
BINARY     = SCRIPT_DIR.parent.parent / "build" / "fluidchen"

# pvpython for ParaView rendering — check all installed versions
_PV_CANDIDATES = [
    Path("/Applications/ParaView-6.0.1.app/Contents/bin/pvpython"),
    Path("/Applications/ParaView-6.1.0.app/Contents/bin/pvpython"),
    Path("/usr/local/bin/pvpython"),
    Path("/usr/bin/pvpython"),
]
PVPYTHON = next((p for p in _PV_CANDIDATES if p.exists()), None)
PLOTS_DIR  = SCRIPT_DIR / "LidDrivenCavity_Output" / "study_plots"

# ── Plot style ────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.grid": True,
    "grid.alpha": 0.35,
})

# ── PIL comparison grid (pixel-perfect, no matplotlib layout clipping) ────────
def _load_font(size):
    for p in ["/System/Library/Fonts/Helvetica.ttc",
              "/System/Library/Fonts/SFNS.ttf",
              "/Library/Fonts/Arial.ttf"]:
        if Path(p).exists():
            try:
                return ImageFont.truetype(p, size)
            except Exception:
                pass
    return ImageFont.load_default()

def make_comparison_grid(img_rows, col_titles, row_labels, suptitle, out_path):
    """
    Assemble a 2-D grid of images (PIL Paths) into a single PNG.
    img_rows : list[list[Path|None]]  — [row][col]
    Uses PIL directly → no matplotlib axes clipping, every pixel is correct.
    """
    # Determine cell size from first valid image
    cell_w, cell_h = 1200, 1000
    for row in img_rows:
        for p in row:
            if p and Path(p).exists():
                cell_w, cell_h = Image.open(p).size
                break
        else:
            continue
        break

    n_rows = len(img_rows)
    n_cols = max(len(r) for r in img_rows)

    PAD      = 10   # px gap between images
    SUP_H    = 64   # suptitle bar height
    COL_H    = 54   # column-title bar height
    ROW_LW   = 170  # row-label column width

    total_w = ROW_LW + n_cols * cell_w + (n_cols + 1) * PAD
    total_h = SUP_H + n_rows * (COL_H + cell_h) + (n_rows + 1) * PAD

    canvas = Image.new("RGB", (total_w, total_h), (245, 245, 245))
    draw   = ImageDraw.Draw(canvas)

    f_sup = _load_font(28)
    f_col = _load_font(20)
    f_row = _load_font(18)

    # ── Suptitle bar ──────────────────────────────────────────────────────────
    draw.rectangle([(0, 0), (total_w, SUP_H)], fill=(40, 40, 50))
    bb = draw.textbbox((0, 0), suptitle, font=f_sup)
    draw.text(((total_w - (bb[2]-bb[0]))//2, (SUP_H - (bb[3]-bb[1]))//2),
              suptitle, fill=(255, 255, 255), font=f_sup)

    # ── Column titles ─────────────────────────────────────────────────────────
    for c, ct in enumerate(col_titles):
        x0 = ROW_LW + PAD + c * (cell_w + PAD)
        draw.rectangle([(x0, SUP_H), (x0 + cell_w, SUP_H + COL_H)],
                       fill=(60, 60, 75))
        bb = draw.textbbox((0, 0), ct, font=f_col)
        draw.text((x0 + (cell_w - (bb[2]-bb[0]))//2,
                   SUP_H + (COL_H - (bb[3]-bb[1]))//2),
                  ct, fill=(220, 220, 220), font=f_col)

    # ── Row labels + images ───────────────────────────────────────────────────
    for r, row in enumerate(img_rows):
        y0 = SUP_H + r * (COL_H + cell_h + PAD) + COL_H + PAD

        # Row label column (dark sidebar)
        draw.rectangle([(0, y0), (ROW_LW, y0 + cell_h)], fill=(60, 60, 75))
        rl = row_labels[r] if r < len(row_labels) else ""
        bb = draw.textbbox((0, 0), rl, font=f_row)
        draw.text((( ROW_LW - (bb[2]-bb[0]))//2,
                   y0 + (cell_h - (bb[3]-bb[1]))//2),
                  rl, fill=(220, 220, 220), font=f_row)

        for c, p in enumerate(row):
            x0 = ROW_LW + PAD + c * (cell_w + PAD)
            if p and Path(p).exists():
                canvas.paste(Image.open(p), (x0, y0))
            else:
                draw.rectangle([(x0, y0), (x0+cell_w, y0+cell_h)],
                                fill=(200, 200, 200))
                draw.text((x0 + cell_w//2, y0 + cell_h//2), "N/A",
                          fill=(100, 100, 100), font=f_col, anchor="mm")

    canvas.save(out_path, dpi=(150, 150))
    return Image.open(out_path).size

# ── Parse the SUMMARY line written by Case::simulate() ────────────────────────
SUMMARY_RE = re.compile(
    r"SUMMARY\s+t=([\d.]+)\s+steps=(\d+)\s+vtk=(\d+)\s+"
    r"avg_sor=([\d.]+)\s+max_sor=(\d+)\s+avg_res=([\d.eE+\-]+)\s+"
    r"avg_dt=([\d.eE+\-]+)\s+(?:solver=\S+\s+)?status=(\w+)"
)

# ── Read solver selection from LidDrivenCavity.dat ────────────────────────────
# The .dat file contains a "solver" field that selects which pressure Poisson
# solver to use.  We read it here so that ALL tasks (4-8) automatically use
# whichever solver is set in the config file.
#
# Available solvers (set in LidDrivenCavity.dat):
#   solver  SOR_MEAN_CORRECTION - Thomas: SOR + Fredholm fix + zero-mean projection
#   solver  SOR_RB              - Red-Black SOR  (checkerboard, parallelisation-ready)
#   solver  SOR_STANDARD        - Iciar: plain SOR (no mean correction)
#   (omit line)                 - defaults to SOR_MEAN_CORRECTION
_DAT_FILE = SCRIPT_DIR / "LidDrivenCavity.dat"

def _read_solver_from_dat(dat_path: Path) -> str:
    """Parse the 'solver' field from a .dat file.  Returns 'SOR_MEAN_CORRECTION' if not found."""
    if not dat_path.exists():
        return "SOR_MEAN_CORRECTION"
    for line in dat_path.read_text().splitlines():
        stripped = line.strip()
        # Skip blank lines and comments
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) >= 2 and parts[0].lower() == "solver":
            return parts[1].upper()   # e.g. "SOR_RB"
    return "SOR_MEAN_CORRECTION"

_ACTIVE_SOLVER = _read_solver_from_dat(_DAT_FILE)
print(f"[run_studies] Solver read from LidDrivenCavity.dat: {_ACTIVE_SOLVER}")

# ── Base configuration (worksheet defaults) ───────────────────────────────────
BASE_CFG = dict(
    xlength=1.0, ylength=1.0,
    imax=50,     jmax=50,
    dt=0.05,     t_end=5.0,   tau=0.5,
    dt_value=999,              # suppress VTK output during studies (speed)
    eps=0.001,   omg=1.7,     gamma=0.5,  itermax=100,
    solver=_ACTIVE_SOLVER,     # forwarded from LidDrivenCavity.dat to every task
    # solver must NOT be the last entry — the C++ parser fails to read the value
    # when solver is at EOF with no tokens following it. Keep nu/GX/... after it.
    nu=0.01,
    GX=0.0,      GY=0.0,
    PI=0.0,      UI=0.0,      VI=0.0,
)

# ── Task-5a SOR omega study settings ─────────────────────────────────────────
# For the omega sweep we deliberately request a residual near double-precision
# roundoff. This is not a physically meaningful convergence target; it is used
# to force the SOR loop to terminate by itermax so omega values can be compared
# by residual decay after the same fixed iteration budget.
SOR_OMEGA_EPS = 1.0e-15
SOR_OMEGA_ITERMAX = 100

# ── Utility: write a .dat file ────────────────────────────────────────────────
def make_dat(cfg: dict, path: Path) -> None:
    with open(path, "w") as f:
        for k, v in cfg.items():
            f.write(f"{k:<12} {v}\n")

# ── Utility: run one case, return parsed stats dict ───────────────────────────
def run_case(cfg: dict, label: str, timeout: int = 120) -> dict:
    tmpdir = Path(tempfile.mkdtemp(prefix=f"fch_{label}_"))
    dat    = tmpdir / f"{label}.dat"
    make_dat(cfg, dat)

    t0 = time.time()
    try:
        proc = subprocess.run(
            [str(BINARY), str(dat)],
            capture_output=True, text=True, timeout=timeout
        )
        stdout = proc.stdout
    except subprocess.TimeoutExpired:
        shutil.rmtree(tmpdir, ignore_errors=True)
        return {"label": label, "status": "TIMEOUT", "wall_time": timeout}
    wall = time.time() - t0

    # cleanup output dir created by solver
    shutil.rmtree(tmpdir / f"{label}_Output", ignore_errors=True)
    shutil.rmtree(tmpdir, ignore_errors=True)

    m = SUMMARY_RE.search(stdout)
    if m:
        t_f, steps, vtk, avg_sor, max_sor, avg_res, avg_dt, status = m.groups()
        return dict(
            label    = label,
            status   = status,
            t_final  = float(t_f),
            steps    = int(steps),
            vtk      = int(vtk),
            avg_sor  = float(avg_sor),
            max_sor  = int(max_sor),
            avg_res  = float(avg_res),
            avg_dt   = float(avg_dt),
            wall_time= round(wall, 1),
        )

    # No SUMMARY → likely diverged before first output
    if "[DIVERGED]" in stdout:
        return {"label": label, "status": "DIVERGED", "wall_time": round(wall, 1)}
    return {"label": label, "status": "ERROR",
            "wall_time": round(wall, 1), "_out": stdout[-400:]}

# ── Utility: run case WITH VTK output kept in out_dir ────────────────────────
def run_case_vtk(cfg: dict, label: str, out_dir: Path, timeout: int = 600) -> dict:
    """Like run_case but keeps VTK output; out_dir must already exist."""
    dat = out_dir / f"{label}.dat"
    make_dat(cfg, dat)
    t0 = time.time()
    try:
        proc = subprocess.run(
            [str(BINARY), str(dat)],
            capture_output=True, text=True, timeout=timeout
        )
        stdout = proc.stdout
    except subprocess.TimeoutExpired:
        return {"label": label, "status": "TIMEOUT", "wall_time": timeout,
                "vtk_dir": out_dir / f"{label}_Output"}
    wall = time.time() - t0
    vtk_dir = out_dir / f"{label}_Output"
    m = SUMMARY_RE.search(stdout)
    if m:
        t_f, steps, vtk, avg_sor, max_sor, avg_res, avg_dt, status = m.groups()
        return dict(label=label, status=status, t_final=float(t_f),
                    steps=int(steps), vtk=int(vtk), avg_sor=float(avg_sor),
                    max_sor=int(max_sor), avg_res=float(avg_res),
                    avg_dt=float(avg_dt), wall_time=round(wall,1), vtk_dir=vtk_dir)
    if "[DIVERGED]" in stdout:
        return {"label": label, "status": "DIVERGED",
                "wall_time": round(wall,1), "vtk_dir": vtk_dir}
    return {"label": label, "status": "ERROR",
            "wall_time": round(wall,1), "vtk_dir": vtk_dir}

# ── Utility: render final-state velocity+streamlines via pvpython ─────────────
_PV_RENDER_SCRIPT = r"""
import sys, glob
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

vtk_dir, out_vel, out_stream = sys.argv[1], sys.argv[2], sys.argv[3]

vtk_files = sorted(glob.glob(vtk_dir + '/*.vtk'))
if not vtk_files:
    sys.exit(1)

# ── Load + go to last timestep ────────────────────────────────────────────────
reader = LegacyVTKReader(FileNames=vtk_files)
scene  = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.GoToLast()

# ── Calculator pipeline (same as visualize.py) ────────────────────────────────
c1 = Calculator(Input=reader); c1.AttributeType='Cell Data'
c1.ResultArrayName='u_comp';   c1.Function='velocity_X'
c2 = Calculator(Input=c1);    c2.AttributeType='Cell Data'
c2.ResultArrayName='v_comp';   c2.Function='velocity_Y'
c3 = Calculator(Input=c2);    c3.AttributeType='Cell Data'
c3.ResultArrayName='vel_mag';  c3.Function='mag(velocity)'

# ── Shared camera (2-D parallel, centered on unit square) ────────────────────
def setup_camera(v):
    v.CameraPosition           = [0.5, 0.5, 3.0]
    v.CameraFocalPoint         = [0.5, 0.5, 0.0]
    v.CameraViewUp             = [0.0, 1.0, 0.0]
    v.CameraParallelProjection = 1
    v.CameraParallelScale      = 0.55

def add_colorbar(lut, view, title):
    sb = GetScalarBar(lut, view)
    sb.Title           = title
    sb.ComponentTitle  = ''
    sb.TitleFontSize   = 13
    sb.LabelFontSize   = 11
    sb.WindowLocation  = 'Lower Right Corner'
    sb.Visibility      = 1

IMG = [1200, 1000]

# ── Image 1: velocity magnitude (Jet 0..1, dark bg) ─────────────────────────
v1 = CreateRenderView()
v1.ViewSize                = IMG
v1.Background              = [0.12, 0.12, 0.18]
v1.OrientationAxesVisibility = 0
d1 = Show(c3, v1)
d1.Representation = 'Surface'
ColorBy(d1, ('CELLS', 'vel_mag'))
lut1 = GetColorTransferFunction('vel_mag')
lut1.ApplyPreset('Jet', True)
lut1.RescaleTransferFunction(0.0, 1.0)
d1.SetScalarBarVisibility(v1, True)
add_colorbar(lut1, v1, 'Velocity Magnitude |u|  [m/s]')
setup_camera(v1)
Render(v1)
SaveScreenshot(out_vel, v1, ImageResolution=IMG)

# ── Image 2: streamlines over vel-mag background (dark blue bg, Jet) ─────────
v2 = CreateRenderView()
v2.ViewSize                = IMG
v2.Background              = [0.05, 0.05, 0.10]
v2.OrientationAxesVisibility = 0
d_bg = Show(c3, v2)
d_bg.Representation = 'Surface'
ColorBy(d_bg, ('CELLS', 'vel_mag'))
lut_bg = GetColorTransferFunction('vel_mag')
lut_bg.ApplyPreset('Jet', True)
lut_bg.RescaleTransferFunction(0.0, 1.0)
d_bg.Opacity = 0.5

stream = StreamTracer(Input=reader, SeedType='Line')
stream.Vectors                 = ['POINTS', 'velocity']
stream.MaximumStreamlineLength = 3.0
stream.IntegrationDirection    = 'BOTH'
stream.SeedType.Point1         = [0.05, 0.5, 0.0]
stream.SeedType.Point2         = [0.95, 0.5, 0.0]
stream.SeedType.Resolution     = 40
d_st = Show(stream, v2)
d_st.Representation = 'Surface'
d_st.LineWidth      = 2.0
ColorBy(d_st, ('POINTS', 'velocity'))
lut_st = GetColorTransferFunction('velocity')
lut_st.ApplyPreset('Jet', True)
lut_st.RescaleTransferFunction(0.0, 1.0)
d_st.SetScalarBarVisibility(v2, True)
add_colorbar(lut_st, v2, 'Velocity |u|  [m/s]')
setup_camera(v2)
Render(v2)
SaveScreenshot(out_stream, v2, ImageResolution=IMG)

print('PVRENDER_OK')
"""

def pvrender_case(vtk_dir: Path, out_vel: Path, out_stream: Path) -> bool:
    """Render velocity magnitude and streamlines for the last timestep."""
    if PVPYTHON is None:
        return False
    tmp = Path(tempfile.mktemp(suffix="_pv.py"))
    tmp.write_text(_PV_RENDER_SCRIPT)
    try:
        proc = subprocess.run(
            [str(PVPYTHON), str(tmp),
             str(vtk_dir), str(out_vel), str(out_stream)],
            capture_output=True, text=True, timeout=120
        )
        return "PVRENDER_OK" in proc.stdout
    except Exception:
        return False
    finally:
        tmp.unlink(missing_ok=True)

# ── Utility: print a bordered ASCII table ────────────────────────────────────
def print_table(headers: list, rows: list) -> None:
    cols   = len(headers)
    widths = [len(h) for h in headers]
    for row in rows:
        for c, v in enumerate(row):
            widths[c] = max(widths[c], len(str(v)))
    sep = "+" + "+".join("-" * (w + 2) for w in widths) + "+"
    fmt = "|" + "|".join(f" {{:<{w}}} " for w in widths) + "|"
    print(sep)
    print(fmt.format(*headers))
    print(sep)
    for row in rows:
        print(fmt.format(*[str(v) for v in row]))
    print(sep)

def section(title: str) -> None:
    print("\n" + "=" * 66)
    print(f"  {title}")
    print("=" * 66)

# ── Task 4: Base LDC simulation ───────────────────────────────────────────────
def task4():
    section("Task 4 — Base Lid-Driven Cavity Simulation")
    print("Grid: 50×50  nu=0.01  Re=100  adaptive dt (tau=0.5)  t_end=50.0")
    print("VTK output every dt_value=0.5 s  →  100 snapshots\n")

    cfg = {**BASE_CFG, "t_end": 50.0, "dt_value": 0.5}
    r   = run_case(cfg, "ldc_base", timeout=300)
    print(f"  done  ({r['status']})\n", flush=True)

    rows = [[
        "50×50",
        "0.01",
        "100",
        f"{r.get('avg_dt', '-'):.2e}" if "avg_dt" in r else "-",
        r.get("steps", "-"),
        r.get("vtk", "-"),
        f"{r.get('avg_sor', '-'):.1f}" if "avg_sor" in r else "-",
        f"{r.get('avg_res', '-'):.3f}" if "avg_res" in r else "-",
        r["status"],
    ]]
    print_table(
        ["grid", "nu", "Re", "avg dt", "steps", "VTK files", "avg SOR iter", "avg residual", "status"],
        rows,
    )
    print()


# ── Task 5: SOR omega and itermax study ───────────────────────────────────────
def task5():
    section("Task 5 — SOR Solver Behavior")
    print("Grid: 50×50  nu=0.01  Re≈100  t_end=50.0  (adaptive dt, tau=0.5)\n")

    # ── 5a: vary omega ──────────────────────────────────────────────────────
    # Deliberately use an unrealistically small residual target so every run is
    # iteration-budget-limited. We then compare omega values by the residual
    # reached after the same number of SOR iterations, not by early convergence.
    print(f"── 5a: Effect of relaxation factor omega "
          f"(itermax={SOR_OMEGA_ITERMAX}, eps={SOR_OMEGA_EPS:.0e}) ──\n")

    omegas = [0.5, 1.0, 1.3, 1.5, 1.7, 1.8, 1.9, 1.95, 1.99]
    rows   = []
    for omg in omegas:
        cfg = {**BASE_CFG,
               "omg": omg,
               "eps": SOR_OMEGA_EPS,
               "itermax": SOR_OMEGA_ITERMAX,
               "t_end": 2.0}
        r   = run_case(cfg, f"omega_{omg}", timeout=300)
        hits = "YES" if r.get("max_sor", 0) >= SOR_OMEGA_ITERMAX else "no"
        rows.append([
            omg,
            f"{r.get('avg_sor', '-'):.0f}" if "avg_sor" in r else "-",
            hits,
            f"{r.get('avg_res', '-'):.3e}" if "avg_res" in r else "-",
            f"{r['wall_time']}s",
            r["status"],
        ])
        print(f"  omega={omg:.2f}  done  ({r['status']})", flush=True)

    print()
    print_table(
        ["omega", "avg SOR iter", "hits itermax", "avg residual", "wall time", "status"],
        rows,
    )

    # Since eps is intentionally unreachable, 'best' means the omega that gives
    # the smallest average residual after the fixed iteration budget.
    valid_rows5a = [r for r in rows if r[3] != "-"]
    if valid_rows5a:
        best_row5a = min(valid_rows5a, key=lambda r: float(r[3]))
        best_omg = float(best_row5a[0])
        print(f"\n  → Best omega ≈ {best_omg}  "
              f"(lowest avg residual after {SOR_OMEGA_ITERMAX} SOR iterations: "
              f"{float(best_row5a[3]):.3e})")
    else:
        best_omg = 1.7
        print("\n  → No valid residuals parsed in 5a; using omega=1.7 for 5b.")

    missed_cap = [r for r in rows if r[2] == "no"]
    if missed_cap:
        print("  WARNING: at least one run did not hit itermax. If the purpose is a")
        print("  strictly fixed-budget comparison, decrease SOR_OMEGA_EPS further or")
        print("  inspect whether the C++ solver is stopping for another reason.")

    print()

    # ── Plot 5a helper ─────────────────────────────────────────────────────────
    def _plot_5a(eps_ref, out_path):
        omg_vals = [float(r[0]) for r in rows if r[3] != "-"]
        res_vals = [float(r[3]) for r in rows if r[3] != "-"]
        colors   = ["tab:green" if o == best_omg else "tab:blue" for o in omg_vals]
        fig, ax  = plt.subplots(figsize=(10, 6))
        bars = ax.bar([str(o) for o in omg_vals], res_vals,
                      color=colors, edgecolor="black", linewidth=0.6)
        y_top = max(res_vals) * 1.35
        for bar, v in zip(bars, res_vals):
            ax.text(bar.get_x() + bar.get_width() / 2, v + y_top * 0.02,
                    f"{v:.3f}", ha="center", va="bottom", fontsize=8)
        ax.axhline(eps_ref, color="purple", linewidth=1.2, linestyle=":",
                   label=f"ε = {eps_ref:.0e}  (unreachable target)")
        ax.set_ylim(0.010, y_top)
        ax.set_xlabel("Relaxation factor ω")
        ax.set_ylabel("Avg SOR residual after fixed iteration budget")
        ax.set_title(f"Task 5a — SOR Residual vs. Relaxation Factor ω\n"
                     f"(itermax={SOR_OMEGA_ITERMAX}, eps={eps_ref:.0e}, "
                     f"50×50 grid, Re=100, solver={_ACTIVE_SOLVER})")
        green_patch = mpatches.Patch(color="tab:green", label=f"Lowest residual: ω = {best_omg}")
        purple_patch = mpatches.Patch(color="purple", label=f"ε = {eps_ref:.0e}")
        ax.legend(handles=[green_patch, purple_patch])
        fig.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)

    _plot_5a(1e-15, PLOTS_DIR / "task5a_omega_residual_unreachable_eps.png")
    _plot_5a(1e-7,  PLOTS_DIR / "task5a_omega_residual.png")
    print(f"\n  → Plot saved: study_plots/task5a_omega_residual_unreachable_eps.png  (ε=1e-15)")
    print(f"  → Plot saved: study_plots/task5a_omega_residual.png  (ε=1e-7)")

    # ── 5b: vary itermax ───────────────────────────────────────────────────
    print(f"\n── 5b: Effect of itermax (omega={best_omg} — optimal from 5a) ──\n")
    itermaxs = [5, 10, 20, 50, 100, 200, 500]
    rows2    = []
    for im in itermaxs:
        cfg = {**BASE_CFG, "omg": best_omg, "itermax": im, "t_end": 50.0}
        r   = run_case(cfg, f"itermax_{im}", timeout=300)
        rows2.append([
            im,
            f"{r.get('avg_sor', '-'):.1f}" if "avg_sor" in r else "-",
            r.get("max_sor", "-"),
            f"{r.get('avg_res', '-'):.4f}" if "avg_res" in r else "-",
            r["status"],
        ])
        print(f"  itermax={im:>3}  done  ({r['status']})", flush=True)

    print()
    print_table(
        ["itermax", "avg SOR iter", "max SOR iter", "avg residual", "status"],
        rows2,
    )
    print("\n  Influence of itermax (with optimal ω, t_end=50):")

    # ── Plot 5b: two side-by-side subplots ────────────────────────────────────
    im_vals   = [r[0] for r in rows2]
    iter_vals = [float(r[1]) if r[1] != "-" else 0.0 for r in rows2]
    res_vals  = [float(r[3]) if r[3] != "-" else 0.0 for r in rows2]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    # Left: avg SOR iterations used vs itermax (efficiency)
    bars5b_1 = ax1.bar([str(v) for v in im_vals], iter_vals, color="tab:orange",
                       edgecolor="black", linewidth=0.6)
    ax1.set_ylim(0, max(iter_vals) * 1.3)
    for bar in bars5b_1:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(iter_vals) * 0.02,
                 f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=8)
    ax1.set_xlabel("itermax")
    ax1.set_ylabel("Avg SOR iterations used")
    ax1.set_title("Iterations used vs. itermax\n(efficiency)")

    # Right: avg SOR residual vs itermax (accuracy)
    ax2.bar([str(v) for v in im_vals], res_vals, color="tab:red",
            edgecolor="black", linewidth=0.6)
    ax2.axhline(0.001, color="green", linewidth=1.5, linestyle="--",
                label="ε = 0.001  (convergence target)")
    # Value labels above each bar — log-aware offset, set ylim BEFORE annotating
    ax2.set_yscale("log")
    ax2.set_ylim(3e-4, max(res_vals) * 12)  # explicit headroom so labels never clip
    for x, v in zip([str(v) for v in im_vals], res_vals):
        ax2.text(x, v * 2.0, f"{v:.4f}", ha="center", va="bottom", fontsize=8)
    ax2.set_xlabel("itermax")
    ax2.set_ylabel("Avg SOR residual  [log scale]")
    ax2.set_title("Residual vs. itermax\n(accuracy — log scale)")
    ax2.legend()

    fig.suptitle(f"Task 5b — Effect of itermax  (ω={best_omg}, 50×50, Re=100, solver={_ACTIVE_SOLVER})", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "task5b_itermax.png")
    plt.close(fig)
    print(f"  → Plot saved: study_plots/task5b_itermax.png")

# ── Task 6: Fixed time step stability ─────────────────────────────────────────
def task6():
    section("Task 6 — Fixed Time Step Stability")
    dx = 1.0 / 50   # dx = 0.02 for 50×50
    print(f"Grid: 50×50  nu=0.01  dx={dx:.4f}  t_end=5.0")
    print(f"Adaptive dt DISABLED (tau=-1 → fixed dt from input file)")
    print(f"CFL stability limit: dt < dx / u_max ≈ {dx:.4f} / 1.0 = {dx:.4f}\n")

    nu   = 0.01
    # Viscous stability limit (Eq.12): dt_visc = dx²·dy²/(2·nu·(dx²+dy²))
    dt_visc = (dx**2 * dx**2) / (2 * nu * (dx**2 + dx**2))
    dts  = [0.001, 0.005, 0.008, 0.010, 0.012, 0.015, 0.020, 0.025, 0.050]
    rows = []
    for dt in dts:
        cfg = {**BASE_CFG, "dt": dt, "tau": -1}   # tau=-1 disables adaptive stepping
        r   = run_case(cfg, f"dt_{dt}")
        cfl  = dt / dx        # CFL with u_max ~ 1
        visc = dt / dt_visc   # fraction of viscous limit
        rows.append([
            dt,
            f"{cfl:.3f}",
            "✓" if cfl < 1.0 else "✗",
            f"{visc:.2f}",
            "✓" if visc < 1.0 else "✗",
            r["status"],
            f"{r['wall_time']}s",
        ])
        print(f"  dt={dt:.3f}  CFL≈{cfl:.2f}  visc={visc:.2f}  {r['status']}", flush=True)

    print()
    print_table(
        ["dt", "CFL=dt/dx", "CFL<1?", "dt/dt_visc", "visc<1?", "status", "time"],
        rows,
    )

    # ── Plot 6 ─────────────────────────────────────────────────────────────────
    dt_vals  = [r[0] for r in rows]
    statuses = [r[5] for r in rows]
    cfl_vals = [float(r[1]) for r in rows]
    visc_vals= [float(r[3]) for r in rows]
    bar_colors = ["tab:green" if s == "OK" else "tab:red" for s in statuses]

    # ── Combined figure: outcome (left) + stability fractions (right) ───────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    green_p = mpatches.Patch(color="tab:green", label="OK (stable)")
    red_p   = mpatches.Patch(color="tab:red",   label="DIVERGED")

    # Left: stability outcome
    ax1.bar([str(d) for d in dt_vals], [1]*len(dt_vals),
            color=bar_colors, edgecolor="black", linewidth=0.6)
    ax1.set_xlabel("Fixed dt")
    ax1.set_ylabel("Stable / Diverged")
    ax1.set_yticks([])
    ax1.set_title("Simulation Outcome per dt")
    ax1.tick_params(axis="x", rotation=45)
    ax1.legend(handles=[green_p, red_p])

    # Right: CFL and viscous stability fractions
    x = range(len(dt_vals))
    w = 0.38
    b_cfl  = ax2.bar([i - w/2 for i in x], cfl_vals,  width=w, label="CFL = dt/dx",  color="tab:blue",   edgecolor="black", linewidth=0.5)
    b_visc = ax2.bar([i + w/2 for i in x], visc_vals, width=w, label="dt/dt_visc",   color="tab:orange", edgecolor="black", linewidth=0.5)
    ax2.axhline(1.0, color="red", linewidth=1.4, linestyle="--", label="limit = 1")
    ymax = max(max(cfl_vals), max(visc_vals)) * 1.22
    ax2.set_ylim(0, ymax)
    offset    = ymax * 0.018
    clearance = ymax * 0.07   # minimum gap above the red limit line at y=1
    def _label_y(h):
        """Place label above bar, but never inside the red limit-line band."""
        y = h + offset
        if h > 0.5 and y < 1.0 + clearance:
            y = 1.0 + clearance
        return y
    for bar in b_cfl:
        ax2.text(bar.get_x() + bar.get_width()/2, _label_y(bar.get_height()),
                 f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=8, color="tab:blue")
    for bar in b_visc:
        ax2.text(bar.get_x() + bar.get_width()/2, _label_y(bar.get_height()),
                 f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=8, color="tab:orange")
    ax2.set_xticks(list(x))
    ax2.set_xticklabels([str(d) for d in dt_vals], rotation=45)
    ax2.set_xlabel("Fixed dt")
    ax2.set_ylabel("Fraction of stability limit")
    ax2.set_title("CFL and Viscous Stability Fractions")
    ax2.legend(fontsize=9)

    fig.suptitle(f"Task 6 — Fixed Time Step Stability  (50×50, nu=0.01, solver={_ACTIVE_SOLVER})", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "task6_dt_stability.png")
    plt.close(fig)
    print(f"  → Plot saved: study_plots/task6_dt_stability.png")

# ── Task 7: Grid refinement study ─────────────────────────────────────────────
def task7():
    section("Task 7 — Grid Refinement (fixed dt=0.05, nu=0.001, Re=1000)")
    print("Fixed dt=0.05  tau=-1 (disabled)  t_end=5.0  nu=0.001")
    print("CFL ≈ dt / dx   (upper bound with u_max ~ 1)   →   stability requires CFL < 1  (dx > dt)\n")

    grids = [16, 32, 64, 128, 256]
    rows  = []
    for n in grids:
        dx  = 1.0 / n
        cfg = {**BASE_CFG, "imax": n, "jmax": n, "dt": 0.05, "tau": -1, "nu": 0.001}
        r   = run_case(cfg, f"grid_{n}", timeout=60)
        cfl = 0.05 / dx
        rows.append([
            f"{n}×{n}",
            f"{dx:.5f}",
            f"{cfl:.2f}",
            "< 1 ✓" if cfl < 1.0 else "> 1 ✗",
            r["status"],
            f"{r.get('avg_sor', '-'):.1f}" if "avg_sor" in r else "-",
            f"{r['wall_time']}s",
        ])
        print(f"  {n:>3}×{n:<3}  dx={dx:.5f}  CFL≈{cfl:.2f}  {r['status']}", flush=True)

    print()
    print_table(
        ["grid", "dx", "CFL≈dt/dx", "criterion", "status", "avg SOR", "wall time"],
        rows,
    )
    print()

    # ── 7b: Visualize stable grid cases ───────────────────────────────────────
    print("\n── 7b: Flow visualization for stable grids ──\n")
    VIS7_DIR = SCRIPT_DIR / "LidDrivenCavity_Output" / "task7_visuals"
    VIS7_DIR.mkdir(parents=True, exist_ok=True)

    vel7_imgs    = {}   # n → Path
    stream7_imgs = {}

    for n in grids:
        dx  = 1.0 / n
        cfl = 0.05 / dx
        label7 = f"grid{n}"
        if cfl >= 1.0:
            label7 += "_marginal"   # CFL>1 but may still run for short t_end
        case_dir = VIS7_DIR / label7
        case_dir.mkdir(exist_ok=True)
        cfg_vis  = {**BASE_CFG, "imax": n, "jmax": n,
                    "dt": 0.05, "tau": -1, "nu": 0.001, "dt_value": 0.5}
        marginal = cfl >= 1.0
        tag = f"CFL={cfl:.2f}>1 — marginal, OK for t_end=5 s only" if marginal else f"CFL={cfl:.2f}<1"
        print(f"  Running {n}×{n} ({tag}) ...", flush=True)
        r = run_case_vtk(cfg_vis, label7, case_dir, timeout=120)
        print(f"    solver: {r['status']}", flush=True)
        if r["status"] != "OK":
            print(f"    → skipping render (diverged)")
            continue
        if r.get("vtk_dir") and r["vtk_dir"].exists():
            out_vel7    = VIS7_DIR / f"{label7}_velocity.png"
            out_stream7 = VIS7_DIR / f"{label7}_streamlines.png"
            ok = pvrender_case(r["vtk_dir"], out_vel7, out_stream7)
            if ok:
                vel7_imgs[n]    = out_vel7
                stream7_imgs[n] = out_stream7
                print(f"    → images saved: {label7}_velocity.png  {label7}_streamlines.png")
            else:
                print(f"    → pvpython rendering failed")

    if vel7_imgs:
        keys     = sorted(vel7_imgs.keys())
        img_rows = [[vel7_imgs.get(n),    stream7_imgs.get(n)] for n in keys]
        img_rows = list(map(list, zip(*img_rows)))  # transpose: rows=type, cols=grid
        col_titles = [
            f"{n}×{n}  (dx={1/n:.4f})" + ("  [marginal]" if 0.05/(1/n) >= 1.0 else "")
            for n in keys
        ]
        out_cmp7 = PLOTS_DIR / "task7_grid_comparison.png"
        sz = make_comparison_grid(
            img_rows, col_titles,
            ["Velocity magnitude", "Streamlines"],
            "Task 7 — Grid Refinement  (fixed dt=0.05, nu=0.001, Re=1000)",
            out_cmp7,
        )
        print(f"\n  → Comparison grid saved: study_plots/task7_grid_comparison.png  {sz}")

    # ── Plot 7 ─────────────────────────────────────────────────────────────────
    grid_labels = [r[0] for r in rows]
    cfl_vals7   = [float(r[2]) for r in rows]
    statuses7   = [r[4] for r in rows]
    bar_colors7 = ["tab:green" if s == "OK" else "tab:red" for s in statuses7]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    bars7_1 = ax1.bar(grid_labels, cfl_vals7, color=bar_colors7, edgecolor="black", linewidth=0.6)
    ax1.axhline(1.0, color="red", linewidth=1.4, linestyle="--", label="CFL limit = 1")
    ax1.set_ylim(0, max(cfl_vals7) * 1.2)
    for bar in bars7_1:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(cfl_vals7) * 0.02,
                 f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=8)
    ax1.set_xlabel("Grid resolution")
    ax1.set_ylabel("CFL ≈ dt / dx")
    ax1.set_title("CFL Number per Grid  (fixed dt=0.05, itermax=100)")
    green_p = mpatches.Patch(color="tab:green", label="OK (stable)")
    red_p   = mpatches.Patch(color="tab:red",   label="DIVERGED")
    ax1.legend(handles=[green_p, red_p, mpatches.Patch(color="none")])
    ax1.legend(handles=[green_p, red_p,
               mpatches.Patch(color="none", label="— CFL limit")])

    sor_vals7 = [float(r[5]) if r[5] != "-" else 0 for r in rows]
    bars7_2 = ax2.bar(grid_labels, sor_vals7, color=bar_colors7, edgecolor="black", linewidth=0.6)
    for bar in bars7_2:
        if bar.get_height() > 0:
            ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                     f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=8)
    ax2.set_xlabel("Grid resolution")
    ax2.set_ylabel("Avg SOR iterations")
    ax2.set_title("Avg SOR Iterations per Grid  (itermax=100)")
    ax2.legend(handles=[green_p, red_p])

    fig.suptitle(f"Task 7 — Grid Refinement Study  (fixed dt=0.05, nu=0.001, Re=1000, itermax=100, solver={_ACTIVE_SOLVER})", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "task7_grid_refinement.png")
    plt.close(fig)
    print(f"  → Plot saved: study_plots/task7_grid_refinement.png")

# ── Task 8: Viscosity / Reynolds number study ─────────────────────────────────
def task8():
    section("Task 8 — Kinematic Viscosity Study (adaptive dt, t_end=50.0)")
    print("Grid: 50×50  adaptive dt (tau=0.5)  t_end=50.0\n")

    nus  = [0.01, 0.002, 0.0005, 0.0001]
    rows = []
    for nu in nus:
        Re  = 1.0 / nu   # U=L=1
        cfg = {**BASE_CFG, "nu": nu, "tau": 0.5, "t_end": 50.0}
        r   = run_case(cfg, f"nu_{nu}", timeout=600)
        avg_dt_str = f"{r['avg_dt']:.2e}" if "avg_dt" in r else "-"
        rows.append([
            nu,
            int(Re),
            avg_dt_str,
            f"{r.get('avg_sor', '-'):.1f}" if "avg_sor" in r else "-",
            r.get("max_sor", "-"),
            "YES" if r.get("max_sor", 0) >= 100 else "no",
            r["status"],
        ])
        print(f"  nu={nu}  Re={int(Re):>5}  {r['status']}", flush=True)

    print()
    print_table(
        ["nu", "Re", "avg dt", "avg SOR iter", "max SOR", "hits itermax", "status"],
        rows,
    )
    print()
    # ── 8c: Visual output per Re (velocity + streamlines) ─────────────────────
    print("\n── 8c: Flow visualization per Re (velocity magnitude + streamlines) ──\n")
    VIS_DIR = SCRIPT_DIR / "LidDrivenCavity_Output" / "task8_visuals"
    VIS_DIR.mkdir(parents=True, exist_ok=True)

    vel_imgs    = {}   # Re → Path of velocity PNG
    stream_imgs = {}   # Re → Path of streamlines PNG

    for nu in nus:
        Re    = int(1.0 / nu)
        label = f"Re{Re}"
        case_dir = VIS_DIR / label
        case_dir.mkdir(exist_ok=True)

        cfg_vis = {**BASE_CFG, "nu": nu, "tau": 0.5, "t_end": 50.0, "dt_value": 0.5}
        print(f"  Running Re={Re:>5} (nu={nu}) for VTK output ...", flush=True)
        r = run_case_vtk(cfg_vis, label, case_dir, timeout=300)
        print(f"    solver: {r['status']}", flush=True)

        if r["status"] == "OK" and r.get("vtk_dir") and r["vtk_dir"].exists():
            out_vel    = VIS_DIR / f"{label}_velocity.png"
            out_stream = VIS_DIR / f"{label}_streamlines.png"
            ok = pvrender_case(r["vtk_dir"], out_vel, out_stream)
            if ok:
                vel_imgs[Re]    = out_vel
                stream_imgs[Re] = out_stream
                print(f"    → images saved: {label}_velocity.png  {label}_streamlines.png")
            else:
                print(f"    → pvpython rendering failed (pvpython={PVPYTHON})")
        else:
            print(f"    → skipping render ({r['status']})")

    # ── Assemble 2-row comparison grid (vel | streamlines) × 4 Re values ──────
    if vel_imgs:
        res_keys   = sorted(vel_imgs.keys())
        img_rows   = [[vel_imgs.get(Re),    stream_imgs.get(Re)] for Re in res_keys]
        img_rows   = list(map(list, zip(*img_rows)))  # transpose
        col_titles = [f"Re = {Re}\n(nu = {1/Re:.4f})" for Re in res_keys]
        out_cmp    = PLOTS_DIR / "task8_re_comparison.png"
        sz = make_comparison_grid(
            img_rows, col_titles,
            ["Velocity magnitude", "Streamlines"],
            "Task 8 — LDC Flow at Different Reynolds Numbers  (t = 10 s)",
            out_cmp,
        )
        print(f"\n  → Comparison grid saved: study_plots/task8_re_comparison.png  {sz}")

    # ── Plot 8 ─────────────────────────────────────────────────────────────────
    re_vals   = [int(1.0/nu) for nu in nus]
    avg_dts   = [r[2] for r in rows]
    dt_floats = []
    for s in avg_dts:
        try:    dt_floats.append(float(s))
        except: dt_floats.append(0.0)
    dt_visc_vals = [(1.0/50)**2 / (4*nu) for nu in nus]   # dx²/(4ν)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    x_pos = range(len(re_vals))
    w = 0.38
    bars1 = ax1.bar([i - w/2 for i in x_pos], dt_floats,    width=w,
                    color="tab:blue",   edgecolor="black", linewidth=0.6, label="avg adaptive dt")
    bars2 = ax1.bar([i + w/2 for i in x_pos], dt_visc_vals, width=w,
                    color="tab:orange", edgecolor="black", linewidth=0.6, label="dt_visc = dx²/(4ν)")
    # Annotate exact values above each bar
    for bar in bars1:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.0003,
                 f"{bar.get_height():.4f}", ha="center", va="bottom", fontsize=8, color="tab:blue")
    for bar in bars2:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.0003,
                 f"{bar.get_height():.4f}", ha="center", va="bottom", fontsize=8, color="tab:orange")
    ax1.set_xticks(list(x_pos))
    ax1.set_xticklabels([str(r) for r in re_vals])
    ax1.set_xlabel("Reynolds number Re")
    ax1.set_ylabel("Time step dt  [s]")
    ax1.set_title("Adaptive dt vs. Re\n(blue = used, orange = viscous limit)")
    ax1.legend()

    sor_vals8 = [float(r[3]) if r[3] != "-" else 0 for r in rows]
    bars3 = ax2.bar([str(r) for r in re_vals], sor_vals8, color="tab:purple",
                    edgecolor="black", linewidth=0.6)
    ax2.axhline(100, color="red", linewidth=1.4, linestyle="--", label="itermax = 100")
    # Annotate exact values
    for bar in bars3:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=9)
    ax2.set_xlabel("Reynolds number Re")
    ax2.set_ylabel("Avg SOR iterations")
    ax2.set_title("SOR Iterations vs. Re")
    ax2.legend()

    fig.suptitle(f"Task 8 — Viscosity / Reynolds Number Study  (50×50, adaptive dt, t_end=50, solver={_ACTIVE_SOLVER})", fontsize=11)
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "task8_viscosity_study.png")
    plt.close(fig)
    print(f"  → Plot saved: study_plots/task8_viscosity_study.png")

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    if not BINARY.exists():
        print(f"[ERROR] Binary not found: {BINARY}")
        print("  Run:  cd fluidchen-skeleton/build && cmake .. && make -j4")
        sys.exit(1)

    # Ensure plots output directory exists
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Decide which tasks to run
    args  = sys.argv[1:]
    tasks = [int(a) for a in args if a.isdigit()] or [4, 5, 6, 7, 8]

    print("=" * 66)
    print("  Worksheet 1 — Simulation Tasks 4–8")
    print(f"  Binary : {BINARY}")
    print(f"  Plots  : {PLOTS_DIR}")
    print(f"  Solver : {_ACTIVE_SOLVER}  (set in LidDrivenCavity.dat)")
    print("=" * 66)

    if 4 in tasks: task4()
    if 5 in tasks: task5()
    if 6 in tasks: task6()
    if 7 in tasks: task7()
    if 8 in tasks: task8()

    print("\n" + "=" * 66)
    print("  All selected studies complete.")
    print("  Results above answer Tasks 4–8 of Worksheet 1.")
    print(f"  Plots saved to: {PLOTS_DIR}")
    print("=" * 66 + "\n")

if __name__ == "__main__":
    main()
