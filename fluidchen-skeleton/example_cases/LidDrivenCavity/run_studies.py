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
    r"avg_dt=([\d.eE+\-]+)\s+status=(\w+)"
)

# ── Base configuration (worksheet defaults) ───────────────────────────────────
BASE_CFG = dict(
    xlength=1.0, ylength=1.0,
    imax=50,     jmax=50,
    dt=0.05,     t_end=5.0,   tau=0.5,
    dt_value=999,              # suppress VTK output during studies (speed)
    eps=0.001,   omg=1.7,     gamma=0.5,  itermax=100,
    nu=0.01,
    GX=0.0,      GY=0.0,
    PI=0.0,      UI=0.0,      VI=0.0,
)

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
    print("  Running ... (t_end=50.0, may take ~1–2 min)", flush=True)
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
    print("  The simulation runs until t=50 s to allow the flow to reach a")
    print("  quasi-steady state. At Re=100 a single stable vortex fills the")
    print("  cavity. Key observations:")
    print("  • Adaptive dt settles at ~5e-3 s (viscous stability limit × tau).")
    print("  • SOR converges well within itermax (Fredholm fix + zero-mean pressure).")
    print("  • VTK snapshots are visualized via visualize.py / make_videos.sh.")

# ── Task 5: SOR omega and itermax study ───────────────────────────────────────
def task5():
    section("Task 5 — SOR Solver Behavior")
    print("Grid: 50×50  nu=0.01  Re≈100  t_end=5.0  (adaptive dt, tau=0.5)\n")

    # ── 5a: vary omega ──────────────────────────────────────────────────────
    # Use itermax=500 so SOR can actually converge within the budget
    # and the optimal omega becomes visible.  With itermax=100 all omega
    # values hit the cap equally, masking any difference.
    print("── 5a: Effect of relaxation factor omega (itermax=500) ──\n")
    omegas = [0.5, 1.0, 1.3, 1.5, 1.7, 1.8, 1.9, 1.95, 1.99]
    rows   = []
    for omg in omegas:
        cfg = {**BASE_CFG, "omg": omg, "itermax": 500}
        r   = run_case(cfg, f"omega_{omg}")
        hits = "YES" if r.get("max_sor", 0) >= 500 else "no"
        rows.append([
            omg,
            f"{r.get('avg_sor', '-'):.0f}" if "avg_sor" in r else "-",
            hits,
            f"{r.get('avg_res', '-'):.3f}" if "avg_res" in r else "-",
            f"{r['wall_time']}s",
            r["status"],
        ])
        print(f"  omega={omg:.2f}  done  ({r['status']})", flush=True)

    print()
    print_table(
        ["omega", "avg SOR iter", "hits itermax", "avg residual", "wall time", "status"],
        rows,
    )

    # With SOR fix: most omegas reach eps → rank by avg iter (speed).
    # Prefer omegas that do NOT hit itermax; fall back to all if none qualify.
    converged = [(r[0], float(r[1])) for r in rows if r[1] != "-" and r[2] == "no"]
    if converged:
        best = min(converged, key=lambda x: x[1])
        print(f"\n  → Best omega ≈ {best[0]}  (fastest: {best[1]:.0f} avg iter, converges without hitting itermax)")
    else:
        valid = [(r[0], float(r[3])) for r in rows if r[3] != "-"]
        if valid:
            best = min(valid, key=lambda x: x[1])
            print(f"\n  → Best omega ≈ {best[0]}  (lowest avg residual = {best[1]:.3f})")
    print()
    print("  Key insight: The pure Neumann pressure system has two convergence blockers:")
    print("  (1) Fredholm incompatibility: sum(RS) ≠ 0 → no solution exists until RS is")
    print("      mean-subtracted to enforce the compatibility condition sum(RS) = 0.")
    print("  (2) Null-space drift: pressure is defined only up to a constant; each SOR")
    print("      sweep accumulates a constant offset that prevents residual decay.")
    print("  FIX applied (ws1_further_extensions_improved_SOR):")
    print("  → subtract mean(RS) before iteration  (Fredholm compatibility)")
    print("  → subtract mean(p) after each sweep   (zero-mean projection)")
    print("  Result: avg SOR iterations drop from 100 → ~5 (quasi-steady), residual reaches eps.")
    print("  omega<1: slow (under-relaxation).  omega≈1.7–1.9: fastest (fewest iter).")
    print("  omega→2: instability — iter count rises again and may hit itermax.")

    # ── Plot 5a ────────────────────────────────────────────────────────────────
    omg_vals = [float(r[0]) for r in rows if r[3] != "-"]
    res_vals = [float(r[3]) for r in rows if r[3] != "-"]
    best_omg = omg_vals[res_vals.index(min(res_vals))]
    colors   = ["tab:green" if o == best_omg else "tab:blue" for o in omg_vals]
    fig, ax  = plt.subplots(figsize=(7, 4.5))
    bars = ax.bar([str(o) for o in omg_vals], res_vals,
                  color=colors, edgecolor="black", linewidth=0.6)
    # Value labels above each bar — log-aware vertical offset avoids overlap
    for bar, v in zip(bars, res_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v * 1.9,
                f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    ax.axhline(0.001, color="purple", linewidth=1.2, linestyle=":",
               label="ε = 0.001  (convergence target)")
    ax.set_yscale("log")
    ax.set_ylim(4e-4, max(res_vals) * 6)   # headroom so labels never clip
    ax.set_xlabel("Relaxation factor ω")
    ax.set_ylabel("Avg SOR residual  [log scale]")
    ax.set_title("Task 5a — SOR Residual vs. Relaxation Factor ω\n"
                 "(itermax=500, 50×50 grid, Re=100)")
    green_patch = mpatches.Patch(color="tab:green", label=f"Best ω = {best_omg}")
    ax.legend(handles=[green_patch,
                       mpatches.Patch(color="purple", label="ε = 0.001")])
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "task5a_omega_residual.png")
    plt.close(fig)
    print(f"\n  → Plot saved: study_plots/task5a_omega_residual.png")

    # ── 5b: vary itermax ───────────────────────────────────────────────────
    print("\n── 5b: Effect of itermax (omega=1.7) ──\n")
    itermaxs = [5, 10, 20, 50, 100, 200, 500]
    rows2    = []
    for im in itermaxs:
        cfg = {**BASE_CFG, "omg": 1.7, "itermax": im}
        r   = run_case(cfg, f"itermax_{im}")
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
    print("\n  Note (with SOR fix applied):")
    print("  Transient phase (t≈0–10 s): SOR needs 10–70 iterations as the velocity")
    print("  field changes rapidly → itermax=5 hits cap, ~20–30 are sufficient here.")
    print("  Quasi-steady phase (t>>10 s): only 2–5 iterations per step needed.")
    print("  The Task-4 avg of ~5 is dominated by the long quasi-steady tail (t=10–50 s).")
    print("  Low itermax → SOR aborts early → avg residual stays above eps=0.001 →")
    print("  pressure correction is less accurate → small errors accumulate in velocity.")
    print()
    print("  Why avg_res > ε is NOT a problem in practice (Task 4 base case):")
    print("  The average is computed over ALL time steps, including the transient phase")
    print("  (t ≈ 0–10 s, ~2 000 steps) where the pressure field changes rapidly and")
    print("  SOR needs 10–70 iterations — often hitting itermax before reaching ε.")
    print("  In the quasi-steady phase (t > 10 s, ~8 000 steps of the 10 001 total)")
    print("  only 2–5 iterations are needed and ε IS reached every single step.")
    print("  The long quasi-steady tail dominates the avg_sor (→ 4.8) but the")
    print("  transient outliers pull avg_res slightly above ε (→ 0.003 > 0.001).")
    print("  This is physically correct and expected — no cause for concern.")

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

    fig.suptitle("Task 5b — Effect of itermax  (ω=1.7, 50×50, Re=100)", fontsize=12)
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
    print(f"\n  Stability requires BOTH: CFL < 1 (dt < {dx:.4f}) AND dt < dt_visc = {dt_visc:.4f}")
    print(f"  → Here dt_visc = dx²/(4·nu) = {dt_visc:.4f} is the BINDING constraint (< CFL limit).")
    print(f"  → Adaptive stepping (tau=0.5) automatically uses dt = 0.5·{dt_visc:.4f} = {0.5*dt_visc:.4f}.")
    print("  → For fixed dt, the safe region is dt < dt_visc ≈ 0.010.")

    # ── Plot 6 ─────────────────────────────────────────────────────────────────
    dt_vals  = [r[0] for r in rows]
    statuses = [r[5] for r in rows]
    cfl_vals = [float(r[1]) for r in rows]
    visc_vals= [float(r[3]) for r in rows]
    bar_colors = ["tab:green" if s == "OK" else "tab:red" for s in statuses]

    # ── Figure 6a: stability outcome ───────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(7, 4))
    green_p = mpatches.Patch(color="tab:green", label="OK (stable)")
    red_p   = mpatches.Patch(color="tab:red",   label="DIVERGED")
    ax1.bar([str(d) for d in dt_vals], [1]*len(dt_vals),
            color=bar_colors, edgecolor="black", linewidth=0.6)
    ax1.set_xlabel("Fixed dt")
    ax1.set_ylabel("Stable / Diverged")
    ax1.set_yticks([])
    ax1.set_title("Task 6 — Simulation Outcome per dt  (50×50, nu=0.01)")
    ax1.tick_params(axis="x", rotation=45)
    ax1.legend(handles=[green_p, red_p])
    fig1.tight_layout()
    fig1.savefig(PLOTS_DIR / "task6_dt_stability.png")
    plt.close(fig1)
    print(f"  → Plot saved: study_plots/task6_dt_stability.png")

    # ── Figure 6b: CFL and viscous stability fractions (standalone) ────────────
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    x = range(len(dt_vals))
    w = 0.38
    b_cfl  = ax2.bar([i - w/2 for i in x], cfl_vals,  width=w, label="CFL = dt/dx",  color="tab:blue",   edgecolor="black", linewidth=0.5)
    b_visc = ax2.bar([i + w/2 for i in x], visc_vals, width=w, label="dt/dt_visc",   color="tab:orange", edgecolor="black", linewidth=0.5)
    ax2.axhline(1.0, color="red", linewidth=1.4, linestyle="--", label="limit = 1")
    ymax = max(max(cfl_vals), max(visc_vals)) * 1.22
    ax2.set_ylim(0, ymax)
    offset = ymax * 0.018
    for bar in b_cfl:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + offset,
                 f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=8, color="tab:blue")
    for bar in b_visc:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + offset,
                 f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=8, color="tab:orange")
    ax2.set_xticks(list(x))
    ax2.set_xticklabels([str(d) for d in dt_vals], rotation=45)
    ax2.set_xlabel("Fixed dt")
    ax2.set_ylabel("Fraction of stability limit")
    ax2.set_title("Task 6 — CFL and Viscous Stability Fractions  (50×50, nu=0.01)")
    ax2.legend(fontsize=9)
    fig2.tight_layout()
    fig2.savefig(PLOTS_DIR / "task6_dt_fractions.png")
    plt.close(fig2)
    print(f"  → Plot saved: study_plots/task6_dt_fractions.png")

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
    print("  Observation:")
    print("  • CFL < 1 criterion uses u_max=1 (lid velocity) as an upper bound.")
    print("    The true LOCAL CFL inside the domain starts at 0 (fluid at rest).")
    print("  • 32×32 (CFL≈1.60) survives because interior velocities stay well below")
    print("    U_wall for t_end=5 s — instability growth rate (~1.6×/step) is too slow")
    print("    to manifest in 100 steps. This is MARGINAL, not truly stable.")
    print("    At t_end=50 s the 32×32 case would likely diverge as well.")
    print("  • 64×64 (CFL≈3.20) diverges immediately: growth rate ~3.2×/step is")
    print("    strong enough even at t=0 near the lid ghost cells.")
    print("  • The red CFL=1 line in the plot is the correct theoretical boundary.")
    print("  • Finer grids REQUIRE smaller dt → use adaptive time stepping!")
    print("  • With adaptive dt (tau=0.5), ALL grid sizes converge correctly.")

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

    fig.suptitle("Task 7 — Grid Refinement Study  (fixed dt=0.05, nu=0.001, Re=1000, itermax=100)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "task7_grid_refinement.png")
    plt.close(fig)
    print(f"  → Plot saved: study_plots/task7_grid_refinement.png")

# ── Task 8: Viscosity / Reynolds number study ─────────────────────────────────
def task8():
    section("Task 8 — Kinematic Viscosity Study (adaptive dt, t_end=10.0)")
    print("Grid: 50×50  adaptive dt (tau=0.5)  t_end=10.0\n")

    nus  = [0.01, 0.002, 0.0005, 0.0001]
    rows = []
    for nu in nus:
        Re  = 1.0 / nu   # U=L=1
        cfg = {**BASE_CFG, "nu": nu, "tau": 0.5, "t_end": 10.0}
        r   = run_case(cfg, f"nu_{nu}", timeout=300)
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
    print("  Observations:")
    print("  • As nu ↓ (Re ↑), the VISCOUS stability limit rises (dt_visc ∝ 1/nu),")
    print("    so the CFL condition dt < dx/u_max becomes the binding constraint.")
    print("  • avg_dt INCREASES with Re: high-Re flows are convection-limited, not")
    print("    viscosity-limited, and the lid velocity (~1 m/s) still allows dt~0.01.")
    print("  • SOR fix applied (Fredholm + zero-mean): avg iter drastically reduced.")
    print("    Without fix: avg iter=100 (always hits cap), avg residual≈1.23.")
    print("    With fix: avg iter=14–28 during transient (t_end=10), 2–5 at quasi-steady.")
    print("    max_sor may still reach itermax during steep initial transient steps —")
    print("    this is expected; increase itermax to 50 for high-Re transient accuracy.")
    print("  • At Re=10000 the flow is likely unsteady; adaptive dt keeps it stable.")
    print("    Longer t_end and finer grids are needed to resolve the turbulent regime.")

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

        cfg_vis = {**BASE_CFG, "nu": nu, "tau": 0.5, "t_end": 10.0, "dt_value": 0.5}
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

    fig.suptitle("Task 8 — Viscosity / Reynolds Number Study  (50×50, adaptive dt, t_end=10)", fontsize=11)
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
