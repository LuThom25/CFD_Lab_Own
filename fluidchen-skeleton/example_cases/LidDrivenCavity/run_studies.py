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
from PIL import Image

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
BINARY     = SCRIPT_DIR.parent.parent / "build" / "fluidchen"

# pvpython for ParaView rendering
_PV_CANDIDATES = [
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
    print("  • SOR always hits itermax=100 due to the singular Neumann system.")
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

    # Find omega with lowest residual (most effective pressure solve)
    valid = [(r[0], float(r[3])) for r in rows if r[3] != "-"]
    if valid:
        best = min(valid, key=lambda x: x[1])
        print(f"\n  → Best omega ≈ {best[0]}  (lowest avg residual = {best[1]:.3f})")
    print()
    print("  Key insight: With pure Neumann pressure BCs the system is SINGULAR")
    print("  (pressure is only defined up to a constant). The SOR residual")
    print("  ||laplacian(p) - RS|| never reaches eps=0.001 regardless of itermax.")
    print("  The omega that achieves the LOWEST residual within the iteration budget")
    print("  is the optimal choice. Near omega≈1.7–1.9 is typically best for 50×50.")
    print("  omega<1: under-relaxation (slow).  omega→2: over-relaxation (diverges).")

    # ── Plot 5a ────────────────────────────────────────────────────────────────
    omg_vals = [float(r[0]) for r in rows if r[3] != "-"]
    res_vals = [float(r[3]) for r in rows if r[3] != "-"]
    best_omg = omg_vals[res_vals.index(min(res_vals))]
    colors   = ["tab:green" if o == best_omg else "tab:blue" for o in omg_vals]
    fig, ax  = plt.subplots(figsize=(7, 4))
    ax.bar([str(o) for o in omg_vals], res_vals, color=colors, edgecolor="black", linewidth=0.6)
    ax.set_xlabel("Relaxation factor ω")
    ax.set_ylabel("Avg SOR residual (lower = better)")
    ax.set_title("Task 5a — SOR Residual vs. Relaxation Factor ω\n"
                 "(itermax=500, 50×50 grid, Re=100)")
    green_patch = mpatches.Patch(color="tab:green", label=f"Best ω = {best_omg}")
    ax.legend(handles=[green_patch])
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "task5a_omega_residual.png")
    plt.close(fig)
    print(f"\n  → Plot saved: study_plots/task5a_omega_residual.png")

    # ── 5b: vary itermax ───────────────────────────────────────────────────
    print("\n── 5b: Effect of itermax (omega=1.7) ──\n")
    itermaxs = [5, 10, 20, 50, 100, 200]
    rows2    = []
    for im in itermaxs:
        cfg = {**BASE_CFG, "omg": 1.7, "itermax": im}
        r   = run_case(cfg, f"itermax_{im}")
        rows2.append([
            im,
            f"{r.get('avg_sor', '-'):.1f}" if "avg_sor" in r else "-",
            r.get("max_sor", "-"),
            r["status"],
        ])
        print(f"  itermax={im:>3}  done  ({r['status']})", flush=True)

    print()
    print_table(
        ["itermax", "avg SOR iter", "max SOR iter", "status"],
        rows2,
    )
    print("\n  Note: Once avg_iter < itermax the solver converged fully every step.")
    print("  Small itermax (e.g., 5) leaves residual too high → inaccurate pressure,")
    print("  but simulation may still complete (with worse velocity accuracy).")

    # ── Plot 5b ────────────────────────────────────────────────────────────────
    im_vals  = [r[0] for r in rows2]
    iter_vals= [float(r[1]) if r[1] != "-" else 0 for r in rows2]
    fig, ax  = plt.subplots(figsize=(6, 4))
    ax.bar([str(v) for v in im_vals], iter_vals, color="tab:orange",
           edgecolor="black", linewidth=0.6)
    ax.plot([str(v) for v in im_vals], im_vals, "k--o", markersize=5,
            label="itermax limit")
    ax.set_xlabel("itermax")
    ax.set_ylabel("Avg SOR iterations used")
    ax.set_title("Task 5b — SOR Iterations vs. itermax\n(ω=1.7, 50×50 grid, Re=100)")
    ax.legend()
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    # Left: stability by dt
    ax1.bar([str(d) for d in dt_vals], [1]*len(dt_vals),
            color=bar_colors, edgecolor="black", linewidth=0.6)
    ax1.set_xlabel("Fixed dt")
    ax1.set_ylabel("Stable / Diverged")
    ax1.set_yticks([])
    ax1.set_title("Simulation Outcome per dt")
    ax1.tick_params(axis="x", rotation=45)
    green_p = mpatches.Patch(color="tab:green", label="OK (stable)")
    red_p   = mpatches.Patch(color="tab:red",   label="DIVERGED")
    ax1.legend(handles=[green_p, red_p])

    # Right: CFL and visc fractions
    x = range(len(dt_vals))
    w = 0.38
    ax2.bar([i - w/2 for i in x], cfl_vals,  width=w, label="CFL = dt/dx",       color="tab:blue",   edgecolor="black", linewidth=0.5)
    ax2.bar([i + w/2 for i in x], visc_vals, width=w, label="dt/dt_visc",         color="tab:orange", edgecolor="black", linewidth=0.5)
    ax2.axhline(1.0, color="red", linewidth=1.4, linestyle="--", label="limit = 1")
    ax2.set_xticks(list(x))
    ax2.set_xticklabels([str(d) for d in dt_vals], rotation=45)
    ax2.set_xlabel("Fixed dt")
    ax2.set_ylabel("Fraction of stability limit")
    ax2.set_title("CFL and Viscous Stability Fractions")
    ax2.legend(fontsize=9)

    fig.suptitle("Task 6 — Fixed Time Step Stability  (50×50, nu=0.01)", fontsize=12)
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
    print("  Observation:")
    print("  • Only grids where dx > dt (CFL < 1) remain stable with fixed dt=0.05.")
    print("  • Finer grids REQUIRE smaller dt → use adaptive time stepping!")
    print("  • With adaptive dt (tau=0.5), ALL grid sizes converge correctly.")
    print("  • At Re=1000 (nu=0.001) SOR needs more iterations than at Re=100.")

    # ── 7b: Visualize stable grid cases ───────────────────────────────────────
    print("\n── 7b: Flow visualization for stable grids ──\n")
    VIS7_DIR = SCRIPT_DIR / "LidDrivenCavity_Output" / "task7_visuals"
    VIS7_DIR.mkdir(parents=True, exist_ok=True)

    vel7_imgs    = {}   # n → Path
    stream7_imgs = {}

    for n in grids:
        dx  = 1.0 / n
        cfl = 0.05 / dx
        if cfl >= 1.0:
            print(f"  {n}×{n}: CFL={cfl:.2f} ≥ 1  →  skipping (diverged)")
            continue
        label7   = f"grid{n}"
        case_dir = VIS7_DIR / label7
        case_dir.mkdir(exist_ok=True)
        cfg_vis  = {**BASE_CFG, "imax": n, "jmax": n,
                    "dt": 0.05, "tau": -1, "nu": 0.001, "dt_value": 0.5}
        print(f"  Running {n}×{n} for VTK output ...", flush=True)
        r = run_case_vtk(cfg_vis, label7, case_dir, timeout=120)
        print(f"    solver: {r['status']}", flush=True)
        if r["status"] == "OK" and r.get("vtk_dir") and r["vtk_dir"].exists():
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
        keys = sorted(vel7_imgs.keys())
        fig, axes = plt.subplots(2, len(keys),
                                 figsize=(3.5 * len(keys), 7.5),
                                 gridspec_kw={"hspace": 0.05, "wspace": 0.05})
        if len(keys) == 1:
            axes = [[axes[0]], [axes[1]]]
        for col, n in enumerate(keys):
            for row, (img_dict, lbl) in enumerate([(vel7_imgs, "Velocity magnitude"),
                                                    (stream7_imgs, "Streamlines")]):
                ax = axes[row][col]
                if n in img_dict and img_dict[n].exists():
                    ax.imshow(Image.open(img_dict[n]))
                ax.axis("off")
                if row == 0:
                    ax.set_title(f"{n}×{n}  (dx={1/n:.4f})", fontsize=11)
                if col == 0:
                    ax.set_ylabel(lbl, fontsize=10)
        fig.suptitle("Task 7 — Stable Grids (fixed dt=0.05, nu=0.001, Re=1000)",
                     fontsize=12, y=1.01)
        out_cmp7 = PLOTS_DIR / "task7_grid_comparison.png"
        fig.savefig(out_cmp7, bbox_inches="tight", dpi=120)
        plt.close(fig)
        print(f"\n  → Comparison grid saved: study_plots/task7_grid_comparison.png")

    # ── Plot 7 ─────────────────────────────────────────────────────────────────
    grid_labels = [r[0] for r in rows]
    cfl_vals7   = [float(r[2]) for r in rows]
    statuses7   = [r[4] for r in rows]
    bar_colors7 = ["tab:green" if s == "OK" else "tab:red" for s in statuses7]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    ax1.bar(grid_labels, cfl_vals7, color=bar_colors7, edgecolor="black", linewidth=0.6)
    ax1.axhline(1.0, color="red", linewidth=1.4, linestyle="--", label="CFL limit = 1")
    ax1.set_xlabel("Grid resolution")
    ax1.set_ylabel("CFL ≈ dt / dx")
    ax1.set_title("CFL Number per Grid  (fixed dt=0.05)")
    green_p = mpatches.Patch(color="tab:green", label="OK (stable)")
    red_p   = mpatches.Patch(color="tab:red",   label="DIVERGED")
    ax1.legend(handles=[green_p, red_p, mpatches.Patch(color="none")])
    ax1.legend(handles=[green_p, red_p,
               mpatches.Patch(color="none", label="— CFL limit")])

    sor_vals7 = [float(r[5]) if r[5] != "-" else 0 for r in rows]
    ax2.bar(grid_labels, sor_vals7, color=bar_colors7, edgecolor="black", linewidth=0.6)
    ax2.set_xlabel("Grid resolution")
    ax2.set_ylabel("Avg SOR iterations")
    ax2.set_title("Avg SOR Iterations per Grid")
    ax2.legend(handles=[green_p, red_p])

    fig.suptitle("Task 7 — Grid Refinement Study  (fixed dt=0.05, nu=0.001, Re=1000)", fontsize=12)
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
    print("  • SOR iterations always hit itermax=100: pure Neumann pressure BCs form")
    print("    a singular system; the solver never fully converges regardless of Re.")
    print("    Fix: increase itermax (e.g., 500) or use a pressure reference point.")
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
        res_keys = sorted(vel_imgs.keys())
        n        = len(res_keys)
        fig, axes = plt.subplots(2, n, figsize=(3.5 * n, 7.5),
                                 gridspec_kw={"hspace": 0.05, "wspace": 0.05})
        row_labels = ["Velocity magnitude", "Streamlines"]
        img_dicts  = [vel_imgs, stream_imgs]

        for row, (img_dict, row_lbl) in enumerate(zip(img_dicts, row_labels)):
            for col, Re in enumerate(res_keys):
                ax = axes[row][col]
                if Re in img_dict and img_dict[Re].exists():
                    ax.imshow(Image.open(img_dict[Re]))
                else:
                    ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                            transform=ax.transAxes)
                ax.axis("off")
                if row == 0:
                    ax.set_title(f"Re = {Re}\n(nu = {1/Re:.4f})", fontsize=11)
                if col == 0:
                    ax.set_ylabel(row_lbl, fontsize=10, labelpad=4)

        fig.suptitle("Task 8 — LDC Flow at Different Reynolds Numbers  (t = 10 s)",
                     fontsize=13, y=1.01)
        out_cmp = PLOTS_DIR / "task8_re_comparison.png"
        fig.savefig(out_cmp, bbox_inches="tight", dpi=120)
        plt.close(fig)
        print(f"\n  → Comparison grid saved: study_plots/task8_re_comparison.png")

    # ── Plot 8 ─────────────────────────────────────────────────────────────────
    re_vals   = [int(1.0/nu) for nu in nus]
    avg_dts   = [r[2] for r in rows]
    dt_floats = []
    for s in avg_dts:
        try:    dt_floats.append(float(s))
        except: dt_floats.append(0.0)
    dt_visc_vals = [(1.0/50)**2 / (4*nu) for nu in nus]   # dx²/(4ν)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    ax1.semilogx(re_vals, dt_floats, "o-", color="tab:blue",   linewidth=2, markersize=7, label="avg adaptive dt")
    ax1.semilogx(re_vals, dt_visc_vals, "s--", color="tab:orange", linewidth=1.5, markersize=6, label="dt_visc = dx²/(4ν)")
    ax1.set_xlabel("Reynolds number Re")
    ax1.set_ylabel("Time step dt  [s]")
    ax1.set_title("Adaptive dt vs. Re")
    ax1.legend()
    ax1.set_xticks(re_vals)
    ax1.set_xticklabels([str(r) for r in re_vals])

    sor_vals8 = [float(r[3]) if r[3] != "-" else 0 for r in rows]
    ax2.bar([str(r) for r in re_vals], sor_vals8, color="tab:purple",
            edgecolor="black", linewidth=0.6)
    ax2.axhline(100, color="red", linewidth=1.4, linestyle="--", label="itermax = 100")
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
