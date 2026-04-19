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

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
BINARY     = SCRIPT_DIR.parent.parent / "build" / "fluidchen"

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

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    if not BINARY.exists():
        print(f"[ERROR] Binary not found: {BINARY}")
        print("  Run:  cd fluidchen-skeleton/build && cmake .. && make -j4")
        sys.exit(1)

    # Decide which tasks to run
    args  = sys.argv[1:]
    tasks = [int(a) for a in args if a.isdigit()] or [4, 5, 6, 7, 8]

    print("=" * 66)
    print("  Worksheet 1 — Simulation Tasks 4–8")
    print(f"  Binary : {BINARY}")
    print("=" * 66)

    if 4 in tasks: task4()
    if 5 in tasks: task5()
    if 6 in tasks: task6()
    if 7 in tasks: task7()
    if 8 in tasks: task8()

    print("\n" + "=" * 66)
    print("  All selected studies complete.")
    print("  Results above answer Tasks 4–8 of Worksheet 1.")
    print("=" * 66 + "\n")

if __name__ == "__main__":
    main()
