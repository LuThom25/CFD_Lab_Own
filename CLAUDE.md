# CLAUDE.md — Session Context & Project Logic

## Project: CFD Lab Worksheet 1 — Lid-Driven Cavity

### Repository Structure
```
CFD_Lab_Own/
├── README.md                          ← project overview & run instructions
├── CLAUDE.md                          ← this file: session context
├── Summary.md                         ← theory & code master guide
├── Report.md                          ← academic lab report
└── fluidchen-skeleton/
    ├── build/fluidchen                ← compiled binary
    ├── src/
    │   ├── Case.cpp                   ← simulate() main loop
    │   ├── Fields.cpp                 ← F, G, RS, dt, velocities
    │   ├── Discretization.cpp         ← laplacian, convection_u/v, interpolate
    │   └── Boundary.cpp               ← no-slip & moving-wall BCs
    ├── include/
    │   ├── Case.hpp / Fields.hpp      ← class declarations
    │   └── ...
    └── example_cases/LidDrivenCavity/
        ├── LidDrivenCavity.dat        ← simulation parameters
        ├── visualize.py               ← pvpython visualization pipeline
        ├── make_videos.sh             ← ffmpeg video assembly
        ├── run_studies.py             ← Tasks 4–8 parameter studies
        └── LidDrivenCavity_Output/
            ├── *.vtk                  ← simulation output (root level)
            ├── task4/                 ← all Task-4 output
            │   ├── final_u.png        ← static final-state images (u,v,p,vel,vec,vec_clean)
            │   ├── frames/            ← animation frames (u,v,p,vel,vec,vec_clean)
            │   └── videos/            ← 6 individual MP4s + video_panel_2x3.mp4
            ├── study_plots/           ← parameter study plots (Tasks 5–8)
            ├── task7_visuals/         ← Task 7 grid-refinement snapshots
            └── task8_visuals/         ← Task 8 Re-number snapshots
```

### Active Branches
| Branch | Purpose |
|--------|---------|
| `ws1` | Core implementation (Tasks 1–4) |
| `ws1_extensions` | Parameter studies (Tasks 5–8) + run_studies.py |
| `ws1_further_extensions` | Vector visuals + README/Summary/Report |
| `ws1_further_extensions_improved_SOR` | SOR convergence fix (Fredholm + zero-mean pressure) |

### Key Design Decisions
1. **Pressure null space — FIXED**: Pure Neumann BCs → singular PPE. Root causes: (a) Fredholm incompatibility (`Σ RS ≠ 0`) and (b) null-space drift (constant accumulates in p each sweep). Fix: subtract `mean(RS)` after `calculate_rs()` and subtract `mean(p)` after each SOR sweep. Result: avg_sor 100→4.8, avg_res 1.23→3.4×10⁻³ (converges to eps=0.001). Visualization: per-frame pressure rescaling still used.
2. **Adaptive dt**: `tau=-1` disables adaptive stepping (used for Tasks 6 & 7 stability studies).
3. **SUMMARY line**: Case::simulate() outputs a machine-readable `SUMMARY ...` line parsed by run_studies.py regex.
4. **Visualization style**: Dark background (0.12/0.18), Jet colormap 0..1 for velocity, parallel camera at (0.5,0.5), 1200×1000 px. `OrientationAxesVisibility=0` suppresses the axes widget globally.
5. **Clean vector plot** (`vec_clean`): Uniform-length white arrows on navy background (Calculator constant=1.0 for magnitude), Stride=3, ScaleFactor=0.060. Shows direction only without Jet color overlay.
6. **CFL marginal stability**: 32×32 at CFL=1.60 (fixed dt, t_end=5s) survives because local CFL starts at 0 and instability grows too slowly (~100 steps). 64×64 at CFL=3.2 diverges immediately. Task 7 marks marginal cases with `[marginal]` label.
7. **openvkl warning**: `Could not find a module for device type "cpu"` is harmless — optional Intel ray-tracing library not installed; ParaView falls back gracefully.

### Full Pipeline (run everything from scratch)
```bash
# 1. Build
cd fluidchen-skeleton/build && make -j4

# 2. Run base simulation (Task 4)
cd ../example_cases/LidDrivenCavity
../../build/fluidchen LidDrivenCavity.dat

# 3. Visualize: static images + animation frames
/Applications/ParaView-6.1.0.app/Contents/bin/pvpython visualize.py

# 4. Assemble videos
bash make_videos.sh

# 5. Parameter studies (Tasks 5–8): runs solver + generates plots
python3 run_studies.py          # all tasks
python3 run_studies.py 5 6 7 8  # specific tasks
```

### Quick Commands
```bash
# Build only
cd fluidchen-skeleton/build && make -j4

# Run simulation only
cd example_cases/LidDrivenCavity && ../../build/fluidchen LidDrivenCavity.dat

# Visualize (ParaView pvpython)
/Applications/ParaView-6.1.0.app/Contents/bin/pvpython visualize.py

# Make videos (requires visualize.py frames first)
bash make_videos.sh

# Parameter studies
python3 run_studies.py 5 6 7 8
```

### Physical Parameters (Base Case)
| Parameter | Value | Meaning |
|-----------|-------|---------|
| `imax`, `jmax` | 50, 50 | Grid cells |
| `nu` | 0.01 | Kinematic viscosity → Re = U·L/ν = 100 |
| `t_end` | 50.0 | End time [s] |
| `tau` | 0.5 | Adaptive dt safety factor |
| `omg` | 1.7 | SOR relaxation factor |
| `eps` | 0.001 | SOR convergence tolerance |
| `itermax` | 100 | Max SOR iterations per step |
| `gamma` | 0.5 | Donor-cell upwind coefficient |
