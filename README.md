# CFD Lab — Worksheet 1: Lid-Driven Cavity

A 2D incompressible Navier-Stokes solver implementing the **Chorin Projection (fractional-step) method** on a staggered Cartesian grid, with SOR pressure solver and donor-cell convection scheme.

## Problem Setup

The lid-driven cavity is a classic CFD benchmark: a unit square filled with viscous fluid, driven by a moving top wall at velocity U=1 m/s. Three walls are fixed (no-slip), the top wall moves horizontally.

```
  ←──── U_lid = 1 m/s ────→
  ┌───────────────────────┐  ← moving wall (Γ_lid)
  │                       │
  │   Viscous fluid Ω     │
  │   Re = U·L/ν          │
  │                       │
  └───────────────────────┘
  ↑ no-slip on 3 walls (Γ_no-slip)
```

## Repository Layout

```
CFD_Lab_Own/
├── README.md                 ← this file
├── CLAUDE.md                 ← session context & quick reference
├── Summary.md                ← comprehensive theory & code guide
├── Report.md                 ← academic lab report (Tasks 1–8)
└── fluidchen-skeleton/
    ├── src/                  ← C++ implementation
    ├── include/              ← headers
    ├── build/                ← compiled binary (after cmake)
    └── example_cases/LidDrivenCavity/
        ├── LidDrivenCavity.dat   ← simulation config
        ├── visualize.py          ← ParaView visualization pipeline
        ├── make_videos.sh        ← ffmpeg video assembly
        ├── run_studies.py        ← Tasks 4–8 parameter studies
        └── LidDrivenCavity_Output/
            ├── task4/            ← all Task-4 output
            │   ├── final_*.png   ← static result images
            │   ├── frames/       ← animation frames
            │   └── videos/       ← MP4 animations
            ├── study_plots/      ← parameter study plots (Tasks 5–8)
            ├── task7_visuals/    ← Task 7 grid-refinement snapshots
            └── task8_visuals/    ← Task 8 Re-number snapshots
```

## Build & Run

### 1. Build the solver
```bash
cd fluidchen-skeleton
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

### 2. Run the base simulation (Task 4)
```bash
cd example_cases/LidDrivenCavity
../../build/fluidchen LidDrivenCavity.dat
```
Expected output: solver header → progress lines every 0.5 s → SUMMARY line.  
Runtime: ~2 min for t_end=50 at 50×50, Re=100.

### 3. Visualize results (requires ParaView)
```bash
/Applications/ParaView-6.0.1.app/Contents/bin/pvpython visualize.py
```
Generates static images (u, v, pressure, velocity magnitude, glyphs, streamlines, Jet-coloured vectors, direction-only white arrows) and animation frames for all quantities. All output goes to `task4/`.

> **Note:** `run_studies.py` uses temporary directories and cleans up all VTK output after each run. Always run `fluidchen LidDrivenCavity.dat` directly before calling `visualize.py`.

### 4. Create videos
```bash
bash make_videos.sh
```
Produces 6 individual videos (`video_u.mp4`, `video_v.mp4`, `video_p.mp4`, `video_vel.mp4`, `video_vec.mp4`, `video_vec_clean.mp4`) and one combined `video_panel_2x3.mp4` — all in `task4/videos/`.

### 5. Run all parameter studies (Tasks 4–8)
```bash
python3 run_studies.py          # all tasks
python3 run_studies.py 5        # single task
python3 run_studies.py 5 6 7 8  # multiple tasks
```
Outputs ASCII tables + matplotlib plots saved to `study_plots/`.

### 6. Run everything at once
```bash
cd fluidchen-skeleton/example_cases/LidDrivenCavity && \
../../build/fluidchen LidDrivenCavity.dat && \
/Applications/ParaView-6.0.1.app/Contents/bin/pvpython visualize.py && \
bash make_videos.sh && \
python3 run_studies.py 5 6 7 8
```

## Branches

| Branch | Content |
|--------|---------|
| `ws1` | Core solver implementation (Tasks 1–4) |
| `ws1_extensions` | Parameter studies, run_studies.py (Tasks 5–8) |
| `ws1_further_extensions` | Vector visualizations, README, Summary, Report |
| `ws1_further_extensions_improved_SOR` | SOR convergence fix (Fredholm + zero-mean pressure) |
| `ws1_further_extensions_improved_SOR_small_improvements` | **Current**: code quality, plot improvements, bug fixes |

## Key Results Summary

| Task | Study | Key Finding |
|------|-------|------------|
| 4 | Base LDC, Re=100, t_end=50 | Quasi-steady vortex; avg_dt=5×10⁻³ s; avg SOR=4.8 iter (Fredholm fix) |
| 5a | SOR ω ∈ [0.5, 1.99], itermax=500, t_end=50 | **Best ω=1.9** (8 avg iter, never hits itermax); ω=1.7 faster on avg (5 iter) but occasionally hits cap |
| 5b | itermax ∈ [5, 500], ω=1.9, t_end=50 | itermax=5 **DIVERGES**; itermax≥10 stable; itermax≥200 reaches ε=0.001 |
| 6 | Fixed dt stability, t_end=5 | dt_visc=0.010 s is binding; stable only for dt < 0.010 s |
| 7 | Grid 16–256, fixed dt=0.05, nu=0.001, t_end=5 | Only 16×16 truly stable (CFL=0.80); 32×32 marginal; 64×64+ diverges. SOR iter: 15.7→35.7→100 |
| 8 | nu=0.01→0.0001 (Re=100→10000), t_end=50 | avg_dt increases with Re (viscous limit relaxes); complex vortex structures at high Re |

## Output Files Overview

All static images and videos are in `LidDrivenCavity_Output/task4/`.

| File | Description |
|------|-------------|
| `task4/final_u.png` | u-velocity, Blue-to-Red, range [-0.5, 1.0] |
| `task4/final_v.png` | v-velocity, Blue-to-Red, range [-0.25, 0.25] |
| `task4/final_pressure.png` | Pressure (relative), Cool-to-Warm, auto-scale |
| `task4/final_velocity.png` | Velocity magnitude, Jet, range [0, 1] |
| `task4/final_glyphs.png` | Arrow glyphs coloured by magnitude |
| `task4/final_streamlines.png` | Streamlines on dark background |
| `task4/final_vectors_bw.png` | Jet-coloured vectors on dark background |
| `task4/final_vectors_clean.png` | **Direction-only** white arrows on navy background |
| `task4/videos/video_vec.mp4` | Jet-coloured vector animation |
| `task4/videos/video_vec_clean.mp4` | Direction-only vector animation |
| `task4/videos/video_panel_2x3.mp4` | Combined panel: u, v, p / vel, vec, vec_clean |
| `study_plots/task8_re_comparison.png` | Flow comparison Re=100/500/2000/10000 |

## Dependencies

- C++17 compiler (clang/gcc)
- CMake ≥ 3.14
- VTK (bundled via CMake)
- ParaView ≥ 5.x with pvpython (tested: 6.0.1, 6.1.0)
- Python 3 with: `matplotlib`, `Pillow`
- ffmpeg (for video assembly)
