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
/Applications/ParaView-6.1.0.app/Contents/bin/pvpython visualize.py
```
Generates static images (u, v, pressure, velocity magnitude, glyphs, streamlines, Jet-coloured vectors, direction-only white arrows) and animation frames for all quantities. All output goes to `task4/`.

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
cd fluidchen-skeleton/build && make -j4 && \
cd ../example_cases/LidDrivenCavity && \
../../build/fluidchen LidDrivenCavity.dat && \
/Applications/ParaView-6.1.0.app/Contents/bin/pvpython visualize.py && \
bash make_videos.sh && \
python3 run_studies.py
```

## Branches

| Branch | Content |
|--------|---------|
| `ws1` | Core solver implementation (Tasks 1–4) |
| `ws1_extensions` | Parameter studies, run_studies.py (Tasks 5–8) |
| `ws1_further_extensions` | Vector visualizations, README, Summary, Report |
| `ws1_further_extensions_improved_SOR` | SOR convergence fix (Fredholm + zero-mean pressure) |

## Key Results Summary

| Task | Study | Key Finding |
|------|-------|------------|
| 4 | Base LDC, Re=100 | Stable vortex at t≈5 s; avg_dt=5×10⁻³ s; **avg SOR=4.8 iter** (Fredholm fix) |
| 5a | SOR omega ∈ [0.5, 1.99] | Pre-fix: best omega≈1.5 (lowest residual). Post-fix: omega≈1.7–1.9 fastest (~5 iter) |
| 5b | itermax ∈ [5, 200] | Pre-fix: always hits cap. Post-fix: converges within budget; low itermax degrades accuracy |
| 6 | Fixed dt stability | dt_visc=0.010 is binding; stable for dt < 0.010 |
| 7 | Grid 16–256, fixed dt=0.05 | 16×16 OK; **32×32 marginal** (CFL=1.6, survives t_end=5 s); 64×64+ diverges |
| 8 | nu=0.01→0.0001 (Re=100→10000) | avg_dt increases with Re (viscous limit relaxes) |

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
- ParaView ≥ 5.x (for `pvpython visualize.py`)
- Python 3 with: `matplotlib`, `Pillow`
- ffmpeg (for video assembly)
