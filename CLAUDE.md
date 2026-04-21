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
            ├── *.vtk                  ← simulation output
            ├── final_*.png            ← static final-state images
            ├── frames/                ← animation frames (u,v,p,vel,vec)
            ├── videos/                ← MP4 videos
            └── study_plots/           ← parameter study plots (Tasks 5–8)
```

### Active Branches
| Branch | Purpose |
|--------|---------|
| `ws1` | Core implementation (Tasks 1–4) |
| `ws1_extensions` | Parameter studies (Tasks 5–8) + run_studies.py |
| `ws1_further_extensions` | Vector visuals + README/Summary/Report |

### Key Design Decisions
1. **Pressure null space**: Pure Neumann BCs → singular PPE. SOR never converges to eps=0.001. Solution: per-frame pressure rescaling for visualization; avg_res metric to compare omega values.
2. **Adaptive dt**: `tau=-1` disables adaptive stepping (used for Tasks 6 & 7 stability studies).
3. **SUMMARY line**: Case::simulate() outputs a machine-readable `SUMMARY ...` line parsed by run_studies.py regex.
4. **Visualization style**: Dark background (0.12/0.18), Jet colormap 0..1 for velocity, parallel camera at (0.5,0.5), 1200×1000 px.

### Quick Commands
```bash
# Build
cd fluidchen-skeleton/build && make -j4

# Run simulation
cd example_cases/LidDrivenCavity
../../build/fluidchen LidDrivenCavity.dat

# Visualize (ParaView)
/Applications/ParaView-6.1.0.app/Contents/bin/pvpython visualize.py

# Make videos
bash make_videos.sh

# Parameter studies (Tasks 4–8)
python3 run_studies.py          # all tasks
python3 run_studies.py 5 6 7 8  # specific tasks
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
