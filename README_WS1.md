# CFD Lab Worksheet 1: Lid-Driven Cavity

In this worksheet the 2D incompressible Navier-Stokes equations in the lid-driven were solved. In our setup, the top wall moves at $(u,v)=(1,0)$ velocity, and the remaining three walls are no-slip. Fields (u, v, p) are stored on a staggered Cartesian grid and advanced in time with explicit Euler + SOR pressure solve.

---

## Repository layout

```
cfd-lab-group-c/
├── src/
│   ├── main.cpp
│   ├── Case.cpp              # Main simulation loop
│   ├── Boundary.cpp          # Wall boundary conditions
│   ├── Discretization.cpp    # FD stencils (convection, diffusion, interpolation)
│   ├── Fields.cpp            # Fluxes, RHS, velocity update, adaptive dt
│   ├── Grid.cpp              # Grid construction and cell access
│   ├── Cell.cpp              # Cell type logic
│   ├── Communication.cpp     # MPI stubs (unused in WS1)
│   └── PressureSolver.cpp    # SOR variants
├── include/
│   ├── Boundary.hpp
│   ├── Case.hpp
│   ├── Cell.hpp
│   ├── Communication.hpp
│   ├── Datastructures.hpp
│   ├── Discretization.hpp
│   ├── Domain.hpp
│   ├── Enums.hpp             # solver_type, cell_type, border_position enums
│   ├── Fields.hpp
│   ├── Grid.hpp
│   └── PressureSolver.hpp
│   └── Util.hpp              # added for the future    
├── example_cases/
│   └── LidDrivenCavity/
│       ├── LidDrivenCavity.dat     # Simulation parameters
│       ├── visualize.py            # ParaView automation (Task 4)
│       ├── run_studies.py          # Parameter sweep runner (Tasks 5–8)
│       ├── make_videos.sh          # ffmpeg video assembly
│       └── LidDrivenCavity_Output/ # All output (git-ignored)
├── docs/
│   └── first-steps.md
├── build/                          # CMake build directory (git-ignored)
└── CMakeLists.txt
```

---

## Implementation

The main implementation of the SOR alsorithm is in `Case::simulate()`, which calls functions from `Fields`, `Discretization`, and `Boundary`. To get the full SOR algorithm working, we needed to complete:

### `src/Case.cpp`
- `Case::simulate()`: the main time-stepping loop: selects adaptive `dt`, applies boundary conditions, computes F/G and the pressure RHS, runs SOR iterations until convergence, updates velocities, and writes VTK output

### `src/Fields.cpp`
- `calculate_fluxes()`: computes the intermediate velocity fluxes F and G
- `calculate_rs()`: copmutes the right-hand side of the pressure Poisson equation
- `calculate_velocities()`: updates u and v after the pressure solve
- `calculate_dt()`: computes next iteration's adaptive time step from the CFL and diffusion stability conditions

### `src/Discretization.cpp`
- `convection_u()`, `convection_v()`: donor-cell convective terms
- `laplacian()`: second-order diffusion stencil
- `interpolate()`: interpolates staggered values to the required locations

### `src/Boundary.cpp`
Both `FixedWallBoundary` and `MovingWallBoundary` implement:
- `applyVelocity()`: sets ghost-cell values to enforce no-slip (fixed walls) or `(u,v) = (1,0)` (moving lid)
- `applyPressure()`: zero-gradient Neumann condition via ghost-cell copy
- `applyFlux()`: sets F/G on boundary cells

Beyond the strictly required implementations, several existing files were also lightly modified to improve code quality and performance. Loop orders in `Grid.cpp` were swapped to match column-major memory layout, magic number literals were replaced with named constants, and error handling was added for malformed geometry files. Performance-oriented changes included unchecked `fast()` accessors in Datastructures.hpp, precomputed reciprocals in Discretization.hpp, and compiler hint macros in the new Util.hpp to reduce overhead in the hot SOR and convection loops. Additionally, const-correctness fixes, constexpr upgrades, and the solver_type enum were introduced to make the codebase more robust and extensible beyond what Worksheet 1 explicitly required.

---

## Pressure solver

Three SOR variants are available and selected directly in the `.dat` file — no recompilation needed:

| Solver | `.dat` value | Description |
|--------|-------------|-------------|
| Standard SOR | `SOR_STANDARD` | Classic row-by-row SOR sweep |
| Red-Black SOR | `SOR_RB` | Checkerboard ordering, better parallelism |
| SOR with mean correction | `SOR_MEAN_CORRECTION` | Adds a mean pressure correction each iteration |

Set it in `LidDrivenCavity.dat`:
```
solver = SOR_STANDARD
```
We chose SOR_STANDARD as the solver for all our simulations. 

The remaining solver parameters are also set in the `.dat` file:

| Parameter | `.dat` keyword | Effect |
|-----------|---------------|--------|
| Relaxation factor | `omg` | ω ∈ (0, 2); 1.7 is a good default |
| Max iterations | `itermax` | SOR stops here even if not converged |
| Residual tolerance | `eps` | SOR stops early when residual < eps |

---

## Visualization and analysis

To automate the ParaView visualization and analysis tasks (sweeping over the relaxation parameter, maximum number of iterations, etc.), we created two Python files which both live in `example_cases/LidDrivenCavity/`.

### `visualize.py` — Task 4 (ParaView)
Invoked via `pvpython`, it loads the VTK output and saves:
- Scalar field plots: `final_u.png`, `final_v.png`, `final_pressure.png`, `final_velocity.png`
- Glyph (arrow) plot: `final_glyphs.png`
- Streamline plot: `final_streamlines.png`
- Animation frames in `frames/` (assembled into videos by `make_videos.sh`)

```bash
"/mnt/c/Program Files/ParaView 6.0.1/bin/pvpython.exe" visualize.py
# or on native Linux:
pvpython visualize.py
```

### `run_studies.py` — Tasks 5–8 (parameter sweeps)
Reruns the simulation with varying parameters and saves comparison plots to `LidDrivenCavity_Output/study_plots/`.

| Task | What is swept | Key parameter |
|------|--------------|--------------|
| 5 | SOR relaxation factor ω | `omg` ∈ (0, 2), effect on convergence speed and iteration count |
| 6 | Fixed time step `dt` | Adaptive stepping disabled; stability boundary identified |
| 7 | Grid resolution | `imax = jmax` ∈ {16, 32, 64, 128, 256} with fixed `dt` |
| 8 | Kinematic viscosity | `nu` ∈ {0.01, 0.002, 0.0005, 0.0001} → Re ∈ {100, 500, 2000, 10000} |

Run all tasks at once or individually:
```bash
python3 run_studies.py 5 6 7 8   # all
python3 run_studies.py 5         # ω sweep only
python3 run_studies.py 6         # fixed dt stability only
python3 run_studies.py 7         # grid refinement only
python3 run_studies.py 8         # viscosity / Re study only
```

---

## Dependencies

| Tool | Purpose |
|------|---------|
| CMake ≥ 3.12 + GCC | Build system |
| ParaView 6.0.1 (`pvpython`) | Visualization |
| Python 3 + `matplotlib` + `pillow` | Study plots |
| `ffmpeg` | Video assembly |

Install Python and ffmpeg in WSL:
```bash
sudo apt install python3-matplotlib python3-pil ffmpeg
```

---

## Compile

```bash
cd ~/CFD_Lab/cfd-lab-group-c
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

The binary is at `build/fluidchen`. Re-run `make -j$(nproc)` after any source change.

---

## Run

### Everything at once (Tasks 4–8)

```bash
cd ~/CFD_Lab/cfd-lab-group-c/example_cases/LidDrivenCavity \
  && ../../build/fluidchen LidDrivenCavity.dat \
  && "/mnt/c/Program Files/ParaView 6.0.1/bin/pvpython.exe" visualize.py \
  && bash make_videos.sh \
  && python3 run_studies.py 5 6 7 8
```

> On native Linux replace the `pvpython.exe` path with just `pvpython`.
> If ParaView is elsewhere: `find /mnt/c -maxdepth 4 -iname "pvpython.exe"`.

### Step by step

| Step | Command | What it does |
|------|---------|-------------|
| Simulate | `../../build/fluidchen LidDrivenCavity.dat` | Writes VTK files to `LidDrivenCavity_Output/` |
| Visualize | `pvpython visualize.py` | Renders PNGs + animation frames (Task 4) |
| Videos | `bash make_videos.sh` | Assembles frames into `.mp4` |
| Studies | `python3 run_studies.py <tasks>` | Runs parameter sweeps and saves plots |

### Running specific tasks only

```bash
python3 run_studies.py 5        # ω sweep (SOR convergence)
python3 run_studies.py 6        # fixed dt stability study
python3 run_studies.py 7        # grid refinement study
python3 run_studies.py 8        # varying viscosity / Reynolds number
python3 run_studies.py 5 6 7 8  # all of the above
```

---

## Output

All results are written to `example_cases/LidDrivenCavity/LidDrivenCavity_Output/`:

```
LidDrivenCavity_Output/
├── task4/
│   ├── final_{u,v,pressure,velocity,glyphs,streamlines}.png
│   ├── frames/       # PNG frames for animation
│   └── videos/       # .mp4 files
└── study_plots/      # Plots for Tasks 5–8
```
