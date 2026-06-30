<div align="center">
  <img width="466" height="492" src="FluidchenLogo.png">
</div>

Fluidchen is a CFD Solver developed for the CFD Lab taught at TUM Informatics, Chair of Scientific Computing in Computer Science.

After forking, use this `README.md` however you want: feel free to remove anything you don't need,
or add any additional details we should know to run the code.

## Working with fluidchen

You will extend this code step-by-step starting from a pure framework to a parallel CFD solver. Please follow these [instructions for work with git and submitting the assignments](docs/first-steps.md).

## Software Requirements

This code is known to work on all currently supported Ubuntu LTS versions (22.04, 20.04).
In particular, you will need:

- A recent version of the GCC compiler (>9.1). Other compilers should also work, but you may need to tweak the CMakeLists.txt file (contributions welcome). GCC 9.3, 11.2 and 13.2 are known to work.
- CMake, to configure the build.
- The VTK library, to generate result files. libvtk7 and libvtk9 are known to work.
- OpenMPI (not for the skeleton, but when you implement parallelization).

Get the dependencies on Ubuntu:

```shell
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install build-essential cmake libvtk7-qt-dev openmpi-bin libopenmpi-dev
```

## Building the code

```shell
git clone https://gitlab.lrz.de/oguzziya/GroupX_CFDLab.git
cd GroupX_CFDLab
mkdir build && cd build
cmake ..
make
```

After `make` completes successfully, you will get an executable `fluidchen` in your `build` directory. Run the code from this directory.

### Build options

By default, **fluidchen** is installed in `DEBUG` mode. To obtain full performance, you can execute cmake as

```shell
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
```

or

```shell
cmake -DCMAKE_CXX_FLAGS="-O3" ..
```

You can see and modify all CMake options with, e.g., `ccmake .` inside `build/` (Ubuntu package `cmake-curses-gui`).

A good idea would be that you setup your computers as runners for [GitLab CI](https://docs.gitlab.com/ee/ci/)
(see the file `.gitlab-ci.yml` here) to check the code building automatically every time you push.

## Running

In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. Navigate to the `build/` directory and run:

```shell
./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `../example_cases/LidDrivenCavity/LidDrivenCavity_Output`, which holds the `.vtk` files of the solution.

If the input file does not contain a geometry file (added later in the course), fluidchen will run the lid-driven cavity case with the given parameters.

## Special systems

### macOS

Students have previously successfully ran **Fluidchen** on macOS. However, this is not a system we regularly test. Your contribution here is, therefore, essential in improving this documentation.

In macOS, you can use the default `clang` compiler. Do not install `gcc`, since it might cause problems with the standard library and VTK. Other dependencies can be installed by using `homebrew` package manager as

```shell
brew install cmake
brew install open-mpi
brew install vtk
```

#### macOS Troubleshooting

- In macOS, the default `g++` command is linked to `clang++` command, which means, `g++` command does not run the GCC compiler but the Clang compiler.
- Setup of GCC compiler is expected to be cumbersome and clashes with lots of other dependencies, therefore please do not use GCC compiler on this project.
- If CMake cannot find the correct C++ binary, you can set it by

  ```shell
  export CXX=`which clang++``
  export CMAKE_CXX_COMPILER=`which clang++``
  ```

which is going to set the corresponding environment variables to the path of Clang compiler. Please note that if you run these commands on a terminal session, they are only going to be valid on that terminal session. In order to make these changes permanent, you can add these lines to your `~/.zshrc` file.

- Although installation of VTK looks simple, sometimes it is possible that CMake cannot find some necessary libraries for VTK, most famous one being Qt5. If you face an error something like:

```shell
CMake Error at /usr/local/lib/cmake/vtk-9.0/VTK-vtk-module-find-packages.cmake:115 (find_package):
 By not providing "FindQt5.cmake" in CMAKE_MODULE_PATH this project has
 asked CMake to find a package configuration file provided by "Qt5", but
 CMake did not find one.

 Could not find a package configuration file provided by "Qt5" (requested
 version 5.15) with any of the following names:

   Qt5Config.cmake
   qt5-config.cmake

 Add the installation prefix of "Qt5" to CMAKE_PREFIX_PATH or set "Qt5_DIR"
 to a directory containing one of the above files.  If "Qt5" provides a
 separate development package or SDK, be sure it has been installed.
```

which means that CMake could not find Qt5. Solution is simple fortunately. First, make sure that you have Qt5 installed:

```shell
brew install qt5
```

Then extend `CMAKE_PREFIX_PATH`, which are the locations where CMake tries to find packages, by adding following lines to your `.zshrc` file

```shell
export CMAKE_PREFIX_PATH="/usr/local/opt/qt5:$CMAKE_PREFIX_PATH"
```

Please not that your installation might be in another location. The most possible another location is `/usr/local/Cellar/qt@5/5.**.*/`, which depends on the Qt5 version.

## Troubleshooting
### VTK not found

You might run into a problem where the VTK library is not found. To fix this, you can try the following steps:

1. Find the installation path of your VTK library
2. Define this path as an environment variable, as e.g. `export VTK_DIR=".../lib/cmake/vtk-8.2"`
3. Start in a clean build folder
4. Run `cmake ..` again

### No rule to make target '/usr/lib/x86_64-linux-gnu/libdl.so'

We are investigating an [issue](https://gitlab.lrz.de/tum-i05/public/fluidchen-skeleton/-/issues/3) that appears on specific systems and combinations of dependencies.

---

## Case Study: NACA 0012 Airfoil Wind Tunnel

A custom case study simulating flow past a **NACA 0012 airfoil** in a 2D wind tunnel channel at **Re = 200**, **α = 0°**.

### Geometry

The airfoil geometry is rasterised onto a 200×80 Cartesian grid using an included Python script. Four angles of attack are provided:

| File | α | Description |
|------|---|-------------|
| `AirfoilWindTunnel.pgm` | 0° | Symmetric, zero lift |
| `AirfoilWindTunnel_alpha5.pgm` | 5° | Mild lift |
| `AirfoilWindTunnel_alpha10.pgm` | 10° | Moderate lift |
| `AirfoilWindTunnel_alpha15.pgm` | 15° | Near stall |

Generate a PGM for any angle of attack:
```shell
python3 tools/generate_airfoil_pgm.py --alpha 0
```

### Domain & Physical Parameters

| Parameter | Value |
|-----------|-------|
| Domain | 4.0 m × 2.0 m (4c × 2c) |
| Grid | 200 × 80 cells |
| Chord length | c = 1.0 m |
| U∞ | 1.0 m/s |
| ν | 0.005 m²/s |
| Re = U∞c/ν | **200** |
| t_end | 50 s |
| dt_value (VTK output) | 5.0 s |

### Running the Airfoil Case

```shell
# Build first (see above)
cd build

# Serial run — EIGEN_CG solver (recommended, ~64 min wall time)
mpirun -n 1 ./fluidchen ../example_cases/AirfoilWindTunnel/AirfoilWindTunnel_Re200_eigen_cg_1_1.dat

# Parallel run — PCG_SSOR solver, 2×2 MPI decomposition
mpirun -n 4 ./fluidchen ../example_cases/AirfoilWindTunnel/AirfoilWindTunnel_Re200_pcg_ssor_2_2.dat
```

Progress is logged to the console (`[step N] t=... iters=... residual=...`) and to `sor_log.csv` in the output directory.

### Analysis & Validation

A Python analysis script generates a 6-page PDF report:

```shell
python3 tools/analyze_airfoil.py \
  example_cases/AirfoilWindTunnel/AirfoilWindTunnel_Re200_eigen_cg_1_1_Output_1_1
```

The report includes:
- Flow field overview (pressure, velocity magnitude, u, v — all on a consistent blue→red scale)
- Close-up view around the airfoil with streamlines
- Full convergence history (residual and CG iterations over all 8,000 time steps)
- Solver cost analysis (iteration distribution, wall-clock breakdown)
- Validation against literature (Cp distribution, CD/CL table)

**Pre-computed results** (EIGEN_CG, Re=200, α=0°) are included in `example_cases/AirfoilWindTunnel/AirfoilWindTunnel_Re200_eigen_cg_1_1_Output_1_1/`:
- `sor_log.csv` — convergence log (8,000 steps, t=0→50)
- `airfoil_analysis.pdf` — full analysis report
- VTK snapshots at t = 0, 5, 10, 15, 20, 25, 30, 35, 40, 45 s

### Validation Results (Re=200, α=0°)

| Quantity | Simulation | Literature | Note |
|----------|-----------|------------|------|
| CD (pressure) | **0.257** | ~0.25–0.30 | Pressure drag only |
| CL | **0.000** | 0.000 | Exact symmetry at α=0° |
| Steady state from | t ≈ 0.6 s | expected | Re=200 → no vortex shedding |
| CG iter/step (median) | **0** | — | Warm-start converges instantly |
| CG iter/step (mean) | **1.7** | — | 99.4% of steps: 0 iterations |

### Pressure Solvers Available

| Solver key | Parallelism | Notes |
|------------|-------------|-------|
| `EIGEN_CG` | Serial (1×1) | Recommended for accuracy |
| `EIGEN_PCG_IC` | Serial (1×1) | Incomplete Cholesky preconditioner |
| `PCG_SSOR` | MPI parallel | Use with 2×2 or 4×1 decomposition |
| `CG_STANDARD` | MPI parallel | Standard CG |
| `SOR_STANDARD` | MPI parallel | Slow convergence on fine grids |
