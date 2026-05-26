# CFD Lab Worksheet 2: Arbitrary geometries and heat transport.

---
The main implementations have been:
1. Including different derived classes for the boundary conditions, as well as ensuring their correct application within the obstacles and heat transport cases.
2. Including heat transport in the system of equations having T as an additional field, and updating it in 'calculate_temperature' inside 'Fields'. Rewriting the intermediate fluxes and add the temperature discretization inside 'Discretization'.

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
├── example_cases/            # all folders contain the .dat files and output files
│   ├── ChannelWithBFS
│   ├── ChannelWithObstacle             
│   ├── FluidTrap         
│   ├── LidDrivenCavity 
│   ├── NaturalConvection  
│   ├── PlaneShearFlow
│   └── RayleighBenard   
│       
├── docs/
│   └── first-steps.md
├── build/                          # CMake build directory (git-ignored)
└── CMakeLists.txt
```

---


## Compile

```bash
cd ~/CFD_Lab/cfd-lab-group-c
cd build
cmake ..
make 
```

---

## Run for the different examples, separately or together:

./fluidchen ../example_cases/ChannelWithBFS/ChannelWithBFS.dat
./fluidchen ../example_cases/ChannelWithObstacle/ChannelWithObstacle.dat
./fluidchen ../example_cases/FluidTrap/FluidTrap.dat
./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat
./fluidchen ../example_cases/NaturalConvection/NaturalConvection_a.dat
./fluidchen ../example_cases/NaturalConvection/NaturalConvection_b.dat
./fluidchen ../example_cases/PlaneShearFlow/PlaneShearFlow.dat
./fluidchen ../example_cases/RayleighBenard/RayleighBenard.dat

```
