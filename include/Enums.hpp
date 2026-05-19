#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 3: fixed wall, 4: moving wall
namespace LidDrivenCavity {
// CQ5: constexpr ensures compile-time evaluation and ODR-safe linkage.
constexpr int moving_wall_id = 8;
constexpr int fixed_wall_id = 4;
constexpr double wall_velocity = 1.0;
} // namespace LidDrivenCavity

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

enum class cell_type {
    FLUID,
    FIXED_WALL,
    MOVING_WALL,
    INFLOW,  // PGM value 1: Dirichlet velocity BC, Neumann pressure BC
    OUTFLOW, // PGM value 2: Neumann velocity BC, Dirichlet pressure BC (p = 0)
    DEFAULT,
};
