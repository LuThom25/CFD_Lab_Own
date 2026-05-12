#include "Boundary.hpp"

namespace {
constexpr double ADIABATIC_WALL_TEMPERATURE = -1.0;

double ghost_temperature(double wall_temperature, double fluid_temperature) {
    if (wall_temperature == ADIABATIC_WALL_TEMPERATURE) return fluid_temperature;
    return 2.0 * wall_temperature - fluid_temperature;
}

void apply_wall_temperature(Fields &field, const std::vector<Cell *> &cells,
                            const std::map<int, double> &wall_temperature) {
    for (auto cell : cells) {
        const auto wall_temperature_it = wall_temperature.find(cell->wall_id());
        if (wall_temperature_it == wall_temperature.end()) continue;

        const double T_wall = wall_temperature_it->second;
        const int i = cell->i();
        const int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.t(i, j) = ghost_temperature(T_wall, field.t(i + 1, j));
        }
        if (cell->is_border(border_position::LEFT)) {
            field.t(i, j) = ghost_temperature(T_wall, field.t(i - 1, j));
        }
        if (cell->is_border(border_position::TOP)) {
            field.t(i, j) = ghost_temperature(T_wall, field.t(i, j + 1));
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.t(i, j) = ghost_temperature(T_wall, field.t(i, j - 1));
        }
    }
}
} // namespace

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyFlux(Fields & /*field*/) {}

void Boundary::applyTemperature(Fields & /*field*/) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocity(Fields &field) {
    // No-slip Dirichlet BCs for all fixed (stationary) walls.
    // For each wall cell, check which sides border a fluid cell and apply:
    //   - Zero normal-face velocity at the shared face (no penetration / no-slip)
    //   - Mirror (ghost) velocity so that interpolated velocity at the wall = 0
    //
    // Staggered grid conventions:
    //   U(i,j) lives at the RIGHT face of cell (i,j)
    //   V(i,j) lives at the TOP  face of cell (i,j)
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            // Fluid is to the right (i+1,j): the shared u-face is u(i,j) → set to 0
            field.u(i, j)     = 0.0;
            // Ghost v at wall cell mirrors fluid v so that (v_wall + v_fluid)/2 = 0
            field.v(i, j)     = -field.v(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            // Fluid is to the left (i-1,j): the shared u-face is u(i-1,j) → set to 0
            field.u(i - 1, j) = 0.0;
            // Ghost v mirrors the fluid v on the left
            field.v(i, j)     = -field.v(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            // Fluid is above (i,j+1): the shared v-face is v(i,j) → set to 0
            field.v(i, j)     = 0.0;
            // Ghost u mirrors the fluid u above
            field.u(i, j)     = -field.u(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            // Fluid is below (i,j-1): the shared v-face is v(i,j-1) → set to 0
            field.v(i, j - 1) = 0.0;
            // Ghost u mirrors the fluid u below
            field.u(i, j)     = -field.u(i, j - 1);
        }
    }
}

void FixedWallBoundary::applyPressure(Fields &field) {
    // Discrete Neumann BC: dp/dn = 0 at all fixed walls.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
}

void FixedWallBoundary::applyTemperature(Fields &field) {
    apply_wall_temperature(field, _cells, _wall_temperature);
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : Boundary(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::applyVelocity(Fields &field) {
    // Moving lid BC: the top wall moves at U_wall in the x-direction (v_wall = 0).
    // For the Lid-Driven Cavity, moving wall cells are at j = jmax+1 and border
    // the fluid below (border_position::BOTTOM).
    //
    //   v(i, j-1) = 0               — no wall-normal penetration
    //   u(i, j)   = 2*U_wall - u(i,j-1)  — ghost u so interpolated u at top face = U_wall
    double U_wall = _wall_velocity.at(LidDrivenCavity::moving_wall_id);

    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::BOTTOM)) {
            // Shared v-face at the top wall: v(i, j-1) = 0 (no penetration)
            field.v(i, j - 1) = 0.0;
            // Ghost u: (u(i,j) + u(i,j-1))/2 = U_wall  =>  u(i,j) = 2*U_wall - u(i,j-1)
            field.u(i, j)     = 2.0 * U_wall - field.u(i, j - 1);
        }
    }
}

void MovingWallBoundary::applyPressure(Fields &field) {
    // Discrete Neumann BC: dp/dn = 0 at the moving wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
}

void MovingWallBoundary::applyTemperature(Fields &field) {
    apply_wall_temperature(field, _cells, _wall_temperature);
}
