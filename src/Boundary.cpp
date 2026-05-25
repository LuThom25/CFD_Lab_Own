#include "Boundary.hpp"

inline constexpr double ADIABATIC_WALL_TEMPERATURE = -1.0;

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
        double ghost_temperature_sum = 0.0;
        int border_count = 0;

        if (cell->is_border(border_position::RIGHT)) {
            ghost_temperature_sum += ghost_temperature(T_wall, field.t(i + 1, j));
            ++border_count;
        }
        if (cell->is_border(border_position::LEFT)) {
            ghost_temperature_sum += ghost_temperature(T_wall, field.t(i - 1, j));
            ++border_count;
        }
        if (cell->is_border(border_position::TOP)) {
            ghost_temperature_sum += ghost_temperature(T_wall, field.t(i, j + 1));
            ++border_count;
        }
        if (cell->is_border(border_position::BOTTOM)) {
            ghost_temperature_sum += ghost_temperature(T_wall, field.t(i, j - 1));
            ++border_count;
        }

        if (border_count > 0) {
            field.t(i, j) = ghost_temperature_sum / static_cast<double>(border_count);
        }
    }
}

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyVelocity(Fields &field) {
    for (auto cell : _cells) {
        applyVelocityToCell(field, cell);
    }
}

void Boundary::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        applyPressureToCell(field, cell);
    }
}

void Boundary::applyFlux(Fields &field) {
    for (auto cell : _cells) {
        applyFluxToCell(field, cell);
    }
}

// This is implemented as an empty function in the base class, since not all boundaries need to apply temperature
// boundary conditions. For example, in the current test cases, only the fixed wall boundaries have temperature boundary
// conditions.
void Boundary::applyTemperature([[maybe_unused]] Fields &field) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocityToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();
    const bool top = cell->is_border(border_position::TOP);
    const bool bottom = cell->is_border(border_position::BOTTOM);
    const bool left = cell->is_border(border_position::LEFT);
    const bool right = cell->is_border(border_position::RIGHT);

    if (top) field.v(i, j) = 0.0;
    if (bottom) field.v(i, j - 1) = 0.0;
    if (left) field.u(i - 1, j) = 0.0;
    if (right) field.u(i, j) = 0.0;

    // At obstacle corners, tangential ghost values can share a staggered slot
    // with the normal velocity of the other wall face. Preserve no-penetration.
    if (top && !right) {
        field.u(i, j) = -field.u(i, j + 1);
    }
    if (bottom && !right) {
        field.u(i, j) = -field.u(i, j - 1);
    }
    if (left && !top) {
        field.v(i, j) = -field.v(i - 1, j);
    }
    if (right && !top) {
        field.v(i, j) = -field.v(i + 1, j);
    }
}

void FixedWallBoundary::applyPressureToCell(Fields &field, Cell *cell) {
    // This code does not support fixed wall cells with one ore two borders. More edges will not cause memory errors,
    // but the pressure will be calculated as the average of the neighboring pressures, which is not correct for fixed
    // walls. However, this is not expected to happen in the current test cases.
    const int i = cell->i(), j = cell->j();
    double pressure = 0.0;
    int count = 0;
    constexpr std::array<double, 4> count_inv = {1.0 / 1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0};

    if (cell->is_border(border_position::TOP)) {
        pressure += field.p(i, j + 1);
        count++;
    }
    if (cell->is_border(border_position::BOTTOM)) {
        pressure += field.p(i, j - 1);
        count++;
    }
    if (cell->is_border(border_position::LEFT)) {
        pressure += field.p(i - 1, j);
        count++;
    }
    if (cell->is_border(border_position::RIGHT)) {
        pressure += field.p(i + 1, j);
        count++;
    }
    field.p(i, j) = pressure * count_inv[count - 1];
}

void FixedWallBoundary::applyTemperature(Fields &field) { apply_wall_temperature(field, _cells, _wall_temperature); }

// For each wall cell, we make zero the flux at the face shared with a fluid neighbour.
// Without this, the pressure solver uses leftover flux values from the previous
// timestep at those faces, which causes wrong pressures and eventually divergence.
void FixedWallBoundary::applyFluxToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();
    if (cell->is_border(border_position::TOP)) field.g(i, j) = 0.0;
    if (cell->is_border(border_position::BOTTOM)) field.g(i, j - 1) = 0.0;
    if (cell->is_border(border_position::LEFT)) field.f(i - 1, j) = 0.0;
    if (cell->is_border(border_position::RIGHT)) field.f(i, j) = 0.0;
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

void MovingWallBoundary::applyVelocityToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();
    const bool top = cell->is_border(border_position::TOP);
    const bool bottom = cell->is_border(border_position::BOTTOM);
    const bool left = cell->is_border(border_position::LEFT);
    const bool right = cell->is_border(border_position::RIGHT);
    const double U_wall = _wall_velocity.at(LidDrivenCavity::moving_wall_id);
    const double V_wall = _wall_velocity.at(LidDrivenCavity::moving_wall_id);

    if (top) field.v(i, j) = 0.0;
    if (bottom) field.v(i, j - 1) = 0.0;
    if (left) field.u(i - 1, j) = 0.0;
    if (right) field.u(i, j) = 0.0;

    if (top && !right) {
        field.u(i, j) = 2.0 * U_wall - field.u(i, j + 1);
    }
    if (bottom && !right) {
        field.u(i, j) = 2.0 * U_wall - field.u(i, j - 1);
    }
    if (left && !top) {
        field.v(i, j) = 2.0 * V_wall - field.v(i - 1, j);
    }
    if (right && !top) {
        field.v(i, j) = 2.0 * V_wall - field.v(i + 1, j);
    }
}

void MovingWallBoundary::applyPressureToCell(Fields &field, Cell *cell) {
    // This code does not support moving wall cells with one ore two borders. More edges will not cause memory errors,
    // but the pressure will be calculated as the average of the neighboring pressures, which is not correct for moving
    // walls. However, this is not expected to happen in the current test cases.
    const int i = cell->i(), j = cell->j();
    double pressure = 0.0;
    int count = 0;
    constexpr std::array<double, 4> count_inv = {1.0 / 1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0};

    if (cell->is_border(border_position::TOP)) {
        pressure += field.p(i, j + 1);
        count++;
    }
    if (cell->is_border(border_position::BOTTOM)) {
        pressure += field.p(i, j - 1);
        count++;
    }
    if (cell->is_border(border_position::LEFT)) {
        pressure += field.p(i - 1, j);
        count++;
    }
    if (cell->is_border(border_position::RIGHT)) {
        pressure += field.p(i + 1, j);
        count++;
    }
    field.p(i, j) = pressure * count_inv[count - 1];
}

void MovingWallBoundary::applyFluxToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();
    if (cell->is_border(border_position::TOP)) field.g(i, j) = 0.0;
    if (cell->is_border(border_position::BOTTOM)) field.g(i, j - 1) = 0.0;
    if (cell->is_border(border_position::LEFT)) field.f(i - 1, j) = 0.0;
    if (cell->is_border(border_position::RIGHT)) field.f(i, j) = 0.0;
}

void MovingWallBoundary::applyTemperature(Fields &field) { apply_wall_temperature(field, _cells, _wall_temperature); }

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double u_in, double v_in)
    : Boundary(cells), _u_in(u_in), _v_in(v_in) {}

void InFlowBoundary::applyVelocityToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();

    if (cell->is_border(border_position::TOP)) {
        field.u(i, j) = 0.0;
        field.v(i, j) = _v_in;
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = 0.0;
        field.v(i, j - 1) = _v_in;
    }
    if (cell->is_border(border_position::LEFT)) {
        field.u(i - 1, j) = _u_in;
        field.v(i, j) = 0.0;
    }
    if (cell->is_border(border_position::RIGHT)) {
        field.u(i, j) = _u_in;
        field.v(i, j) = 0.0;
    }
}

void InFlowBoundary::applyPressureToCell(Fields &field, Cell *cell) {
    // This code does not support inflow cells with more than one border, but this is not expected to happen in the
    // current test cases
    const int i = cell->i(), j = cell->j();

    if (cell->is_border(border_position::TOP)) {
        field.p(i, j) = field.p(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.p(i, j) = field.p(i, j - 1);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.p(i, j) = field.p(i - 1, j);
    }
    if (cell->is_border(border_position::RIGHT)) {
        field.p(i, j) = field.p(i + 1, j);
    }
}

void InFlowBoundary::applyFluxToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();

    if (cell->is_border(border_position::TOP)) {
        field.g(i, j) = _v_in;
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.g(i, j - 1) = _v_in;
    }
    if (cell->is_border(border_position::LEFT)) {
        field.f(i - 1, j) = _u_in;
    }
    if (cell->is_border(border_position::RIGHT)) {
        field.f(i, j) = _u_in;
    }
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

void OutFlowBoundary::applyVelocityToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();

    if (cell->is_border(border_position::TOP)) {
        field.u(i, j) = field.u(i, j + 1);
        field.v(i, j) = field.v(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = field.u(i, j - 1);
        field.v(i, j - 1) = field.v(i, j - 2);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.u(i - 1, j) = field.u(i - 2, j);
        field.v(i, j) = field.v(i - 1, j);
    }
    if (cell->is_border(border_position::RIGHT)) {
        field.u(i, j) = field.u(i + 1, j);
        field.v(i, j) = field.v(i + 1, j);
    }
}

void OutFlowBoundary::applyPressureToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();

    field.p(i, j) = 0.0;
}

void OutFlowBoundary::applyFluxToCell(Fields &field, Cell *cell) {
    const int i = cell->i(), j = cell->j();

    if (cell->is_border(border_position::TOP)) {
        field.g(i, j) = field.v(i, j);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.g(i, j - 1) = field.v(i, j - 1);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.f(i - 1, j) = field.u(i - 1, j);
    }
    if (cell->is_border(border_position::RIGHT)) {
        field.f(i, j) = field.u(i, j);
    }
}
