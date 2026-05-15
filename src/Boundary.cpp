#include "Boundary.hpp"

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyVelocity(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::TOP)) {
            applyVelocityTop(field, i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            applyVelocityBottom(field, i, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            applyVelocityLeft(field, i, j);
        }
        if (cell->is_border(border_position::RIGHT)) {
            applyVelocityRight(field, i, j);
        }
    }
}

void Boundary::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::TOP)) {
            applyPressureTop(field, i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            applyPressureBottom(field, i, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            applyPressureLeft(field, i, j);
        }
        if (cell->is_border(border_position::RIGHT)) {
            applyPressureRight(field, i, j);
        }
    }
}

void Boundary::applyFlux(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::TOP)) {
            applyFluxTop(field, i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            applyFluxBottom(field, i, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            applyFluxLeft(field, i, j);
        }
        if (cell->is_border(border_position::RIGHT)) {
            applyFluxRight(field, i, j);
        }
    }
}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocityTop(Fields &field, int i, int j) {
    // No-slip BC at the top wall: v(i,j) = 0 and u(i,j) = -u(i,j+1)
    field.v(i, j) = 0.0;
    field.u(i, j) = -field.u(i, j + 1);
}

void FixedWallBoundary::applyVelocityBottom(Fields &field, int i, int j) {
    // No-slip BC at the bottom wall: v(i,j-1) = 0 and u(i,j) = -u(i,j-1)
    field.v(i, j - 1) = 0.0;
    field.u(i, j) = -field.u(i, j - 1);
}

void FixedWallBoundary::applyVelocityLeft(Fields &field, int i, int j) {
    // No-slip BC at the left wall: u(i-1,j) = 0 and v(i,j) = -v(i-1,j)
    field.u(i - 1, j) = 0.0;
    field.v(i, j) = -field.v(i - 1, j);
}

void FixedWallBoundary::applyVelocityRight(Fields &field, int i, int j) {
    // No-slip BC at the right wall: u(i,j) = 0 and v(i,j) = -v(i+1,j)
    field.u(i, j) = 0.0;
    field.v(i, j) = -field.v(i + 1, j);
}

void FixedWallBoundary::applyPressureTop(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the top wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i, j + 1);
}

void FixedWallBoundary::applyPressureBottom(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the bottom wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i, j - 1);
}

void FixedWallBoundary::applyPressureLeft(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the left wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i - 1, j);
}

void FixedWallBoundary::applyPressureRight(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the right wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i + 1, j);
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : Boundary(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::applyVelocityTop(Fields &field, int i, int j) {
    // Shared v-face at the top wall: v(i,j) = 0 (no penetration)
    field.v(i, j) = 0.0;
    // Ghost u: (u(i,j) + u(i,j+1))/2 = 0  =>  u(i,j) = -u(i,j+1)
    double U_wall = _wall_velocity.at(LidDrivenCavity::moving_wall_id);
    field.u(i, j) = 2.0 * U_wall - field.u(i, j + 1);
}

void MovingWallBoundary::applyVelocityBottom(Fields &field, int i, int j) {
    // Shared v-face at the top wall: v(i, j-1) = 0 (no penetration)
    field.v(i, j - 1) = 0.0;
    // Ghost u: (u(i,j) + u(i,j-1))/2 = U_wall  =>  u(i,j) = 2*U_wall - u(i,j-1)
    double U_wall = _wall_velocity.at(LidDrivenCavity::moving_wall_id);
    field.u(i, j) = 2.0 * U_wall - field.u(i, j - 1);
}

void MovingWallBoundary::applyVelocityLeft(Fields &field, int i, int j) {
    // Shared u-face at the left wall: u(i-1,j) = 0 (no penetration)
    field.u(i - 1, j) = 0.0;
    // Ghost v: (v(i,j) + v(i-1,j))/2 = 0  =>  v(i,j) = -v(i-1,j)
    double V_wall = _wall_velocity.at(LidDrivenCavity::moving_wall_id);
    field.v(i, j) = 2.0 * V_wall - field.v(i - 1, j);
}

void MovingWallBoundary::applyVelocityRight(Fields &field, int i, int j) {
    // Shared u-face at the right wall: u(i,j) = 0 (no penetration)
    field.u(i, j) = 0.0;
    // Ghost v: (v(i,j) + v(i+1,j))/2 = 0  =>  v(i,j) = -v(i+1,j)
    double V_wall = _wall_velocity.at(LidDrivenCavity::moving_wall_id);
    field.v(i, j) = 2.0 * V_wall - field.v(i + 1, j);
}

void MovingWallBoundary::applyPressureTop(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the top wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i, j + 1);
}

void MovingWallBoundary::applyPressureBottom(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the bottom wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i, j - 1);
}

void MovingWallBoundary::applyPressureLeft(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the left wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i - 1, j);
}

void MovingWallBoundary::applyPressureRight(Fields &field, int i, int j) {
    // Discrete Neumann BC: dp/dn = 0 at the right wall.
    // Ghost cell pressure equals the pressure of the adjacent fluid cell.
    field.p(i, j) = field.p(i + 1, j);
}

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double u_in, double v_in)
    : Boundary(cells), _u_in(u_in), _v_in(v_in) {}

void InFlowBoundary::applyVelocityTop(Fields &field, int i, int j) {
    field.u(i, j) = 0.0;
    field.v(i, j) = _v_in;
}

void InFlowBoundary::applyVelocityBottom(Fields &field, int i, int j) {
    field.u(i, j) = 0.0;
    field.v(i, j - 1) = _v_in;
}

void InFlowBoundary::applyVelocityLeft(Fields &field, int i, int j) {
    field.u(i - 1, j) = _u_in;
    field.v(i, j) = 0.0;
}

void InFlowBoundary::applyVelocityRight(Fields &field, int i, int j) {
    field.u(i, j) = _u_in;
    field.v(i, j) = 0.0;
}

void InFlowBoundary::applyPressureTop(Fields &field, int i, int j) { field.p(i, j) = field.p(i, j + 1); }

void InFlowBoundary::applyPressureBottom(Fields &field, int i, int j) { field.p(i, j) = field.p(i, j - 1); }

void InFlowBoundary::applyPressureLeft(Fields &field, int i, int j) { field.p(i, j) = field.p(i - 1, j); }

void InFlowBoundary::applyPressureRight(Fields &field, int i, int j) { field.p(i, j) = field.p(i + 1, j); }

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

void OutFlowBoundary::applyVelocityTop(Fields &field, int i, int j) {
    field.u(i, j) = 0.0;
    field.v(i, j) = field.v(i, j + 1);
}

void OutFlowBoundary::applyVelocityBottom(Fields &field, int i, int j) {
    field.u(i, j) = 0.0;
    field.v(i, j - 1) = field.v(i, j - 1);
}

void OutFlowBoundary::applyVelocityLeft(Fields &field, int i, int j) {
    field.u(i - 1, j) = field.u(i - 1, j);
    field.v(i, j) = 0.0;
}

void OutFlowBoundary::applyVelocityRight(Fields &field, int i, int j) {
    field.u(i, j) = field.u(i + 1, j);
    field.v(i, j) = 0.0;
}

void OutFlowBoundary::applyPressureTop(Fields &field, int i, int j) { field.p(i, j) = 0.0; }

void OutFlowBoundary::applyPressureBottom(Fields &field, int i, int j) { field.p(i, j) = 0.0; }

void OutFlowBoundary::applyPressureLeft(Fields &field, int i, int j) { field.p(i, j) = 0.0; }

void OutFlowBoundary::applyPressureRight(Fields &field, int i, int j) { field.p(i, j) = 0.0; }
