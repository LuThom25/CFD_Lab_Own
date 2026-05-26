
#include <algorithm>
#include <iostream>
#include <limits>

#include "Communication.hpp"
#include "Fields.hpp"

Fields::Fields(double nu, double dt, double tau, double alpha, double beta, double gx, double gy, int imax, int jmax,
               double UI, double VI, double PI, double TI)
    : _nu(nu), _gx(gx), _gy(gy), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);
    _T = Matrix<double>(imax + 2, jmax + 2, TI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {
    // Computation of intermediate velocity fluxes F and G tilda
    for (auto cell : grid.fluid_cells()) {
        int i = cell->i();
        int j = cell->j();

        // F_tilda(i,j)
        _F(i, j) = _U(i, j) +
                   _dt * (_nu * Discretization::laplacian(_U, i, j) - Discretization::convection_u(_U, _V, i, j)) -
                   _beta * _dt * Discretization::interpolate(_T, i, j, 1, 0) * _gx;

        // G_tilda(i,j)
        _G(i, j) = _V(i, j) +
                   _dt * (_nu * Discretization::laplacian(_V, i, j) - Discretization::convection_v(_U, _V, i, j)) -
                   _beta * _dt * Discretization::interpolate(_T, i, j, 0, 1) * _gy;
    }
}

void Fields::calculate_rs(Grid &grid) {
    // Computation of the Right-Hand Side of the Pressure Poisson Equation (PPE)
    // RS(i,j) = (1/dt) * [ (F_tilda(i,j) - F_tilda(i-1,j))/dx + (G_tilda(i,j) - G_tilda(i,j-1))/dy ]
    for (auto cell : grid.fluid_cells()) {
        int i = cell->i();
        int j = cell->j();

        _RS(i, j) = (1.0 / _dt) * ((_F(i, j) - _F(i - 1, j)) / grid.dx() + (_G(i, j) - _G(i, j - 1)) / grid.dy());
    }

    // Fredholm compatibility condition: only for closed domains (all-Neumann pressure BCs).
    // If outflow cells exist, a Dirichlet pressure reference is prescribed at the boundary
    // and the mean subtraction must NOT be applied, subtracting the mean would shift the entire
    // pressure field away from that reference, forcing the SOR to spend extra iterations
    // readjusting the pressure level each timestep.
    if (grid.outflow_cells().empty()) {
        const auto &cells = grid.fluid_cells();
        const size_t N = cells.size();
        if (N == 0) return;
        double sum = 0.0;
        for (auto cell : cells)
            sum += _RS(cell->i(), cell->j());
        const double mean = sum / static_cast<double>(N);
        for (auto cell : cells)
            _RS(cell->i(), cell->j()) -= mean;
    }
}

void Fields::calculate_velocities(Grid &grid) {
    // Implementing Eq. 7 & 8: correct velocities using the updated pressure gradient
    for (auto cell : grid.fluid_cells()) {
        int i = cell->i();
        int j = cell->j();

        // Eq. 7: U(i,j) = F(i,j) - dt * (P(i+1,j) - P(i,j)) / dx
        _U(i, j) = _F(i, j) - _dt * (_P(i + 1, j) - _P(i, j)) / grid.dx();

        // Eq. 8: V(i,j) = G(i,j) - dt * (P(i,j+1) - P(i,j)) / dy
        _V(i, j) = _G(i, j) - _dt * (_P(i, j + 1) - _P(i, j)) / grid.dy();
    }
}

void Fields::calculate_temperature(Grid &grid) {
    // Start from the current field so non-fluid/ghost values are preserved.
    Matrix<double> T_next = _T;

    for (auto cell : grid.fluid_cells()) {
        int i = cell->i();
        int j = cell->j();

        double convectionT = Discretization::convection_T(_T, _U, _V, i, j);
        double laplacianT = Discretization::laplacian(_T, i, j);

        T_next(i, j) = _T(i, j) + _dt * (_alpha * laplacianT - convectionT);
    }
    _T = T_next; // Assignment operator
}
double Fields::calculate_dt(Grid &grid) {
    // Adaptive time step control based on viscous, thermal, and convective stability limits.
    // Only recompute if tau > 0 (otherwise use the fixed dt from the input file).
    if (_tau <= 0.0) {
        return _dt;
    }

    double dx = grid.dx();
    double dy = grid.dy();

    // Eq. 12: viscous stability condition
    double dt_visc = (dx * dx * dy * dy) / (2.0 * _nu * (dx * dx + dy * dy));

    // Explicit temperature diffusion stability condition. When the energy
    // equation is disabled, alpha stays zero and this condition is inactive.
    double dt_temp = (_alpha > 0.0) ? (dx * dx * dy * dy) / (2.0 * _alpha * (dx * dx + dy * dy))
                                    : std::numeric_limits<double>::max();

    // Eq. 13: convective (CFL) stability conditions — find max |u| and |v| over fluid cells
    double umax = 0.0;
    double vmax = 0.0;
    for (auto cell : grid.fluid_cells()) {
        int i = cell->i();
        int j = cell->j();
        umax = std::max(umax, std::abs(_U(i, j)));
        vmax = std::max(vmax, std::abs(_V(i, j)));
    }

    double dt_u = (umax > 0.0) ? dx / umax : std::numeric_limits<double>::max();
    double dt_v = (vmax > 0.0) ? dy / vmax : std::numeric_limits<double>::max();

    _dt = _tau * std::min({dt_visc, dt_temp, dt_u, dt_v});
    return _dt;
}

double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::t(int i, int j) { return _T(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
double Fields::tau() const { return _tau; }
double Fields::alpha() const { return _alpha; }
double Fields::beta() const { return _beta; }
