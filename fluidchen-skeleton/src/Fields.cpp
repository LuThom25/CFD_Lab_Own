
#include <algorithm>
#include <iostream>
#include <limits>

#include "Communication.hpp"
#include "Fields.hpp"

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {
    // Implementing Eq. 9 & 10: compute intermediate velocity fluxes F and G
    // using explicit Euler time integration. F and G incorporate diffusion (Laplacian)
    // and convection (donor-cell scheme) terms from the momentum equations.
    for (auto cell : grid.fluid_cells()) {
        int i = cell->i();
        int j = cell->j();

        // Eq. 9: F(i,j) = U + dt * (nu * laplacian(U) - convection_u + gx)
        _F(i, j) = _U(i, j) + _dt * (
            _nu * Discretization::laplacian(_U, i, j)
            - Discretization::convection_u(_U, _V, i, j)
            + _gx
        );

        // Eq. 10: G(i,j) = V + dt * (nu * laplacian(V) - convection_v + gy)
        _G(i, j) = _V(i, j) + _dt * (
            _nu * Discretization::laplacian(_V, i, j)
            - Discretization::convection_v(_U, _V, i, j)
            + _gy
        );
    }
}

void Fields::calculate_rs(Grid &grid) {
    // Implementing Eq. 11: compute the Right-Hand Side of the Pressure Poisson Equation (PPE)
    // RS(i,j) = (1/dt) * [ (F(i,j) - F(i-1,j))/dx + (G(i,j) - G(i,j-1))/dy ]
    for (auto cell : grid.fluid_cells()) {
        int i = cell->i();
        int j = cell->j();

        _RS(i, j) = (1.0 / _dt) * (
            (_F(i, j) - _F(i - 1, j)) / grid.dx()
          + (_G(i, j) - _G(i, j - 1)) / grid.dy()
        );
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

double Fields::calculate_dt(Grid &grid) {
    // Implementing Eqs. 12 & 13: adaptive time step control based on three stability criteria.
    // Only recompute if tau > 0 (otherwise use the fixed dt from the input file).
    if (_tau <= 0.0) {
        return _dt;
    }

    double dx = grid.dx();
    double dy = grid.dy();

    // Eq. 12: viscous stability condition
    double dt_visc = (dx * dx * dy * dy) / (2.0 * _nu * (dx * dx + dy * dy));

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

    _dt = _tau * std::min({dt_visc, dt_u, dt_v});
    return _dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
