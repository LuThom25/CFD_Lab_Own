#include <cmath>

#include "Communication.hpp"
#include "PressureSolver.hpp"

// THOMAS IMPLEMENTATION OF 'PressureSolver.cpp'
SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> & /*boundaries*/) {

    double dx    = grid.dx();
    double dy    = grid.dy();
    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = ω·h²/4 when dx==dy==h

    // ── SOR sweep ────────────────────────────────────────────────────────────────
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

    // ── Zero-mean pressure projection ────────────────────────────────────────────
    // The all-Neumann Poisson system is singular: pressure is defined only up to
    // a constant (null space = constant vectors).  Each SOR sweep drifts the
    // solution along the null space, preventing convergence.  Subtracting the
    // mean after every sweep projects the iterate onto the space orthogonal to
    // the null space, where the discrete Poisson operator IS invertible.
    const auto &cells = grid.fluid_cells();
    const std::size_t N = cells.size();
    if (N > 0) {
        double p_sum = 0.0;
        for (auto cell : cells) p_sum += field.p(cell->i(), cell->j());
        const double p_mean = p_sum / static_cast<double>(N);
        for (auto cell : cells) field.p(cell->i(), cell->j()) -= p_mean;
    }

    // ── Residual (L2 norm on the zero-mean field) ─────────────────────────────────
    double rloc = 0.0;
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += val * val;
    }
    return std::sqrt(rloc / static_cast<double>(N > 0 ? N : 1));
}

// ICIAR'S IMPLEMENTATION OF 'PressureSolver.cpp'
SOR::SOR_Iciar(double omega) : _omega(omega) {}

double SOR_Iciar::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double dx = grid.dx();
    double dy = grid.dy();

    // Pre-computed coefficient for the update formula
    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h
    
    // Iteration of the SOR going from it to it+1
    // Loop over all fluid cells to compute p(i,j) at it+1 from it values
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

    double res = 0.0;   // residual initialization
    double rloc = 0.0;  // accumulates val^2 for every cell

    // Cummulative sum of squared residuals (Laplacian(p(i,j)) - rhs(i,j))^2 over all cells
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    {
        res = rloc / (grid.fluid_cells().size()); // mean of squared residuals
        res = std::sqrt(res);                     // L2 norm of the residual
    }

    return res;
}

