#include <cmath>

#include "Communication.hpp"
#include "PressureSolver.hpp"

SOR_Standard::SOR_Standard(double omega) : _omega(omega) {}

// We divide into update method and residual calculation method. This way the residual 
// can be computed with updated BCs and halo values. Previously it also worked, since the difference
// was probaly negligible. Now that we can have way more non-updated cells thanks to the halo, 
// it may become relevant. 

void SOR_Standard::iterate(Fields &field, Grid &grid) {
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
}

double SOR_Standard::calculate_residual(Fields &field, Grid &grid) {
    double res = 0.0;  // residual initialization
    double rloc = 0.0; // accumulates val^2 for every cell

    // Cummulative sum of squared residuals (Laplacian(p(i,j)) - rhs(i,j))^2 over all cells
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    {
        res = rloc / static_cast<double>(grid.fluid_cells().size()); // mean of squared residuals
        res = std::sqrt(res);                                        // L2 norm of the residual
    }

    return res;
}
