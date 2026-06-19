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

// ---------------------------------------------------------------------------
// CG_Solver
// ---------------------------------------------------------------------------

// Auxiliary vectors are sized identically to the pressure matrix (size_x+2,
// size_y+2) so that the Laplacian stencil can be applied directly via (i,j).
// nc and nr are obtained from field.p_matrix().num_cols()/num_rows() in
// Case.cpp, where the grid is already fully constructed before the solver.
CG_Solver::CG_Solver(double tolerance, int max_iter, int nc, int nr)
    : _tolerance(tolerance), _max_iter(max_iter),
      _r(nc, nr, 0.0), _d(nc, nr, 0.0), _q(nc, nr, 0.0) {}

// Local dot product over fluid_cells() — no MPI reduce.
// The global reduce is performed explicitly inside iterate() where alpha and
// beta are computed, matching the pattern in Case.cpp for the residual reduce.
double CG_Solver::dot_local(const Matrix<double> &a, const Matrix<double> &b,
                             const std::vector<Cell *> &cells) {
    double s = 0.0;
    for (auto cell : cells)
        s += a(cell->i(), cell->j()) * b(cell->i(), cell->j());
    return s;
}

// AXPY: y(i,j) += alpha * x(i,j)  over fluid_cells()
void CG_Solver::axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                     const std::vector<Cell *> &cells) {
    for (auto cell : cells)
        y(cell->i(), cell->j()) += alpha * x(cell->i(), cell->j());
}

void CG_Solver::iterate(Fields &field, Grid &grid) {
    const auto &cells = grid.fluid_cells();
    const double n_local  = static_cast<double>(cells.size());
    const double n_global = Communication::reduce_sum(n_local);

    // --- Initial residual: r = Δₕp - RS  (= b_cg - A_cg*p, A_cg = -Δₕ) ---
    // Using the pressure field as-is (warm start): p already holds the solution
    // from the previous timestep, which is close to the new solution and
    // reduces CG iterations.
    for (auto cell : cells) {
        int i = cell->i(), j = cell->j();
        _r(i, j) = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
    }

    // Initial search direction d = r
    for (auto cell : cells) {
        int i = cell->i(), j = cell->j();
        _d(i, j) = _r(i, j);
    }

    // α = (r,r)/(d,q) must be identical on every MPI rank — a rank-local dot
    // product would give a different α per rank and break the algorithm.
    // reduce_sum aggregates the partial sums into one global value.
    // SOR does not need this because each cell update is purely local:
    // p(i,j) depends only on its neighbors, never on a globally shared scalar.
    double rr = Communication::reduce_sum(dot_local(_r, _r, cells));

    // --- CG iteration loop ---
    // Stopping criterion uses the same normalised L2-norm as SOR (via
    // calculate_residual below), but checked internally so that iterate()
    // delivers a converged p before Case.cpp reapplies BCs.
    _last_iter_count = 0;
    for (int iter = 0; iter < _max_iter && std::sqrt(rr / n_global) > _tolerance; ++iter) {

        // Halo of d must be current before the stencil is applied.
        // Case.cpp communicates p (not d), so this call stays here.
        // With iproc=jproc=1 this is a no-op.
        Communication::communicate_field(_d, grid.domain());

        // Matvec: q = A_cg * d = -Δₕd
        // Sign convention: A_cg = -Δₕ is positive semidefinite → CG converges.
        for (auto cell : cells) {
            int i = cell->i(), j = cell->j();
            _q(i, j) = -Discretization::laplacian(_d, i, j);
        }

        // Step size: α = (r,r)_global / (d,q)_global
        double dq    = Communication::reduce_sum(dot_local(_d, _q, cells));
        double alpha = rr / dq;

        // p += α d,   r -= α q
        axpy( alpha, _d, field.p_matrix(), cells);
        axpy(-alpha, _q, _r,               cells);

        // Conjugacy coefficient and new search direction: d = r_new + β d
        double rr_new = Communication::reduce_sum(dot_local(_r, _r, cells));
        double beta   = rr_new / rr;
        for (auto cell : cells) {
            int i = cell->i(), j = cell->j();
            _d(i, j) = _r(i, j) + beta * _d(i, j);
        }
        rr = rr_new;
        ++_last_iter_count;
    }
}

// Identical to SOR_Standard::calculate_residual(): recomputes the local L2-norm
// from the current pressure field. Case.cpp performs the global MPI reduce and
// checks convergence — no change to Case.cpp needed.
// No mean-pressure subtraction: the Fredholm correction on RS in
// Fields::calculate_rs() is sufficient, matching SOR behaviour.
double CG_Solver::calculate_residual(Fields &field, Grid &grid) {
    double rloc = 0.0;
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    return std::sqrt(rloc / static_cast<double>(grid.fluid_cells().size()));
}

// ---------------------------------------------------------------------------
// PCG_SSOR — Preconditioned CG with Red-Black SSOR preconditioner
// ---------------------------------------------------------------------------

PCG_SSOR::PCG_SSOR(double tolerance, int max_iter, int nc, int nr)
    : _tolerance(tolerance), _max_iter(max_iter),
      _r(nc, nr, 0.0), _d(nc, nr, 0.0), _q(nc, nr, 0.0), _z(nc, nr, 0.0) {}

double PCG_SSOR::dot_local(const Matrix<double> &a, const Matrix<double> &b,
                            const std::vector<Cell *> &cells) {
    double s = 0.0;
    for (auto c : cells)
        s += a(c->i(), c->j()) * b(c->i(), c->j());
    return s;
}

void PCG_SSOR::axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                    const std::vector<Cell *> &cells) {
    for (auto c : cells)
        y(c->i(), c->j()) += alpha * x(c->i(), c->j());
}

void PCG_SSOR::apply_ssor(const Matrix<double> &r, Matrix<double> &z, const Grid &grid) {
    const double dx   = grid.dx(), dy = grid.dy();
    const double d_ii = 2.0 / (dx * dx) + 2.0 / (dy * dy);

    for (auto c : _red_cells)   z(c->i(), c->j()) = 0.0;
    for (auto c : _black_cells) z(c->i(), c->j()) = 0.0;

    // Forward sweep: red then black
    for (auto c : _red_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    Communication::communicate_field(z, grid.domain());

    for (auto c : _black_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    Communication::communicate_field(z, grid.domain());

    // Backward sweep: red only (black is skipped — it is a no-op).
    // Proof: backward-black reads only red neighbours via sor_helper. Those
    // red values were set in Pass 1 and have not changed since (Pass 2 only
    // writes black cells). So backward-black would produce the same z_black
    // as Pass 2, wasting 1 cell loop + 1 communicate_field per apply_ssor
    // call (= 2 MPI operations per outer PCG iteration).
    for (auto c : _red_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    Communication::communicate_field(z, grid.domain());
}

void PCG_SSOR::iterate(Fields &field, Grid &grid) {
    const auto  &cells    = grid.fluid_cells();
    const double n_local  = static_cast<double>(cells.size());
    const double n_global = Communication::reduce_sum(n_local);

    // Build red/black cell lists once — coloring uses global (i,j), so it is
    // identical across all MPI decompositions, guaranteeing M = Mᵀ.
    if (!_cells_partitioned) {
        for (auto c : cells) {
            if ((c->i() + c->j()) % 2 == 0) _red_cells.push_back(c);
            else                              _black_cells.push_back(c);
        }
        _cells_partitioned = true;
    }

    // Warm start: r = Δₕp - RS
    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _r(i, j) = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
    }

    apply_ssor(_r, _z, grid);

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _d(i, j) = _z(i, j);
    }

    // ρ = (r,z) — M⁻¹-weighted dot product (replaces (r,r) from plain CG)
    double rz = Communication::reduce_sum(dot_local(_r, _z, cells));
    double rr = Communication::reduce_sum(dot_local(_r, _r, cells));

    _last_iter_count = 0;
    for (int iter = 0; iter < _max_iter && std::sqrt(rr / n_global) > _tolerance; ++iter) {

        Communication::communicate_field(_d, grid.domain());

        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _q(i, j) = -Discretization::laplacian(_d, i, j);
        }

        const double dq    = Communication::reduce_sum(dot_local(_d, _q, cells));
        const double alpha = rz / dq;

        axpy( alpha, _d, field.p_matrix(), cells);
        axpy(-alpha, _q, _r,               cells);

        apply_ssor(_r, _z, grid);

        const double rz_new = Communication::reduce_sum(dot_local(_r, _z, cells));
        const double beta   = rz_new / rz;

        // d = z + β·d  (z, not r — key PCG difference from plain CG)
        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _d(i, j) = _z(i, j) + beta * _d(i, j);
        }

        rz = rz_new;
        rr = Communication::reduce_sum(dot_local(_r, _r, cells));
        ++_last_iter_count;
    }
}

double PCG_SSOR::calculate_residual(Fields &field, Grid &grid) {
    double rloc = 0.0;
    for (auto c : grid.fluid_cells()) {
        int i = c->i(), j = c->j();
        const double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += val * val;
    }
    return std::sqrt(rloc / static_cast<double>(grid.fluid_cells().size()));
}
