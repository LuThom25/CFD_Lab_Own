#include <cmath>
#include <mpi.h>

#include "Communication.hpp"
#include "PressureSolver.hpp"

#ifdef USE_EIGEN
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#endif

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

    // Forward red
    for (auto c : _red_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    Communication::communicate_field(z, grid.domain());

    // Forward black
    for (auto c : _black_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    Communication::communicate_field(z, grid.domain());

    // Backward red (Fix A: backward-black omitted — no-op)
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

    // Fix C: compute rz and rr together in one MPI_Allreduce (saves 1 AllReduce).
    double buf0[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
    MPI_Allreduce(MPI_IN_PLACE, buf0, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double rz = buf0[0];
    double rr = buf0[1];

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

        // Fix C: rz_new and rr_new in one MPI_Allreduce (saves 1 AllReduce/iter).
        double buf[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
        MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double rz_new = buf[0];
        rr = buf[1];
        const double beta = rz_new / rz;

        // d = z + β·d  (z, not r — key PCG difference from plain CG)
        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _d(i, j) = _z(i, j) + beta * _d(i, j);
        }

        rz = rz_new;
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

// ---------------------------------------------------------------------------
// PCG_SSOR_Serial — no-MPI version for pure-compute benchmarking
// Algorithmically identical to PCG_SSOR; all communicate_field and
// MPI_Allreduce calls omitted. Correct only for single-rank runs.
// ---------------------------------------------------------------------------

PCG_SSOR_Serial::PCG_SSOR_Serial(double tolerance, int max_iter, int nc, int nr)
    : _tolerance(tolerance), _max_iter(max_iter),
      _r(nc, nr, 0.0), _d(nc, nr, 0.0), _q(nc, nr, 0.0), _z(nc, nr, 0.0) {}

double PCG_SSOR_Serial::dot_local(const Matrix<double> &a, const Matrix<double> &b,
                                   const std::vector<Cell *> &cells) {
    double s = 0.0;
    for (auto c : cells)
        s += a(c->i(), c->j()) * b(c->i(), c->j());
    return s;
}

void PCG_SSOR_Serial::axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                            const std::vector<Cell *> &cells) {
    for (auto c : cells)
        y(c->i(), c->j()) += alpha * x(c->i(), c->j());
}

void PCG_SSOR_Serial::apply_ssor(const Matrix<double> &r, Matrix<double> &z, const Grid &grid) {
    const double dx   = grid.dx(), dy = grid.dy();
    const double d_ii = 2.0 / (dx * dx) + 2.0 / (dy * dy);

    for (auto c : _red_cells)   z(c->i(), c->j()) = 0.0;
    for (auto c : _black_cells) z(c->i(), c->j()) = 0.0;

    // Forward red — no communicate after (serial: halo already consistent)
    for (auto c : _red_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    // Forward black
    for (auto c : _black_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    // Backward red (backward-black omitted — no-op)
    for (auto c : _red_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) / d_ii;
    }
    // No communicate calls anywhere in this function
}

void PCG_SSOR_Serial::iterate(Fields &field, Grid &grid) {
    const auto  &cells   = grid.fluid_cells();
    const double n_total = static_cast<double>(cells.size());

    if (!_cells_partitioned) {
        for (auto c : cells) {
            if ((c->i() + c->j()) % 2 == 0) _red_cells.push_back(c);
            else                              _black_cells.push_back(c);
        }
        _cells_partitioned = true;
    }

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _r(i, j) = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
    }

    apply_ssor(_r, _z, grid);

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _d(i, j) = _z(i, j);
    }

    // No MPI_Allreduce — dot products are global (correct for 1-rank run)
    double rz = dot_local(_r, _z, cells);
    double rr = dot_local(_r, _r, cells);

    _last_iter_count = 0;
    for (int iter = 0; iter < _max_iter && std::sqrt(rr / n_total) > _tolerance; ++iter) {

        // No communicate_field(_d) — halo consistent on single rank
        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _q(i, j) = -Discretization::laplacian(_d, i, j);
        }

        const double dq    = dot_local(_d, _q, cells);
        const double alpha = rz / dq;

        axpy( alpha, _d, field.p_matrix(), cells);
        axpy(-alpha, _q, _r,               cells);

        apply_ssor(_r, _z, grid);

        const double rz_new = dot_local(_r, _z, cells);
        rr                  = dot_local(_r, _r, cells);
        const double beta   = rz_new / rz;

        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _d(i, j) = _z(i, j) + beta * _d(i, j);
        }

        rz = rz_new;
        ++_last_iter_count;
    }
}

double PCG_SSOR_Serial::calculate_residual(Fields &field, Grid &grid) {
    double rloc = 0.0;
    for (auto c : grid.fluid_cells()) {
        int i = c->i(), j = c->j();
        const double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += val * val;
    }
    return std::sqrt(rloc / static_cast<double>(grid.fluid_cells().size()));
}

// ---------------------------------------------------------------------------
// PCG_Jacobi — diagonal (Jacobi) preconditioner
// ---------------------------------------------------------------------------

PCG_Jacobi::PCG_Jacobi(double tolerance, int max_iter, int nc, int nr)
    : _tolerance(tolerance), _max_iter(max_iter),
      _r(nc, nr, 0.0), _d(nc, nr, 0.0), _q(nc, nr, 0.0), _z(nc, nr, 0.0) {}

double PCG_Jacobi::dot_local(const Matrix<double> &a, const Matrix<double> &b,
                              const std::vector<Cell *> &cells) {
    double s = 0.0;
    for (auto c : cells)
        s += a(c->i(), c->j()) * b(c->i(), c->j());
    return s;
}

void PCG_Jacobi::axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                      const std::vector<Cell *> &cells) {
    for (auto c : cells)
        y(c->i(), c->j()) += alpha * x(c->i(), c->j());
}

void PCG_Jacobi::apply_jacobi(const Matrix<double> &r, Matrix<double> &z,
                               const std::vector<Cell *> &cells) {
    // M = D_A: z(i,j) = r(i,j) / d_ii = r(i,j) * _d_inv
    // No communicate_field needed — diagonal preconditioner is purely local.
    for (auto c : cells)
        z(c->i(), c->j()) = r(c->i(), c->j()) * _d_inv;
}

void PCG_Jacobi::iterate(Fields &field, Grid &grid) {
    const auto  &cells    = grid.fluid_cells();
    const double n_local  = static_cast<double>(cells.size());
    const double n_global = Communication::reduce_sum(n_local);

    if (!_setup_done) {
        const double dx = grid.dx(), dy = grid.dy();
        _d_inv     = 1.0 / (2.0 / (dx * dx) + 2.0 / (dy * dy));
        _setup_done = true;
    }

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _r(i, j) = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
    }

    apply_jacobi(_r, _z, cells);

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _d(i, j) = _z(i, j);
    }

    double buf[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
    MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double rz = buf[0];
    double rr = buf[1];

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

        apply_jacobi(_r, _z, cells);    // no MPI in apply_jacobi

        double buf2[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
        MPI_Allreduce(MPI_IN_PLACE, buf2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double rz_new = buf2[0];
        rr                  = buf2[1];
        const double beta   = rz_new / rz;

        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _d(i, j) = _z(i, j) + beta * _d(i, j);
        }
        rz = rz_new;
        ++_last_iter_count;
    }
}

double PCG_Jacobi::calculate_residual(Fields &field, Grid &grid) {
    double rloc = 0.0;
    for (auto c : grid.fluid_cells()) {
        int i = c->i(), j = c->j();
        const double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += val * val;
    }
    return std::sqrt(rloc / static_cast<double>(grid.fluid_cells().size()));
}

// ---------------------------------------------------------------------------
// PCG_RBLU — Incomplete Cholesky IC(0) with Red-Black ordering
// ---------------------------------------------------------------------------

PCG_RBLU::PCG_RBLU(double tolerance, int max_iter, int nc, int nr)
    : _tolerance(tolerance), _max_iter(max_iter),
      _r(nc, nr, 0.0), _d(nc, nr, 0.0), _q(nc, nr, 0.0), _z(nc, nr, 0.0) {}

double PCG_RBLU::dot_local(const Matrix<double> &a, const Matrix<double> &b,
                            const std::vector<Cell *> &cells) {
    double s = 0.0;
    for (auto c : cells)
        s += a(c->i(), c->j()) * b(c->i(), c->j());
    return s;
}

void PCG_RBLU::axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                    const std::vector<Cell *> &cells) {
    for (auto c : cells)
        y(c->i(), c->j()) += alpha * x(c->i(), c->j());
}

void PCG_RBLU::apply_rblu(const Matrix<double> &r, Matrix<double> &z, const Grid &grid) {
    for (auto c : _red_cells)   z(c->i(), c->j()) = 0.0;
    for (auto c : _black_cells) z(c->i(), c->j()) = 0.0;

    // Pass 1: forward red — z_R = (r_R + nb_sum) * _d_inv
    // sor_helper = 0 on first pass since z=0
    for (auto c : _red_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) * _d_inv;
    }
    Communication::communicate_field(z, grid.domain());

    // Pass 2: forward black — Schur complement: z_B = (r_B + nb_sum) * _s_inv
    // Key difference vs apply_ssor: _s_inv (< _d_inv) → stronger scaling
    for (auto c : _black_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) * _s_inv;
    }
    Communication::communicate_field(z, grid.domain());

    // Pass 3: backward red — z_R updated using (now current) black neighbours
    // Fix A: backward-black omitted (no-op, same as apply_ssor)
    for (auto c : _red_cells) {
        int i = c->i(), j = c->j();
        z(i, j) = (r(i, j) + Discretization::sor_helper(z, i, j)) * _d_inv;
    }
    Communication::communicate_field(z, grid.domain());
}

void PCG_RBLU::iterate(Fields &field, Grid &grid) {
    const auto  &cells    = grid.fluid_cells();
    const double n_local  = static_cast<double>(cells.size());
    const double n_global = Communication::reduce_sum(n_local);

    if (!_cells_partitioned) {
        const double dx = grid.dx(), dy = grid.dy();
        const double dii = 2.0 / (dx * dx) + 2.0 / (dy * dy);
        // IC(0) Schur complement diagonal:
        // s_ii = d_ii - 2*(1/dx²)²/d_ii - 2*(1/dy²)²/d_ii
        const double sii = dii
                         - 2.0 / (dx * dx * dx * dx) / dii
                         - 2.0 / (dy * dy * dy * dy) / dii;
        _d_inv = 1.0 / dii;
        _s_inv = 1.0 / sii;
        for (auto c : cells) {
            if ((c->i() + c->j()) % 2 == 0) _red_cells.push_back(c);
            else                              _black_cells.push_back(c);
        }
        _cells_partitioned = true;
    }

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _r(i, j) = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
    }

    apply_rblu(_r, _z, grid);

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _d(i, j) = _z(i, j);
    }

    double buf[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
    MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double rz = buf[0];
    double rr = buf[1];

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

        apply_rblu(_r, _z, grid);

        double buf2[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
        MPI_Allreduce(MPI_IN_PLACE, buf2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double rz_new = buf2[0];
        rr                  = buf2[1];
        const double beta   = rz_new / rz;

        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _d(i, j) = _z(i, j) + beta * _d(i, j);
        }
        rz = rz_new;
        ++_last_iter_count;
    }
}

double PCG_RBLU::calculate_residual(Fields &field, Grid &grid) {
    double rloc = 0.0;
    for (auto c : grid.fluid_cells()) {
        int i = c->i(), j = c->j();
        const double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += val * val;
    }
    return std::sqrt(rloc / static_cast<double>(grid.fluid_cells().size()));
}

// ---------------------------------------------------------------------------
// PCG_MG — Multigrid V-cycle preconditioner
// ---------------------------------------------------------------------------

PCG_MG::PCG_MG(double tolerance, int max_iter, int nc, int nr)
    : _tolerance(tolerance), _max_iter(max_iter),
      _r(nc, nr, 0.0), _d(nc, nr, 0.0), _q(nc, nr, 0.0), _z(nc, nr, 0.0) {}

double PCG_MG::dot_local(const Matrix<double> &a, const Matrix<double> &b,
                          const std::vector<Cell *> &cells) {
    double s = 0.0;
    for (auto c : cells)
        s += a(c->i(), c->j()) * b(c->i(), c->j());
    return s;
}

void PCG_MG::axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                  const std::vector<Cell *> &cells) {
    for (auto c : cells)
        y(c->i(), c->j()) += alpha * x(c->i(), c->j());
}

void PCG_MG::setup_hierarchy(const Grid &grid) {
    _levels.clear();
    int nx = grid.size_x();
    int ny = grid.size_y();
    double dx = grid.dx();
    double dy = grid.dy();
    // Build levels until grid becomes too small or max 6 levels
    while (nx >= 4 && ny >= 4 && (int)_levels.size() < 6) {
        _levels.push_back({nx, ny, dx, dy,
                           Matrix<double>(nx + 2, ny + 2, 0.0),
                           Matrix<double>(nx + 2, ny + 2, 0.0),
                           Matrix<double>(nx + 2, ny + 2, 0.0)});
        nx /= 2; ny /= 2; dx *= 2.0; dy *= 2.0;
    }
    // Coarsest level (direct solve)
    _levels.push_back({nx, ny, dx, dy,
                       Matrix<double>(nx + 2, ny + 2, 0.0),
                       Matrix<double>(nx + 2, ny + 2, 0.0),
                       Matrix<double>(nx + 2, ny + 2, 0.0)});
}

// nu Red-Black Gauss-Seidel sweeps on _levels[lv].
// Uses level-specific dx,dy (not global Discretization parameters).
// Solves: A_lv * z = r  where A_lv is Laplace with spacing dx_lv, dy_lv.
// Sign convention: A_cg = -Δₕ (positive semidefinite) → update:
//   z(i,j) = [r(i,j) + Σ_neighbours z(nb)/(d²)] / d_ii
void PCG_MG::smooth(int lv, int sweeps) {
    auto &L = _levels[lv];
    const double dxx = L.dx * L.dx;
    const double dyy = L.dy * L.dy;
    const double dii = 2.0 / dxx + 2.0 / dyy;

    for (int s = 0; s < sweeps; ++s) {
        // Red pass: (i+j) even
        for (int i = 1; i <= L.nx; ++i)
            for (int j = 1; j <= L.ny; ++j) {
                if ((i + j) % 2 != 0) continue;
                const double nb_sum = L.z(i-1,j)/dxx + L.z(i+1,j)/dxx
                                    + L.z(i,j-1)/dyy + L.z(i,j+1)/dyy;
                L.z(i, j) = (L.r(i, j) + nb_sum) / dii;
            }
        // Black pass: (i+j) odd
        for (int i = 1; i <= L.nx; ++i)
            for (int j = 1; j <= L.ny; ++j) {
                if ((i + j) % 2 != 1) continue;
                const double nb_sum = L.z(i-1,j)/dxx + L.z(i+1,j)/dxx
                                    + L.z(i,j-1)/dyy + L.z(i,j+1)/dyy;
                L.z(i, j) = (L.r(i, j) + nb_sum) / dii;
            }
    }
}

// Full-Weighting restriction: residual from level fine to level coarse.
// First computes tmp = r_fine - A_fine * z_fine (the actual residual),
// then restricts tmp into coarse.r. Resets coarse.z = 0 for correction equation.
void PCG_MG::restrict_residual(int fine, int coarse) {
    auto &F = _levels[fine];
    auto &C = _levels[coarse];
    const double dxx = F.dx * F.dx;
    const double dyy = F.dy * F.dy;
    const double dii = 2.0 / dxx + 2.0 / dyy;

    // tmp = r - A_cg*z = r - (-Δₕz) = r + Δₕz
    // Δₕz(i,j) = (z(i-1,j)+z(i+1,j))/dx² + (z(i,j-1)+z(i,j+1))/dy² - dii*z(i,j)
    for (int i = 1; i <= F.nx; ++i)
        for (int j = 1; j <= F.ny; ++j) {
            const double lap = (F.z(i-1,j)+F.z(i+1,j))/dxx
                             + (F.z(i,j-1)+F.z(i,j+1))/dyy
                             - dii * F.z(i, j);
            F.tmp(i, j) = F.r(i, j) + lap;   // r - A_cg*z  (A_cg = -Δₕ)
        }

    // Full-Weighting: 1/4 center, 1/8 edges, 1/16 corners
    for (int I = 1; I <= C.nx; ++I)
        for (int J = 1; J <= C.ny; ++J) {
            const int i = 2 * I, j = 2 * J;
            C.r(I, J) = 0.25   * F.tmp(i,   j  )
                      + 0.125  *(F.tmp(i-1, j  ) + F.tmp(i+1, j  )
                               + F.tmp(i,   j-1) + F.tmp(i,   j+1))
                      + 0.0625 *(F.tmp(i-1, j-1) + F.tmp(i+1, j-1)
                               + F.tmp(i-1, j+1) + F.tmp(i+1, j+1));
        }

    // Reset correction for coarse grid solve
    for (int I = 0; I <= C.nx + 1; ++I)
        for (int J = 0; J <= C.ny + 1; ++J)
            C.z(I, J) = 0.0;
}

// Bilinear prolongation: correction from coarse to fine level (additive).
// Coincident: z_fine(2I,2J) += z_coarse(I,J)
// Edge midpoints: average of two surrounding coarse points
// Cell centers: average of four surrounding coarse points
void PCG_MG::prolongate_correction(int coarse, int fine) {
    auto &C = _levels[coarse];
    auto &F = _levels[fine];

    for (int I = 1; I <= C.nx; ++I)
        for (int J = 1; J <= C.ny; ++J) {
            const int i = 2 * I, j = 2 * J;
            // Coincident points
            F.z(i,   j  ) += C.z(I, J);
            // Edge midpoints (clamp to fine grid interior)
            if (i + 1 <= F.nx) F.z(i+1, j  ) += 0.5*(C.z(I,J) + C.z(I+1,J));
            if (j + 1 <= F.ny) F.z(i,   j+1) += 0.5*(C.z(I,J) + C.z(I,J+1));
            // Cell center
            if (i + 1 <= F.nx && j + 1 <= F.ny)
                F.z(i+1, j+1) += 0.25*(C.z(I,J) + C.z(I+1,J)
                                       + C.z(I,J+1) + C.z(I+1,J+1));
        }
}

// Direct solve on coarsest level: SOR until residual < 1e-12 or 5000 iterations.
void PCG_MG::direct_solve(int lv) {
    auto &L = _levels[lv];
    const double dxx  = L.dx * L.dx;
    const double dyy  = L.dy * L.dy;
    const double dii  = 2.0 / dxx + 2.0 / dyy;
    const double coeff = 1.0 / dii;
    for (int iter = 0; iter < 5000; ++iter) {
        double res = 0.0;
        for (int i = 1; i <= L.nx; ++i)
            for (int j = 1; j <= L.ny; ++j) {
                const double nb_sum = L.z(i-1,j)/dxx + L.z(i+1,j)/dxx
                                    + L.z(i,j-1)/dyy + L.z(i,j+1)/dyy;
                const double z_new = (L.r(i, j) + nb_sum) * coeff;
                res += (z_new - L.z(i, j)) * (z_new - L.z(i, j));
                L.z(i, j) = z_new;
            }
        if (std::sqrt(res) < 1e-12) break;
    }
}

// V-Cycle: recursive multigrid solver.
// Pre-smooth, restrict, recurse, prolongate, post-smooth.
void PCG_MG::vcycle(int lv) {
    if (lv == (int)_levels.size() - 1) {
        direct_solve(lv);
        return;
    }
    smooth(lv, 2);                      // 2 pre-smoothing sweeps
    restrict_residual(lv, lv + 1);      // full-weighting restriction
    vcycle(lv + 1);                     // recursive call
    prolongate_correction(lv + 1, lv);  // bilinear prolongation
    smooth(lv, 2);                      // 2 post-smoothing sweeps
}

void PCG_MG::iterate(Fields &field, Grid &grid) {
    const auto  &cells    = grid.fluid_cells();
    const double n_local  = static_cast<double>(cells.size());
    const double n_global = Communication::reduce_sum(n_local);

    if (!_setup_done) {
        setup_hierarchy(grid);
        _setup_done = true;
    }

    // Warm start: r = Δₕp - RS
    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _r(i, j) = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
    }

    // Apply V-cycle as preconditioner: copy r into levels[0].r, run vcycle, extract z
    auto &L0 = _levels[0];
    for (int i = 1; i <= L0.nx; ++i)
        for (int j = 1; j <= L0.ny; ++j) {
            L0.r(i, j) = _r(i, j);
            L0.z(i, j) = 0.0;
        }
    vcycle(0);
    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _z(i, j) = L0.z(i, j);
    }

    for (auto c : cells) {
        int i = c->i(), j = c->j();
        _d(i, j) = _z(i, j);
    }

    double buf[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
    MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double rz = buf[0];
    double rr = buf[1];

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

        // Apply V-cycle to new residual
        for (int i = 1; i <= L0.nx; ++i)
            for (int j = 1; j <= L0.ny; ++j) {
                L0.r(i, j) = _r(i, j);
                L0.z(i, j) = 0.0;
            }
        vcycle(0);
        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _z(i, j) = L0.z(i, j);
        }

        double buf2[2] = { dot_local(_r, _z, cells), dot_local(_r, _r, cells) };
        MPI_Allreduce(MPI_IN_PLACE, buf2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double rz_new = buf2[0];
        rr                  = buf2[1];
        const double beta   = rz_new / rz;

        for (auto c : cells) {
            int i = c->i(), j = c->j();
            _d(i, j) = _z(i, j) + beta * _d(i, j);
        }
        rz = rz_new;
        ++_last_iter_count;
    }
}

double PCG_MG::calculate_residual(Fields &field, Grid &grid) {
    double rloc = 0.0;
    for (auto c : grid.fluid_cells()) {
        int i = c->i(), j = c->j();
        const double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += val * val;
    }
    return std::sqrt(rloc / static_cast<double>(grid.fluid_cells().size()));
}

// =============================================================================
// Eigen_CG — Eigen::ConjugateGradient reference solver (serial, 1×1 only)
//
// Purpose: identical algorithm to CG_Solver, but delegated to Eigen's
// optimised sparse CG so we can benchmark our hand-written implementation
// against an established reference library.
// Restriction: Eigen has no MPI support → only valid for iproc=jproc=1.
// =============================================================================
#ifdef USE_EIGEN

// ---------------------------------------------------------------------------
// Constructor
// Same signature as CG_Solver(tolerance, max_iter, nc, nr).
// nc/nr are the padded p_matrix dimensions; Eigen_CG does not pre-allocate
// auxiliary matrices from them — the system size N is determined from
// fluid_cells() in setup(), exactly as CG_Solver's loops skip non-fluid cells.
// ---------------------------------------------------------------------------
Eigen_CG::Eigen_CG(double tolerance, int max_iter, int /*nc*/, int /*nr*/)
    : _tolerance(tolerance), _max_iter(max_iter) {}

// ---------------------------------------------------------------------------
// setup() — called once on the first iterate(), never again.
//
// CG_Solver has no setup phase: it applies the Laplacian stencil on-the-fly
// via Discretization::laplacian() inside the iteration loop.  Eigen's CG
// requires an explicit matrix object, so we assemble A once here and reuse
// it every timestep.  The stencil coefficients are identical to those inside
// Discretization::laplacian() — only the representation differs.
// ---------------------------------------------------------------------------
void Eigen_CG::setup(Grid &grid) {
    const double dx2    = grid.dx() * grid.dx();
    const double dy2    = grid.dy() * grid.dy();
    _dx2 = dx2;
    _dy2 = dy2;
    // Diagonal entry of A_cg = -Δₕ.  Same value that Discretization::laplacian
    // would subtract from the centre cell when computing A_cg * p.
    const double center = 2.0 / dx2 + 2.0 / dy2;

    const auto &cells = grid.fluid_cells();
    _N = static_cast<int>(cells.size());  // size of the linear system

    // --- Index map: grid (i,j) → row index in A ---
    // CG_Solver works directly with Matrix<double>(i,j) and iterates over
    // fluid_cells() — no explicit numbering needed.  Eigen needs a dense
    // vector of length N, so we must assign each fluid cell a unique index.
    //
    // _idx is a flat array sized (imax+2) × (jmax+2), matching the padded
    // p_matrix layout.  _stride = jmax+2 so _idx[i*_stride+j] maps (i,j)
    // to a row — same formula as Matrix<double>'s internal indexing, so no
    // magic numbers and no hash collisions.  Non-fluid entries stay at -1.
    _stride = grid.domain().domain_jmax + 2;
    _idx.assign(static_cast<size_t>((grid.domain().domain_imax + 2) * _stride), -1);
    int row = 0;
    for (auto c : cells)
        _idx[static_cast<size_t>(c->i() * _stride + c->j())] = row++;

    // --- Assemble sparse matrix A in triplet (COO) format ---
    // Each fluid cell contributes up to 5 entries: the diagonal and up to 4
    // off-diagonals.  This is the explicit form of q = A_cg * d = -Δₕd that
    // CG_Solver computes implicitly via Discretization::laplacian(_d, i, j).
    //
    // Non-fluid neighbours (obstacle/ghost cells) have _idx = -1 and are
    // skipped — equivalent to CG_Solver reading p=0 from those cells because
    // BCs hold them at zero, so they contribute nothing to the stencil sum.
    _A.resize(_N, _N);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(5 * _N);  // at most 5 non-zeros per row

    for (auto c : cells) {
        int i = c->i();
        int j = c->j();
        int r = _idx[static_cast<size_t>(i * _stride + j)];

        // Diagonal: centre coefficient of -Δₕ (positive → A_cg is SPD → CG converges)
        triplets.emplace_back(r, r, center);

        // Off-diagonals: -1/dx² or -1/dy² for each fluid neighbour.
        // The sign is negative because A_cg = -Δₕ flips the off-diagonal sign
        // relative to the raw finite-difference stencil.
        auto add_nb = [&](int ni, int nj, double coeff) {
            int col = _idx[static_cast<size_t>(ni * _stride + nj)];
            if (col >= 0)  // skip if neighbour is obstacle or ghost cell
                triplets.emplace_back(r, col, -coeff);
        };
        add_nb(i - 1, j,     1.0 / dx2);  // left
        add_nb(i + 1, j,     1.0 / dx2);  // right
        add_nb(i,     j - 1, 1.0 / dy2);  // bottom
        add_nb(i,     j + 1, 1.0 / dy2);  // top
    }

    // Convert triplets → compressed sparse column format (CSC) used by Eigen.
    // Duplicate entries are summed automatically (none here, but good practice).
    _A.setFromTriplets(triplets.begin(), triplets.end());

    // compute() analyses the sparsity pattern and prepares the preconditioner.
    // For IdentityPreconditioner this is essentially a no-op, but it must be
    // called before solveWithGuess().  Keeping _solver as a member means this
    // call happens only once — not every timestep, unlike if solver were local.
    _solver.setMaxIterations(_max_iter);
    _solver.compute(_A);
    _setup_done = true;
}

// ---------------------------------------------------------------------------
// iterate() — called every timestep, mirrors CG_Solver::iterate() step by step.
// ---------------------------------------------------------------------------
void Eigen_CG::iterate(Fields &field, Grid &grid) {
    if (!_setup_done) setup(grid);  // lazy init on first call

    const auto &cells = grid.fluid_cells();

    // --- Build RHS b and initial guess x ---
    // CG_Solver works in-place: it computes r = Δₕp - RS directly into _r.
    // Eigen needs flat VectorXd.  The system being solved is A_cg · p = b
    // with A_cg = -Δₕ (positive definite).  The residual identity gives:
    //   r = b - A_cg · p  =  b + Δₕp
    // CG_Solver's initial residual is r = Δₕp - RS, so b = -RS.
    // Using +RS would solve A_cg · p = +RS → solution is -p_true (sign-flipped).
    // That bug is hidden by warm-start (0 iterations → p unchanged) but would
    // produce wrong results on cold start.
    //
    // x(k) = p(i,j): warm start — reuse pressure from previous timestep.
    // With the correct b = -RS: r₀ = -RS + Δₕp_prev ≈ 0 (since Δₕp_prev ≈ RS
    // from the previous solve), so Eigen exits in 0 iterations, identical to
    // CG_Solver which also checks sqrt(rr/N) < eps before the first iteration.
    // Build RHS b and initial guess x.
    //
    // CG_Solver computes r₀ = Δₕ(field.p_matrix(), i, j) - RS, which reads ALL
    // neighbours including ghost/boundary cells whose p is set by BCs (e.g. Neumann:
    // p_ghost = p_adjacent).  Eigen's internal r₀ = b - A*x skips non-fluid cells
    // because they are absent from A.  Without correction these residuals differ and
    // warm-start fails despite equivalent tolerances.
    //
    // Fix: for every non-fluid neighbour of a fluid cell add p(nb)/h² to b.
    // This makes Eigen's r₀ = b - A*x identical to CG_Solver's laplacian(p) - RS:
    //   b(k) - (A*x)(k) = [-RS + Σ_nonfluid p(nb)/h²] - [center*p - Σ_fluid p(nb)/h²]
    //                    = -RS + Σ_ALL p(nb)/h² - center*p
    //                    = Δₕ(p, i, j) - RS  ✓
    //
    // b must be rebuilt each timestep because ghost/BC pressures change.
    Eigen::VectorXd b(_N), x(_N);
    for (auto c : cells) {
        int i = c->i(), j = c->j();
        int k = _idx[static_cast<size_t>(i * _stride + j)];

        double rhs = -field.rs(i, j);

        // Add Neumann/obstacle-cell contributions that Discretization::laplacian()
        // picks up naturally but Eigen's matrix skips.
        auto add_bc = [&](int ni, int nj, double inv_h2) {
            if (_idx[static_cast<size_t>(ni * _stride + nj)] < 0)
                rhs += field.p(ni, nj) * inv_h2;
        };
        const double idx2 = 1.0 / _dx2;
        const double idy2 = 1.0 / _dy2;
        add_bc(i - 1, j,     idx2);
        add_bc(i + 1, j,     idx2);
        add_bc(i,     j - 1, idy2);
        add_bc(i,     j + 1, idy2);

        b(k) = rhs;
        x(k) = field.p(i, j);   // warm start
    }

    // --- Tolerance conversion ---
    // CG_Solver stops when sqrt(rr / N_global) < _tolerance  (absolute L2-norm).
    // Eigen stops when ‖r‖₂ < tol * ‖b‖₂                   (relative L2-norm).
    //
    // Mathematically equivalent conversion:
    //   CG_Solver: ‖r‖₂ / sqrt(N) < eps  →  ‖r‖₂ < eps * sqrt(N)
    //   Eigen:     ‖r‖₂ < tol * ‖b‖₂
    //   → tol = eps * sqrt(N) / ‖b‖₂
    //
    // This is the best possible match given Eigen's API.  A ±1 iteration
    // difference vs CG_Solver is possible because Eigen places its convergence
    // check internally (after the update) while CG_Solver checks at the top of
    // the loop (before the update).  This is a known limitation when comparing
    // two CG implementations via a black-box tolerance interface.
    //
    // Edge case: ‖b‖₂ ≈ 0 means RS ≈ 0 (trivially converged system).
    // Pass _tolerance directly so Eigen exits in 0 iterations.
    const double b_norm = b.norm();
    _solver.setTolerance(b_norm > 1e-300
        ? _tolerance * std::sqrt(static_cast<double>(_N)) / b_norm
        : _tolerance);

    // solveWithGuess(b, x): plain CG — no SOR, no preconditioner.
    // The second argument x is the warm-start vector (field.p from the previous
    // timestep), so CG starts from x₀ = p_old instead of x₀ = 0.
    // This keeps the initial residual r₀ = b - A·x₀ tiny and lets CG exit
    // in 0–2 iterations once the flow is established.
    x = _solver.solveWithGuess(b, x);

    // Store iteration count for CSV logging (same as CG_Solver's _last_iter_count).
    _last_iter_count = static_cast<int>(_solver.iterations());

    // --- Write solution back to field.p() ---
    // CG_Solver updates field.p_matrix() in-place via axpy(), so no copy-back
    // is needed there.  Here we must copy x → field.p() explicitly because
    // Eigen operates on its own VectorXd, not on the Matrix<double> directly.
    // Only fluid cells are written; obstacle/ghost cells are left unchanged —
    // BCs are reapplied by Case.cpp after iterate() returns, same as CG_Solver.
    for (auto c : cells) {
        int k = _idx[static_cast<size_t>(c->i() * _stride + c->j())];
        field.p(c->i(), c->j()) = x(k);
    }
}

// ---------------------------------------------------------------------------
// calculate_residual() — byte-identical to CG_Solver::calculate_residual().
// Recomputes sqrt(rr/N) from the current field.p() via Discretization::laplacian().
// Case.cpp calls this after iterate() and applies BCs first, so field.p() is
// already the post-solve, post-BC pressure when this runs.
// ---------------------------------------------------------------------------
double Eigen_CG::calculate_residual(Fields &field, Grid &grid) {
    double rloc = 0.0;
    for (auto c : grid.fluid_cells()) {
        int i = c->i(), j = c->j();
        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += val * val;
    }
    return std::sqrt(rloc / static_cast<double>(grid.fluid_cells().size()));
}

#endif // USE_EIGEN
