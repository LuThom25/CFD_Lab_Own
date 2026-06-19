#pragma once

#include <utility>

#include "Boundary.hpp"
#include "Fields.hpp"
#include "Grid.hpp"

/**
 * @brief Abstract class for pressure Poisson equation solver
 *
 */
class PressureSolver {
  public:
    PressureSolver() = default;
    virtual ~PressureSolver() = default;

    virtual void iterate(Fields &field, Grid &grid) = 0;
    virtual double calculate_residual(Fields &field, Grid &grid) = 0;
};

/**
 * @brief Standard SOR without null-space correction.
 *
 * Iciar's implementation. Plain SOR sweep followed by L2-norm residual
 * computation. No mean subtraction applied.
 */
class SOR_Standard : public PressureSolver {
  public:
    SOR_Standard() = default;

    /**
     * @brief Constructor of SOR_Standard solver
     *
     * @param[in] omega  SOR relaxation factor
     */
    SOR_Standard(double omega);

    virtual ~SOR_Standard() = default;

    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;

  private:
    double _omega;
};

/**
 * @brief Plain (unpreconditioned) Conjugate Gradient solver for the pressure
 *        Poisson equation.
 *
 * Design mirrors SOR_Standard: iterate() performs the full CG solve per
 * timestep, calculate_residual() recomputes the L2-norm from the current
 * pressure field — identical formula to SOR_Standard::calculate_residual().
 * No mean-pressure subtraction is applied; the Fredholm correction on RS
 * in Fields::calculate_rs() is sufficient, matching SOR behaviour.
 *
 * Auxiliary vectors (_r, _d, _q) are allocated once in the constructor
 * with the same dimensions as the pressure matrix (size_x+2, size_y+2).
 *
 * Algorithm executed inside iterate() each timestep:
 *   1.  r = Δₕp - RS          initial residual (warm start from previous timestep)
 *   2.  d = r                  first search direction
 *   3.  rr = reduce_sum(r·r)   global ‖r‖²                          [MPI]
 *   Loop until converged:
 *   a.  communicate(d)         halo exchange so stencil sees neighbours [MPI]
 *   b.  q = -Δₕd               matrix-vector product (A_cg = -Δₕ)
 *   c.  dq = reduce_sum(d·q)   global dot product for step size        [MPI]
 *   d.  α = rr / dq
 *   e.  p += α·d               solution update
 *   f.  r -= α·q               residual update
 *   g.  rr_new = reduce_sum(r·r)                                       [MPI]
 *   h.  β = rr_new / rr        conjugacy coefficient
 *   i.  d = r + β·d            new search direction
 *   j.  rr = rr_new
 */
class CG_Solver : public PressureSolver {
  public:
    CG_Solver() = default;

    /**
     * @brief Constructor of CG_Solver
     *
     * @param[in] tolerance  Convergence tolerance (eps from .dat file)
     * @param[in] max_iter   Maximum CG iterations (itermax from .dat file)
     * @param[in] nc         Number of matrix columns  (p_matrix().num_cols())
     * @param[in] nr         Number of matrix rows     (p_matrix().num_rows())
     */
    CG_Solver(double tolerance, int max_iter, int nc, int nr);

    virtual ~CG_Solver() = default;

    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;

    // Number of CG iterations performed in the last iterate() call.
    int last_iter_count() const { return _last_iter_count; }

  private:
    double _tolerance{0.0};
    int    _max_iter{0};
    int    _last_iter_count{0};

    Matrix<double> _r; // CG residual:      r = Δₕp - RS  (= b_cg - A_cg*p)
    Matrix<double> _d; // search direction: communicated before each matvec
    Matrix<double> _q; // matvec result:    q = A_cg*d = -Δₕd

    // Sums a(i,j)*b(i,j) over this rank's fluid cells only — no MPI.
    // Caller must follow up with Communication::reduce_sum() to get the
    // global dot product across all ranks (done in iterate() for alpha/beta).
    double dot_local(const Matrix<double> &a, const Matrix<double> &b,
                     const std::vector<Cell *> &cells);

    // AXPY: y(i,j) += alpha * x(i,j)  over fluid_cells()
    void axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
              const std::vector<Cell *> &cells);
};

/**
 * @brief Preconditioned Conjugate Gradient solver with Red-Black SSOR preconditioner.
 *
 * Design mirrors CG_Solver: iterate() performs the full PCG solve per
 * timestep, calculate_residual() recomputes the L2-norm from the current
 * pressure field — identical formula to CG_Solver::calculate_residual().
 * No mean-pressure subtraction is applied; the Fredholm correction on RS
 * in Fields::calculate_rs() is sufficient, matching SOR behaviour.
 *
 * Auxiliary vectors (_r, _d, _q, _z) are allocated once in the constructor
 * with the same dimensions as the pressure matrix (size_x+2, size_y+2).
 * Red/black cell lists are built once on the first iterate() call.
 *
 * The SSOR preconditioner M = (D+L)D⁻¹(D+Lᵀ) reduces κ(A) from O(N) to
 * O(√N), cutting iterations from O(√N) to O(N^{1/4}).  Red-Black ordering
 * ensures M = Mᵀ regardless of MPI domain decomposition — a requirement
 * for PCG correctness.  apply_ssor() executes 4 passes (forward: red→black,
 * backward: black→red), each followed by a halo exchange.
 *
 * Algorithm executed inside iterate() each timestep:
 *   1.  r  = Δₕp - RS         initial residual (warm start from previous timestep)
 *   2.  z  = apply_ssor(r)    preconditioned residual                [4 halos]
 *   3.  d  = z                 first search direction
 *   4.  ρ  = reduce_sum(r·z)  global M⁻¹-weighted dot product       [MPI]
 *   5.  rr = reduce_sum(r·r)  global ‖r‖² for stopping criterion    [MPI]
 *   Loop until √(rr/N) < eps:
 *   a.  communicate(d)         halo exchange so stencil sees neighbours [MPI]
 *   b.  q  = -Δₕd              matrix-vector product (A_cg = -Δₕ)
 *   c.  dq = reduce_sum(d·q)  global dot product for step size       [MPI]
 *   d.  α  = ρ / dq
 *   e.  p += α·d               solution update
 *   f.  r -= α·q               residual update
 *   g.  z  = apply_ssor(r)    preconditioned new residual            [4 halos]
 *   h.  ρ_new = reduce_sum(r·z)                                      [MPI]
 *   i.  β  = ρ_new / ρ        conjugacy coefficient
 *   j.  d  = z + β·d          new search direction  (z, not r — PCG)
 *   k.  ρ  = ρ_new
 *   l.  rr = reduce_sum(r·r)                                         [MPI]
 */
class PCG_SSOR : public PressureSolver {
  public:
    PCG_SSOR() = default;
    PCG_SSOR(double tolerance, int max_iter, int nc, int nr);
    virtual ~PCG_SSOR() = default;
    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;
    int last_iter_count() const { return _last_iter_count; }

  private:
    double _tolerance{0.0};
    int    _max_iter{0};
    int    _last_iter_count{0};

    Matrix<double> _r; // residual:             r = Δₕp - RS
    Matrix<double> _d; // search direction:     communicated before matvec
    Matrix<double> _q; // matvec result:        q = -Δₕd
    Matrix<double> _z; // preconditioned resid: z = M⁻¹r

    std::vector<Cell *> _red_cells;    // cells where (i+j)%2 == 0
    std::vector<Cell *> _black_cells;  // cells where (i+j)%2 == 1
    bool _cells_partitioned{false};

    // 4-pass Red-Black SSOR: solves M·z = r, stores result in z
    void apply_ssor(const Matrix<double> &r, Matrix<double> &z, const Grid &grid);

    double dot_local(const Matrix<double> &a, const Matrix<double> &b,
                     const std::vector<Cell *> &cells);
    void   axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                const std::vector<Cell *> &cells);
};

/**
 * @brief Serial (no-MPI) version of PCG_SSOR for pure-compute benchmarking.
 *
 * Algorithmically identical to PCG_SSOR but all communicate_field and
 * MPI_Allreduce calls are removed. Dot products use local values directly
 * (correct for 1-rank runs where all cells are local). Intended for
 * measuring the raw compute cost of PCG_SSOR without any MPI overhead.
 *
 * Use solver key: PCG_SSOR_SERIAL
 */
class PCG_SSOR_Serial : public PressureSolver {
  public:
    PCG_SSOR_Serial() = default;
    PCG_SSOR_Serial(double tolerance, int max_iter, int nc, int nr);
    virtual ~PCG_SSOR_Serial() = default;
    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;
    int last_iter_count() const { return _last_iter_count; }

  private:
    double _tolerance{0.0};
    int    _max_iter{0};
    int    _last_iter_count{0};

    Matrix<double> _r, _d, _q, _z;
    std::vector<Cell *> _red_cells, _black_cells;
    bool _cells_partitioned{false};

    void   apply_ssor(const Matrix<double> &r, Matrix<double> &z, const Grid &grid);
    double dot_local(const Matrix<double> &a, const Matrix<double> &b,
                     const std::vector<Cell *> &cells);
    void   axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                const std::vector<Cell *> &cells);
};
