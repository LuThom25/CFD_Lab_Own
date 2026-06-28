#pragma once

#include <utility>
#include <vector>

#ifdef USE_EIGEN
#include <Eigen/Sparse>
#endif

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

#ifdef USE_EIGEN
/**
 * @brief Eigen::ConjugateGradient wrapper — identical interface to CG_Solver.
 *
 * Delegates the linear solve to Eigen's built-in CG with IdentityPreconditioner
 * (no preconditioning), making it a 1-to-1 algorithmic reference for CG_Solver.
 * Serial only (Eigen has no MPI support); valid for 1×1 decomposition.
 *
 * Matrix assembly (setup): the 5-point Laplace operator is assembled once as an
 * Eigen::SparseMatrix<double> using the same stencil coefficients as
 * Discretization::laplacian().  Fluid cells are numbered 0..N-1 via _idx map;
 * obstacle cells are excluded from the system (same convention as CG_Solver
 * which skips non-fluid cells in all loops).
 *
 * iterate() mirrors CG_Solver::iterate() step by step:
 *   1.  Build RHS b from field.rs()           — same as CG_Solver initial r
 *   2.  Build initial guess x from field.p()  — warm start, same as CG_Solver
 *   3.  Set tolerance: solver.setTolerance(_tolerance / sqrt(N))
 *       so Eigen's ‖r‖₂ < tol*‖b‖₂ matches our sqrt(rr/N) < eps criterion
 *   4.  solver.solveWithGuess(b, x)           — Eigen CG loop
 *   5.  Write x back to field.p()
 *
 * calculate_residual(): identical formula to CG_Solver::calculate_residual(),
 * recomputed from the current field.p() via Discretization::laplacian().
 */
class Eigen_CG : public PressureSolver {
  public:
    // Same constructor signature as CG_Solver — nc/nr are p_matrix dimensions.
    Eigen_CG(double tolerance, int max_iter, int nc, int nr);
    ~Eigen_CG() = default;

    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;

    int last_iter_count() const { return _last_iter_count; }

  private:
    double _tolerance;
    int    _max_iter;
    int    _last_iter_count{0};
    bool   _setup_done{false};

    // Assembled once in setup(); grid is static so A never changes.
    Eigen::SparseMatrix<double> _A;

    // Solver kept as member so compute() (symbolic+numeric factorisation of the
    // preconditioner) is called only once in setup(), not every timestep.
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                             Eigen::Lower | Eigen::Upper,
                             Eigen::IdentityPreconditioner> _solver;

    // Flat index map: _idx[i * _stride + j] = row in _A.
    // _stride = number of j-indices per row (domain_jmax + 2).
    // Avoids magic numbers; O(1) lookup without hash overhead.
    std::vector<int> _idx;
    int _stride{0};

    // Total number of fluid cells = size of the linear system.
    int _N{0};

    // Grid spacing squared — needed every timestep to add BC contributions to b.
    double _dx2{0.0};
    double _dy2{0.0};

    void setup(Grid &grid);
};
#endif // USE_EIGEN

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

/**
 * @brief PCG with diagonal (Jacobi) preconditioner.
 *
 * M = D_A where d_ii = 2/dx² + 2/dy² (constant on uniform grid).
 * apply_jacobi: z(i,j) = r(i,j) / d_ii — purely local, no communicate_field.
 * Convergence: ~10-20% fewer iterations than plain CG.
 *
 * Use solver key: PCG_JACOBI
 */
class PCG_Jacobi : public PressureSolver {
  public:
    PCG_Jacobi() = default;
    PCG_Jacobi(double tolerance, int max_iter, int nc, int nr);
    virtual ~PCG_Jacobi() = default;
    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;
    int last_iter_count() const { return _last_iter_count; }

  private:
    double _tolerance{0.0};
    int    _max_iter{0};
    int    _last_iter_count{0};
    double _d_inv{0.0};    // 1/(2/dx²+2/dy²), computed once on first iterate()
    bool   _setup_done{false};

    Matrix<double> _r, _d, _q, _z;

    void   apply_jacobi(const Matrix<double> &r, Matrix<double> &z,
                        const std::vector<Cell *> &cells);
    double dot_local(const Matrix<double> &a, const Matrix<double> &b,
                     const std::vector<Cell *> &cells);
    void   axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                const std::vector<Cell *> &cells);
};

/**
 * @brief PCG with Incomplete Cholesky IC(0) preconditioner in Red-Black ordering.
 *
 * Red-Black ordering gives block structure [D_R, A_RB; A_BR, D_B].
 * Schur complement diagonal (IC(0) approximation):
 *   s_ii = d_ii - 2*(1/dx²)²/d_ii - 2*(1/dy²)²/d_ii
 * For dx=dy=h: s_ii = 3/h² (25% below d_ii → stronger preconditioning).
 * apply_rblu: 3 passes identical to apply_ssor but black pass uses s_inv.
 * Convergence: ~30-50% fewer iterations than PCG_SSOR.
 *
 * Use solver key: PCG_RBLU
 */
class PCG_RBLU : public PressureSolver {
  public:
    PCG_RBLU() = default;
    PCG_RBLU(double tolerance, int max_iter, int nc, int nr);
    virtual ~PCG_RBLU() = default;
    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;
    int last_iter_count() const { return _last_iter_count; }

  private:
    double _tolerance{0.0};
    int    _max_iter{0};
    int    _last_iter_count{0};
    double _d_inv{0.0};    // 1/d_ii = 1/(2/dx²+2/dy²)
    double _s_inv{0.0};    // 1/s_ii  (Schur complement diagonal)
    bool   _cells_partitioned{false};

    Matrix<double> _r, _d, _q, _z;
    std::vector<Cell *> _red_cells, _black_cells;

    void   apply_rblu(const Matrix<double> &r, Matrix<double> &z, const Grid &grid);
    double dot_local(const Matrix<double> &a, const Matrix<double> &b,
                     const std::vector<Cell *> &cells);
    void   axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                const std::vector<Cell *> &cells);
};

/**
 * @brief PCG with Multigrid V-cycle preconditioner.
 *
 * V-cycle: nu1=nu2=2 Red-Black GS pre/post sweeps, full-weighting restriction,
 * bilinear prolongation, direct solve (SOR) on coarsest grid.
 * Grid hierarchy: halve nx,ny until min(nx,ny) < 4 (max 6 levels).
 * Galerkin coarsening: coarse operator = standard Laplace with doubled spacing.
 * Convergence: h-independent, typically < 10 iterations for Poisson.
 *
 * Use solver key: PCG_MG
 */
class PCG_MG : public PressureSolver {
  public:
    PCG_MG() = default;
    PCG_MG(double tolerance, int max_iter, int nc, int nr);
    virtual ~PCG_MG() = default;
    void iterate(Fields &field, Grid &grid) override;
    double calculate_residual(Fields &field, Grid &grid) override;
    int last_iter_count() const { return _last_iter_count; }

  private:
    double _tolerance{0.0};
    int    _max_iter{0};
    int    _last_iter_count{0};
    bool   _setup_done{false};

    Matrix<double> _r, _d, _q, _z;    // CG vectors on finest level

    struct MG_Level {
        int nx, ny;        // interior grid points (without halos)
        double dx, dy;
        Matrix<double> r;  // right-hand side / residual  (nx+2, ny+2)
        Matrix<double> z;  // solution / correction        (nx+2, ny+2)
        Matrix<double> tmp; // workspace for residual      (nx+2, ny+2)
    };
    std::vector<MG_Level> _levels;    // _levels[0]=finest, _levels.back()=coarsest

    void setup_hierarchy(const Grid &grid);
    void vcycle(int level);
    void smooth(int level, int sweeps);
    void restrict_residual(int fine, int coarse);
    void prolongate_correction(int coarse, int fine);
    void direct_solve(int level);

    double dot_local(const Matrix<double> &a, const Matrix<double> &b,
                     const std::vector<Cell *> &cells);
    void   axpy(double alpha, const Matrix<double> &x, Matrix<double> &y,
                const std::vector<Cell *> &cells);
};
