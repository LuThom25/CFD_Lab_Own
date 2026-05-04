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

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) = 0;
};

/**
 * @brief SOR with Fredholm compatibility fix and zero-mean pressure projection.
 *
 * Thomas's implementation. Subtracts mean(RS) after each RHS assembly and
 * mean(p) after every sweep so the iterate stays in the subspace orthogonal
 * to the null space of the all-Neumann Poisson system, enabling convergence.
 */
class SOR_Mean_Correction : public PressureSolver {
  public:
    SOR_Mean_Correction() = default;

    /**
     * @brief Constructor of SOR_Mean_Correction solver
     *
     * @param[in] omega  SOR relaxation factor
     */
    SOR_Mean_Correction(double omega);

    virtual ~SOR_Mean_Correction() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega;
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

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega;
};

/**
 * @brief Red-Black Successive Over-Relaxation solver.
 *        Cells are updated in a checkerboard pattern (red first, then black),
 *        which decouples neighbouring updates within each colour pass and
 *        allows straightforward OpenMP parallelisation in future.
 */
class SOR_RB : public PressureSolver {
  public:
    SOR_RB() = default;

    /**
     * @brief Constructor of SOR_RB solver
     *
     * @param[in] omega  relaxation factor
     */
    SOR_RB(double omega);

    virtual ~SOR_RB() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid  to be used
     * @param[in] boundaries to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega;
};
