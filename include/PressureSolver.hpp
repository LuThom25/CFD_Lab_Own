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
