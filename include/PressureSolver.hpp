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
