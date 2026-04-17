#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"

/**
 * @brief Abstract of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Method to patch the velocity boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyVelocity(Fields &field) = 0;

    /**
     * @brief Method to patch the pressure boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyPressure(Fields &field) = 0;

    /**
     * @brief Method to patch the flux (F & G) boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyFlux(Fields &field);

    virtual ~Boundary() = default;

  protected:
    Boundary(std::vector<Cell *> cells);
    std::vector<Cell *> _cells;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature);
    virtual ~FixedWallBoundary() = default;
    virtual void applyVelocity(Fields &field);
    virtual void applyPressure(Fields &field);

  private:
    std::map<int, double> _wall_temperature;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                       std::map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void applyVelocity(Fields &field);
    virtual void applyPressure(Fields &field);

  private:
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};
