#pragma once

#include <map>
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
    void applyVelocity(Fields &field);

    /**
     * @brief Method to patch the velocity boundary conditions to the given field for a specific cell.
     *
     * @param[in] Field to be applied
     * @param[in] Cell to which the boundary condition is applied
     */
    virtual void applyVelocityToCell(Fields &field, Cell *cell) = 0;

    /**
     * @brief Method to patch the pressure boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    void applyPressure(Fields &field);

    /**
     * @brief Method to patch the pressure boundary conditions to the given field for a specific cell.
     *
     * @param[in] Field to be applied
     * @param[in] Cell to which the boundary condition is applied
     */
    virtual void applyPressureToCell(Fields &field, Cell *cell) = 0;

    /**
     * @brief Method to patch the flux (F & G) boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    void applyFlux(Fields &field);

    /**
     * @brief Method to patch the flux (F & G) boundary conditions to the given field for a specific cell.
     *
     * @param[in] Field to be applied
     * @param[in] Cell to which the boundary condition is applied
     */
    virtual void applyFluxToCell([[maybe_unused]] Fields &field, [[maybe_unused]] Cell *cell) {};

    /**
     * @brief Method to patch the temperature boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyTemperature(Fields &field);

    /**
     * @brief Virtual destructor for the boundary class.
     */
    virtual ~Boundary() = default;

  protected:
    /**
     * @brief Constructor for the boundary class.
     */
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
    ~FixedWallBoundary() = default;

    void applyVelocityToCell(Fields &field, Cell *cell) override;

    void applyPressureToCell(Fields &field, Cell *cell) override;

    void applyTemperature(Fields &field) override;

    void applyFluxToCell(Fields &field, Cell *cell) override;

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
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity, std::map<int, double> wall_temperature);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                       std::map<int, double> wall_temperature);
    ~MovingWallBoundary() = default;

    void applyVelocityToCell(Fields &field, Cell *cell) override;

    void applyPressureToCell(Fields &field, Cell *cell) override;

    void applyFluxToCell(Fields &field, Cell *cell) override;

    void applyTemperature(Fields &field) override;

  private:
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};

/**
 * @brief Inflow boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, Neumann for pressure
 */
class InFlowBoundary : public Boundary {
  public:
    InFlowBoundary(std::vector<Cell *> cells, double u_in, double v_in);
    ~InFlowBoundary() = default;

    void applyVelocityToCell(Fields &field, Cell *cell) override;

    void applyPressureToCell(Fields &field, Cell *cell) override;

    void applyFluxToCell(Fields &field, Cell *cell) override;

  private:
    double _u_in;
    double _v_in;
};

/**
 * @brief Outflow boundary condition for the outer boundaries of the domain.
 * Neumann for velocities, Dirichlet(0) for pressure
 */
class OutFlowBoundary : public Boundary {
  public:
    OutFlowBoundary(std::vector<Cell *> cells);
    ~OutFlowBoundary() = default;

    void applyVelocityToCell(Fields &field, Cell *cell) override;

    void applyPressureToCell(Fields &field, Cell *cell) override;

    void applyFluxToCell(Fields &field, Cell *cell) override;
};
