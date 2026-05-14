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

    virtual void applyVelocityTop(Fields &field, int i, int j) = 0;
    virtual void applyVelocityBottom(Fields &field, int i, int j) = 0;
    virtual void applyVelocityLeft(Fields &field, int i, int j) = 0;
    virtual void applyVelocityRight(Fields &field, int i, int j) = 0;

    /**
     * @brief Method to patch the pressure boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    void applyPressure(Fields &field);

    virtual void applyPressureTop(Fields &field, int i, int j) = 0;
    virtual void applyPressureBottom(Fields &field, int i, int j) = 0;
    virtual void applyPressureLeft(Fields &field, int i, int j) = 0;
    virtual void applyPressureRight(Fields &field, int i, int j) = 0;

    /**
     * @brief Method to patch the flux (F & G) boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    void applyFlux(Fields &field);

    virtual void applyFluxTop(Fields &, int, int) {};
    virtual void applyFluxBottom(Fields &, int, int) {};
    virtual void applyFluxLeft(Fields &, int, int) {};
    virtual void applyFluxRight(Fields &, int, int) {};

    /**
     * @brief Method to patch the temperature boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyTemperature(Fields &field);

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
    ~FixedWallBoundary() = default;

    void applyVelocityTop(Fields &field, int i, int j) override;
    void applyVelocityBottom(Fields &field, int i, int j) override;
    void applyVelocityLeft(Fields &field, int i, int j) override;
    void applyVelocityRight(Fields &field, int i, int j) override;

    void applyPressureTop(Fields &field, int i, int j) override;
    void applyPressureBottom(Fields &field, int i, int j) override;
    void applyPressureLeft(Fields &field, int i, int j) override;
    void applyPressureRight(Fields &field, int i, int j) override;

    void applyFluxTop(Fields &field, int i, int j) override;
    void applyFluxBottom(Fields &field, int i, int j) override;
    void applyFluxLeft(Fields &field, int i, int j) override;
    void applyFluxRight(Fields &field, int i, int j) override;

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

    void applyVelocityTop(Fields &field, int i, int j) override;
    void applyVelocityBottom(Fields &field, int i, int j) override;
    void applyVelocityLeft(Fields &field, int i, int j) override;
    void applyVelocityRight(Fields &field, int i, int j) override;

    void applyPressureTop(Fields &field, int i, int j) override;
    void applyPressureBottom(Fields &field, int i, int j) override;
    void applyPressureLeft(Fields &field, int i, int j) override;
    void applyPressureRight(Fields &field, int i, int j) override;

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

    void applyVelocityTop(Fields &field, int i, int j) override;
    void applyVelocityBottom(Fields &field, int i, int j) override;
    void applyVelocityLeft(Fields &field, int i, int j) override;
    void applyVelocityRight(Fields &field, int i, int j) override;

    void applyPressureTop(Fields &field, int i, int j) override;
    void applyPressureBottom(Fields &field, int i, int j) override;
    void applyPressureLeft(Fields &field, int i, int j) override;
    void applyPressureRight(Fields &field, int i, int j) override;

  private:
    double _u_in;
    double _v_in;
};
