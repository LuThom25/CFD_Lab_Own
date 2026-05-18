#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"

/**
 * @brief Class to hold and orchestrate the simulation flow.
 *
 */
class Case {
  public:
    /**
     * @brief Parallel constructor for the Case.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters. Creates output directory.
     *
     * @param[in] Input file name
     */
    Case(std::string file_name, int argn, char **args);

    /**
     * @brief Main function to simulate the flow until the end time.
     *
     * Calculates the fluxes
     * Calculates the right hand side
     * Solves pressure
     * Calculates velocities
     * Outputs the solution files
     */
    void simulate();

  private:
    /// Plain case name without paths
    std::string _case_name;
    /// Output directory name
    std::string _dict_name;
    /// Geometry file name
    std::string _geom_name{"NONE"};
    /// Relative input file path
    std::string _prefix;

    /// Simulation time
    double _t_end;
    /// Solution file outputting frequency
    double _output_freq;

    Fields _field;
    Grid _grid;
    Discretization _discretization;
    std::unique_ptr<PressureSolver> _pressure_solver;
    /// Collection of all boundaries for the simulated case
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    /// Solver convergence tolerance
    double _tolerance;

    /// Maximum number of iterations for the solver
    int _max_iter;

    /// Kinematic viscosity (stored for console output)
    double _nu{0.0};
    /// SOR relaxation factor (stored for console output)
    double _omg{0.0};
    /// Name of the active pressure solver (stored for console output)
    std::string _solver_name{"SOR"};

    // ── WS2: geometry and boundary parameters ─────────────────────────────────
    /// Inlet x-velocity read from "UIN" in .dat (used by Patrick's InFlowBoundary)
    double _UIN{0.0};
    /// Inlet y-velocity read from "VIN" in .dat (used by Patrick's InFlowBoundary)
    double _VIN{0.0};
    /// True when "energy_eq on" appears in .dat (used by Dani/Iciar)
    bool _energy_eq{false};
    /// Wall temperatures keyed by PGM cell ID (3–7); -1.0 = adiabatic.
    /// Read from wall_temp_3/4/5 entries in .dat. Used by Dani for temperature BCs.
    std::map<int, double> _wall_temperatures;

    // WS2 energy transport variables
    double TI{};    // initial temperature
    double alpha{}; // thermal conductivity
    double beta{};  // thermal expansion

    /**
     * @brief Creating file names from given input data file
     *
     * Extracts path of the case file and creates code-readable file names
     * for outputting directory and geometry file.
     *
     * @param[in] input data file
     */
    void set_file_names(std::string file_name);

    /**
     * @brief Solution file outputter
     *
     * Outputs the solution files in .vtk format. Ghost cells are excluded.
     * Pressure is cell variable while velocity is point variable while being
     * interpolated to the cell faces
     *
     * @param[in] Timestep of the solution
     * @param[in] Current rank of the executing process
     */
    void output_vtk(int t, int my_rank = 0);

    /**
     * @brief Fill out domain object
     *
     * Fills the Domain object with the geometrical information about the domain size. In case of running the code in
     * parallel, the information should correspond to the subdomain belonging to the executing MPI-process after
     * decomposition.
     *
     * @param[in] Reference to the domain object
     * @param[in] Number of cells in x-direction for this MPI rank
     * @param[in] Number of cells in y-direction for this MPI rank
     */
    void build_domain(Domain &domain, int imax_domain, int jmax_domain);
};
