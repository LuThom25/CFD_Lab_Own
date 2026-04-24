#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace filesystem = std::filesystem;

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

#include "Case.hpp"
#include "Enums.hpp"

Case::Case(std::string file_name, int argn, char **args) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu{};      /* viscosity   */
    double UI{};      /* velocity x-direction */
    double VI{};      /* velocity y-direction */
    double PI{};      /* pressure */
    double GX{};      /* gravitation x-direction */
    double GY{};      /* gravitation y-direction */
    double xlength{}; /* length of the domain x-dir.*/
    double ylength{}; /* length of the domain y-dir.*/
    double dt{};      /* time step */
    int imax{};       /* number of cells x-direction*/
    int jmax{};       /* number of cells y-direction*/
    double gamma{};   /* uppwind differencing factor*/
    double omg{};     /* relaxation factor */
    double tau{};     /* safety factor for time step*/
    int itermax{};    /* max. number of iterations for pressure per time step */
    double eps{};     /* accuracy bound for pressure*/

    // R1: fail fast if the input file cannot be opened.
    if (!file.is_open()) {
        std::cerr << "Error: could not open input file '" << file_name << "'.\n";
        std::exit(1);
    }

    {
        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            // R1: guard against empty token (e.g. trailing whitespace at EOF).
            if (var.empty()) continue;
            if (var[0] == '#') { /* ignore comment line */
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
            }
        }
    }
    file.close();

    // R1: validate that all required physical parameters were actually set and
    //     are physically meaningful — catch missing/misspelled keys early.
    if (nu <= 0.0)        { std::cerr << "Error: nu must be > 0 (got "      << nu      << ").\n"; std::exit(1); }
    if (imax <= 0)        { std::cerr << "Error: imax must be > 0 (got "    << imax    << ").\n"; std::exit(1); }
    if (jmax <= 0)        { std::cerr << "Error: jmax must be > 0 (got "    << jmax    << ").\n"; std::exit(1); }
    if (xlength <= 0.0)   { std::cerr << "Error: xlength must be > 0 (got " << xlength << ").\n"; std::exit(1); }
    if (ylength <= 0.0)   { std::cerr << "Error: ylength must be > 0 (got " << ylength << ").\n"; std::exit(1); }
    if (_t_end  <= 0.0)   { std::cerr << "Error: t_end must be > 0 (got "   << _t_end  << ").\n"; std::exit(1); }
    if (eps     <= 0.0)   { std::cerr << "Error: eps must be > 0 (got "     << eps     << ").\n"; std::exit(1); }
    if (itermax <= 0)     { std::cerr << "Error: itermax must be > 0 (got " << itermax << ").\n"; std::exit(1); }

    std::map<int, double> wall_vel;
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / static_cast<double>(imax);
    domain.dy = ylength / static_cast<double>(jmax);
    domain.domain_imax = imax;
    domain.domain_jmax = jmax;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;
    _nu  = nu;
    _omg = omg;

    // Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
    }
    if (not _grid.fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
    }
}

void Case::set_file_names(std::string file_name) {
    // CQ2: use std::filesystem::path for portable, readable path decomposition
    //      instead of the previous manual character-by-character reverse loop.
    filesystem::path fp(file_name);

    // Stem = filename without extension, e.g. "LidDrivenCavity"
    _case_name = fp.stem().string();

    // Parent directory + trailing separator, e.g. "../../example_cases/LidDrivenCavity/"
    // Used to resolve a relative geometry file path supplied in the .dat file.
    _prefix = fp.parent_path().string();
    if (!_prefix.empty()) _prefix += '/';

    // Output directory: <parent>/<case_name>_Output
    filesystem::path output_dir = fp.parent_path() / (_case_name + "_Output");
    _dict_name = output_dir.string();

    // Prepend parent directory to geometry file path if one was specified.
    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    try {
        filesystem::create_directory(output_dir);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply velocity boundary conditions for all the boundaries in _boundaries container
 *   using applyVelocity() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Apply Flux boundary conditions using applyFlux()
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver
 * - Update pressure boundary conditions after each iteration of the SOR solver
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - calculate the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface and/or common functions. You need to define functions with individual functionalities in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate() {
    double t             = 0.0;
    double dt            = _field.dt();
    int    timestep      = 0;
    int    vtk_count     = 0;
    double output_counter = 0.0;

    // Bookkeeping for console output and the machine-readable SUMMARY line
    int    last_sor_iter  = 0;
    double last_sor_res   = 0.0;
    long   total_sor_iter = 0;
    int    max_sor_iter   = 0;
    double total_res_sum  = 0.0;   // sum of achieved residuals for avg_res
    double total_dt_sum   = 0.0;
    bool   diverged       = false;

    // ── Start-up header ────────────────────────────────────────────────────────
    // Re = U*L/nu with U=1, L=1 (lid-driven cavity scaling)
    double Re = (_nu > 0.0) ? (1.0 / _nu) : std::numeric_limits<double>::infinity();
    std::cout
        << "\n============================================================\n"
        << " Fluidchen CFD Solver  —  " << _case_name << "\n"
        << "============================================================\n"
        << std::fixed << std::setprecision(4)
        << "  Grid    : " << _grid.size_x() << " x " << _grid.size_y()
        << "  (dx=" << _grid.dx() << "  dy=" << _grid.dy() << ")\n"
        << "  nu      : " << _nu << "   Re ~ " << std::setprecision(0) << Re << "\n"
        << std::setprecision(2)
        << "  t_end   : " << _t_end << "   output every " << _output_freq << " time units\n";

    if (_field.tau() > 0.0)
        std::cout << "  dt      : adaptive (tau=" << _field.tau()
                  << ")   initial dt = " << std::setprecision(6) << dt << "\n";
    else
        std::cout << "  dt      : FIXED = " << std::setprecision(6) << dt
                  << "   (adaptive disabled — tau <= 0)\n";

    std::cout
        << "  SOR     : omega=" << std::setprecision(2) << _omg
        << "   itermax=" << _max_iter
        << "   eps=" << std::scientific << std::setprecision(2) << _tolerance << "\n"
        << "  Output  : " << _dict_name << "/\n"
        << "------------------------------------------------------------\n"
        << std::flush;

    // ── Initial state ──────────────────────────────────────────────────────────
    output_vtk(timestep);
    vtk_count++;

    // ── Main loop (Chorin Projection / fractional-step method) ────────────────
    while (t < _t_end) {

        // Step 1: Velocity BCs (ghost cells, moving lid)
        for (auto &boundary : _boundaries)
            boundary->applyVelocity(_field);

        // Step 2: Intermediate fluxes F and G (Eq. 9 & 10)
        _field.calculate_fluxes(_grid);

        // Step 3: Flux BCs at walls
        for (auto &boundary : _boundaries)
            boundary->applyFlux(_field);

        // Step 4: RHS of pressure Poisson equation (Eq. 11)
        _field.calculate_rs(_grid);

        // Step 5: SOR pressure solve — iterate until res < eps or itermax reached
        int    iter     = 0;
        double residual = std::numeric_limits<double>::max();
        while (iter < _max_iter && residual > _tolerance) {
            residual = _pressure_solver->solve(_field, _grid, _boundaries);
            for (auto &boundary : _boundaries)
                boundary->applyPressure(_field);
            ++iter;
        }
        last_sor_iter  = iter;
        last_sor_res   = residual;
        total_sor_iter += iter;
        total_res_sum  += residual;
        if (iter > max_sor_iter) max_sor_iter = iter;

        // Divergence check — NaN/Inf residual or runaway values signal instability
        if (std::isnan(residual) || std::isinf(residual) || residual > 1.0e8) {
            std::cout
                << "\n[DIVERGED] t=" << std::fixed << std::setprecision(4) << t
                << "  step=" << timestep
                << "  residual=" << std::scientific << std::setprecision(2) << residual << "\n"
                << "           Possible cause: dt too large (CFL violation).\n"
                << "           Try smaller dt, or enable adaptive stepping (tau > 0).\n";
            diverged = true;
            break;
        }

        // Step 6: Correct velocities using updated pressure (Eq. 7 & 8)
        _field.calculate_velocities(_grid);

        // Step 7: Compute adaptive dt for the next step (Eqs. 12 & 13)
        dt = _field.calculate_dt(_grid);
        total_dt_sum += dt;

        t += dt;
        ++timestep;
        output_counter += dt;

        // Step 8: VTK output + progress line at each output interval
        if (output_counter >= _output_freq) {
            output_vtk(timestep);
            vtk_count++;
            output_counter -= _output_freq;

            std::cout
                << "  [vtk=" << std::setw(3) << vtk_count
                << " | t=" << std::fixed << std::setprecision(3) << std::setw(8) << t
                << " | step=" << std::setw(6) << timestep
                << " | dt=" << std::scientific << std::setprecision(2) << dt
                << " | SOR: iter=" << std::setw(3) << last_sor_iter
                << "  res=" << std::setprecision(2) << last_sor_res
                << "]\n" << std::flush;
        }
    }

    // ── Final summary ─────────────────────────────────────────────────────────
    double avg_sor = (timestep > 0) ? static_cast<double>(total_sor_iter) / timestep : 0.0;
    double avg_res = (timestep > 0) ? total_res_sum / timestep : last_sor_res;
    double avg_dt  = (timestep > 0) ? total_dt_sum / timestep : dt;

    std::cout
        << "------------------------------------------------------------\n"
        << "  Done: t=" << std::fixed << std::setprecision(3) << t
        << "  steps=" << timestep << "  VTK files=" << vtk_count << "\n"
        << "  SOR : avg=" << std::fixed << std::setprecision(1) << avg_sor
        << "  max=" << max_sor_iter
        << "  avg_dt=" << std::scientific << std::setprecision(2) << avg_dt
        << "\n\n";

    // One-line machine-readable summary (parsed by run_studies.py)
    std::cout
        << "SUMMARY"
        << " t="       << std::fixed     << std::setprecision(3) << t
        << " steps="   << timestep
        << " vtk="     << vtk_count
        << " avg_sor=" << std::fixed      << std::setprecision(1) << avg_sor
        << " max_sor=" << max_sor_iter
        << " avg_res=" << std::scientific << std::setprecision(2) << avg_res
        << " avg_dt="  << std::scientific << std::setprecision(2) << avg_dt
        << " status="  << (diverged ? "DIVERGED" : "OK")
        << "\n";
}

void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().iminb * dx;
    double y = _grid.domain().jminb * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().iminb * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkSmartPointer<vtkDoubleArray> Pressure = vtkSmartPointer<vtkDoubleArray>::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array for cell data
    vtkSmartPointer<vtkDoubleArray> Velocity = vtkSmartPointer<vtkDoubleArray>::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print pressure, velocity and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
            vel[0] = (_field.u(i - 1, j) + _field.u(i, j)) * 0.5;
            vel[1] = (_field.v(i, j - 1) + _field.v(i, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Velocity Array for point data
    vtkSmartPointer<vtkDoubleArray> VelocityPoints = vtkSmartPointer<vtkDoubleArray>::New();
    VelocityPoints->SetName("velocity");
    VelocityPoints->SetNumberOfComponents(3);

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            VelocityPoints->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetCellData()->AddArray(Velocity);
    structuredGrid->GetPointData()->AddArray(VelocityPoints);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();

    // R3: verify the file was actually created — catches silent failures due to
    //     disk-full conditions or permission errors (vtkStructuredGridWriter
    //     returns void and does not throw on failure).
    if (!filesystem::exists(outputname)) {
        std::cerr << "Warning: VTK output file was not created: " << outputname << "\n"
                  << "         Check available disk space and write permissions.\n";
    }
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    domain.iminb = 0;
    domain.jminb = 0;
    domain.imaxb = imax_domain + 2;
    domain.jmaxb = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
}
