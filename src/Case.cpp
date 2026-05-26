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
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

#include "Case.hpp"
#include "Enums.hpp"

Case::Case(std::string file_name, int /*argn*/, char ** /*args*/) {
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
    // WS2 inlet velocities: used by the InFlowBoundary
    double UIN{0.0};
    double VIN{0.0};
    double TI{0.0};    /* initial temperature */
    double alpha{0.0}; /* thermal diffusivity */
    double beta{0.0};  /* thermal expansion coefficient */

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
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

                /* WS2 parameters
                 *  Geometry file: read directly into _geom_name so set_file_names() can
                 *  prepend the directory prefix. If absent, _geom_name stays "NONE" and
                 *  the hardcoded lid-driven cavity fallback is used.
                 */
                if (var == "geo_file") file >> _geom_name; /* */
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                if (var == "energy_eq") {
                    std::string s;
                    file >> s;
                    _energy_eq = (s == "on");
                }
                if (var == "TI") {
                    file >> TI;
                }
                if (var == "alpha") {
                    file >> alpha;
                }
                if (var == "beta") {
                    file >> beta;
                }
                // Read wall temperatures
                if (var.size() > 10 && var.rfind("wall_temp_", 0) == 0) {
                    int wall_id = std::stoi(var.substr(10));
                    double temp;
                    file >> temp;
                    _wall_temperatures[wall_id] = temp;
                }
                // Both spellings appear across the provided .dat files; just consume the value.
                if (var == "num_walls" || var == "num_of_walls") {
                    int n;
                    file >> n;
                }
            }
        }
    }
    file.close();
    _UIN = UIN;
    _VIN = VIN;

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
    _field = Fields(nu, dt, tau, alpha, beta, GX, GY, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, TI);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR_Standard>(omg);
    _max_iter = itermax;
    _tolerance = eps;
    _nu = nu;
    _omg = omg;

    // Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(),
                                                                   LidDrivenCavity::wall_velocity, _wall_temperatures));
    }
    if (not _grid.fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells(), _wall_temperatures));
    }
    if (not _grid.inflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<InFlowBoundary>(_grid.inflow_cells(), _UIN, _VIN));
    }
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutFlowBoundary>(_grid.outflow_cells()));
    }

}

void Case::set_file_names(std::string file_name) {
    /* Use std::filesystem::path for portable, readable path decomposition
     * instead of the previous manual character-by-character reverse loop.
     */
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
void Case::simulate() { // Inialize variables
    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    int vtk_count = 0;
    double output_counter = 0.0;

    std::cout << "Starting simulation: " << _case_name << " t=" << std::fixed << std::setprecision(3) << t
              << std::endl;

    // Initial state
    if (_energy_eq) {
        for (auto &boundary : _boundaries)
            boundary->applyTemperature(_field);
    }
    output_vtk(timestep);
    vtk_count++;

    // Loop over time
    while (t < _t_end) {

        // Step 1: Velocity BCs (ghost cells, moving lid)
        for (auto &boundary : _boundaries)
            boundary->applyVelocity(_field);

        // Step 2: Compute new temperature values (dependent on u and v) and apply BCs
        if (_energy_eq) {
            _field.calculate_temperature(_grid);
            for (auto &boundary : _boundaries)
                boundary->applyTemperature(_field);
        }

        // Step 3: Compute intermediate fluxes F and G (tilda) and their BCs
        _field.calculate_fluxes(_grid);
        for (auto &boundary : _boundaries)
            boundary->applyFlux(_field);

        // Step 4: Calculate the RHS of pressure Poisson equation
        _field.calculate_rs(_grid);

        // Step 5: SOR pressure solve: iterate until res < eps or itermax reached
        // Initialize SOR stopping criteria
        int iter = 0;
        double residual = std::numeric_limits<double>::max();
        // Call SOR solve outputting the residual after solving the PPE
        while (iter < _max_iter && residual > _tolerance) {
            residual = _pressure_solver->solve(_field, _grid, _boundaries);
            // Apply BCs on the pressure field
            for (auto &boundary : _boundaries)
                boundary->applyPressure(_field);
            ++iter;
        }

        // Divergence check
        if (std::isnan(residual) || std::isinf(residual) || residual > 1.0e8) {
            ++timestep; // count this step so avg_sor = total_iter / total_steps is bounded by itermax
            std::cout << "\n[DIVERGED] t=" << std::fixed << std::setprecision(4) << t << "  step=" << timestep
                      << "  residual=" << std::scientific << std::setprecision(2) << residual << "\n";
            break;
        }

        // Step 6: Correct velocities using updated pressure
        _field.calculate_velocities(_grid);

        // Step 7: Compute adaptive dt for the next timestep
        dt = _field.calculate_dt(_grid);

        // Update time, timestep and counter
        t += dt;
        ++timestep;
        output_counter += dt;

        // Step 8: VTK output
        // only output when enough simulated time has passed to avoid having thousands of ouput vtk
        if (output_counter >= _output_freq) {
            output_vtk(timestep);
            vtk_count++;
            output_counter -= _output_freq;
        }
    }

    std::cout << "Finished simulation: " << _case_name << " t=" << std::fixed << std::setprecision(3) << t
              << std::endl;
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

    // Temperature array
    vtkSmartPointer<vtkDoubleArray> Temperature = vtkSmartPointer<vtkDoubleArray>::New();
    Temperature->SetName("temperature");
    Temperature->SetNumberOfComponents(1);

    // Temp Velocity
    std::array<double, 3> vel;
    vel[2] = 0; // Set z component to 0

    // Print pressure, velocity and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
            vel[0] = (_field.u(i - 1, j) + _field.u(i, j)) * 0.5;
            vel[1] = (_field.v(i, j - 1) + _field.v(i, j)) * 0.5;
            Velocity->InsertNextTuple(vel.data());
            double temp = _field.t(i, j);
            Temperature->InsertNextTuple(&temp);
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
            VelocityPoints->InsertNextTuple(vel.data());
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetCellData()->AddArray(Velocity);
    structuredGrid->GetPointData()->AddArray(VelocityPoints);
    structuredGrid->GetCellData()->AddArray(Temperature);

    /* Geometry / obstacle field
     *  Encode cell type as an integer per cell so ParaView can colour obstacles without
     *  needing a separate geometry file. Domain boundary ghost cells are omitted (loop
     *  starts at 1), matching the same range used for Pressure and Velocity above.
     *  0 = FLUID, 1 = FIXED_WALL, 2 = MOVING_WALL, 3 = INFLOW, 4 = OUTFLOW
     */
    vtkSmartPointer<vtkIntArray> Obstacle = vtkSmartPointer<vtkIntArray>::New();
    Obstacle->SetName("obstacle");
    Obstacle->SetNumberOfComponents(1);
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            int flag = 0;
            switch (_grid.cells()(i, j).type()) {
            case cell_type::FIXED_WALL:
                flag = 1;
                break;
            case cell_type::MOVING_WALL:
                flag = 2;
                break;
            case cell_type::INFLOW:
                flag = 3;
                break;
            case cell_type::OUTFLOW:
                flag = 4;
                break;
            default:
                flag = 0;
                break;
            }
            Obstacle->InsertNextValue(flag);
        }
    }
    structuredGrid->GetCellData()->AddArray(Obstacle);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    domain.iminb = 0;
    domain.jminb = 0;
    domain.imaxb = imax_domain + 2;
    domain.jmaxb = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
}
