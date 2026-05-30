#include <iostream>
#include <string>

#include "Case.hpp"
#include "Communication.hpp"

int main(int argn, char **args) {
    Communication::init_parallel(argn, args); // Initialize MPI
    
    if (argn > 1) {
        std::string file_name{args[1]};
        // Need to pass to the constructor: size and my_rank, to know which 
        // subdomain it owns and who its neighbohrs are
        Case problem(file_name, argn, args, Communication::get_size(), Communication::get_rank());
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
    Communication::finalize(); // Finalize MPI
    
}
