#include <mpi.h>
#include <iostream>
#include <string>

#include "Case.hpp"

int main(int argn, char **args) {
    // Initialize MPI: assign rank numbers and define global communicator
    MPI_Init(&argn, &args); 
    
    // Get the number of processes in the communicator
    // -> MPI_COMM_WORLD is the global communicator comprising all processes 
    int size; // total number of processes, each of which owns one domain
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    
    // Get the assigned rank number of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (argn > 1) {
        std::string file_name{args[1]};
        // Need to pass to the constructor: size and my_rank, to know which 
        // subdomain it owns and who its neighbohrs are
        Case problem(file_name, argn, args, size, rank);
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
    // Wait for all processes to exit
    MPI_Finalize();
}
