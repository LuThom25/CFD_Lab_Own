// Communication.cpp — Placeholder for Worksheet 2 (MPI parallelisation)
//
// This file is intentionally empty in Worksheet 1, which runs on a single
// process only.  In a parallel extension (WS2), this class would implement
// halo-exchange routines that synchronise ghost-cell values across MPI ranks
// after each SOR sweep and after computing fluxes/velocities.
//
// Typical methods to implement here:
//   void Communication::communicate(Fields &field, const Domain &domain);
//       → MPI_Sendrecv calls to exchange pressure/velocity boundary layers
//         between neighbouring subdomains.
//
// The #include "Communication.hpp" directives already present in
// Fields.cpp and PressureSolver.cpp serve as call-site placeholders.
#include "Communication.hpp"

#include <vector>
void Communication::init_parallel(int& argn, char **&args){
    // Initialize MPI: assign rank numbers and define global communicator
    MPI_Init(&argn, &args); 
    
    // Get the number of processes in the communicator
    // -> MPI_COMM_WORLD is the global communicator comprising all processes  
    MPI_Comm_size(MPI_COMM_WORLD, &_size); 
    
    // Get the assigned rank number of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
}

void Communication::finalize(){
    // Wait for all processes to exit
    MPI_Finalize();
}

int Communication::get_size(){
    return _size;
}

int Communication::get_rank(){
    return _rank;
}

void Communication::communicate_field(Matrix<double> &field, const Domain &domain){
    // MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
    //     int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source,
    //     int recvtag, MPI_Comm comm, MPI_Status *status)
    int n_col = field.num_cols();
    int n_row = field.num_rows();

    // Send upwards and receive downwards
    if (domain.bottom_physical && !domain.top_physical) {
        MPI_Send(field.data() + (n_row - 2) * n_col + 1, n_col - 2, MPI_DOUBLE, domain.rank_top, 0, MPI_COMM_WORLD);
    }
    if (domain.top_physical && !domain.bottom_physical) {
        MPI_Recv(field.data() + 1, n_col - 2, MPI_DOUBLE, domain.rank_bottom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (!domain.top_physical && !domain.bottom_physical) {
        MPI_Sendrecv(field.data() + (n_row - 2) * n_col + 1, n_col - 2, MPI_DOUBLE, domain.rank_top, 0,
                    field.data() + 1, n_col - 2, MPI_DOUBLE, domain.rank_bottom, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Send downwards and receive upwards
    if (domain.bottom_physical && !domain.top_physical) {
        MPI_Recv(field.data() + (n_row - 1) * n_col + 1, n_col - 2, MPI_DOUBLE, domain.rank_top, 0, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (domain.top_physical && !domain.bottom_physical) {
        MPI_Send(field.data() + n_col + 1, n_col - 2, MPI_DOUBLE, domain.rank_bottom, 0, MPI_COMM_WORLD);
    }
    if (!domain.top_physical && !domain.bottom_physical) {
        MPI_Sendrecv(field.data() + n_col + 1, n_col - 2, MPI_DOUBLE, domain.rank_bottom, 0,
                    field.data() + (n_row - 1) * n_col + 1, n_col - 2, MPI_DOUBLE, domain.rank_top, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    // Send right and receive left
    if (domain.left_physical && !domain.right_physical) {
        std::vector<double> col_vec = field.get_col(n_col - 2);
        MPI_Send(col_vec.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_right, 0, MPI_COMM_WORLD);
    }
    if (domain.right_physical && !domain.left_physical) {
        std::vector<double> col_vec(n_row, 0);
        MPI_Recv(col_vec.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_left, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        field.set_col(col_vec, 0);
    }
    if (!domain.left_physical && !domain.right_physical) {
        std::vector<double> col_vec_send = field.get_col(n_col - 2);
        std::vector<double> col_vec_receive(n_row, 0);

        MPI_Sendrecv(col_vec_send.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_right, 0,
                     col_vec_receive.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_left, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        field.set_col(col_vec_receive, 0);
    }

    // Send right and receive left
    if (domain.left_physical && !domain.right_physical) {
        std::vector<double> col_vec = field.get_col(n_col - 2);
        MPI_Send(col_vec.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_right, 0, MPI_COMM_WORLD);
    }
    if (domain.right_physical && !domain.left_physical) {
        std::vector<double> col_vec(n_row, 0);
        MPI_Recv(col_vec.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_left, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        field.set_col(col_vec, 0);
    }
    if (!domain.left_physical && !domain.right_physical) {
        std::vector<double> col_vec_send = field.get_col(n_col - 2);
        std::vector<double> col_vec_receive(n_row, 0);

        MPI_Sendrecv(col_vec_send.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_right, 0,
                     col_vec_receive.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_left, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        field.set_col(col_vec_receive, 0);
    }

    // Send left and receive right
    if (!domain.left_physical && domain.right_physical) {
        std::vector<double> col_vec = field.get_col(1);
        MPI_Send(col_vec.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_left, 0, MPI_COMM_WORLD);
    }
    if (!domain.right_physical && domain.left_physical) {
        std::vector<double> col_vec(n_row, 0);
        MPI_Recv(col_vec.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_right, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        field.set_col(col_vec, n_col - 1);
    }
    if (!domain.left_physical && !domain.right_physical) {
        std::vector<double> col_vec_send = field.get_col(1);
        std::vector<double> col_vec_receive(n_row, 0);

        MPI_Sendrecv(col_vec_send.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_left, 0,
                     col_vec_receive.data() + 1, n_row - 2, MPI_DOUBLE, domain.rank_right, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        field.set_col(col_vec_receive, n_col - 1);
    }
}
double Communication::reduce_min(double value){
    // We want to calculate the minimum of all values in different processes
    // We use all reduce to "collect" values from all ranks, we want to read a single addres. 
    // The output type will be double and the operation performed is the minimum of all values
    // and returns this minimum to all processes and overwrites in the initial address
    MPI_Allreduce(&value, &value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return value;
}
double Communication::reduce_sum(double value){
    MPI_Allreduce(&value, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return value;
}
