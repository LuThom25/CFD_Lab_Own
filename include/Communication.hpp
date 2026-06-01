// Communication.hpp — Placeholder for Worksheet 2 (MPI parallelisation)
//
// This header is intentionally empty in Worksheet 1, which runs on a single
// process only.  In a parallel extension (WS2), this class would declare
// communication routines for domain decomposition:

// Domain.hpp already defines the Domain struct (iminb/imaxb/jminb/jmaxb,
// domain_imax/jmax) and includes <mpi.h> in preparation for this extension.
// Fields.cpp and PressureSolver.cpp already include this header as
// call-site placeholders.
#pragma once
#include <mpi.h>
#include "Datastructures.hpp"
#include "Domain.hpp"
class Communication {
    private:
        static inline int _size; // total number of processes, each of which owns one domain
        static inline int _rank;

    public:
        static void init_parallel(int& argn, char **&args);
        static void finalize();
        static int get_size();
        static int get_rank();
        // Exchange ghost-cell (halo) values between neighbouring MPI ranks.
        // Called after each field update (either explicit or SOR iteration)
        static void communicate_field(Matrix<double> &field, const Domain &domain);

        static double reduce_min(double value);
        static double reduce_sum(double value);
};
