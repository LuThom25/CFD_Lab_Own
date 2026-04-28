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
