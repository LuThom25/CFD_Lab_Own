// Communication.hpp — Placeholder for Worksheet 2 (MPI parallelisation)
//
// This header is intentionally empty in Worksheet 1, which runs on a single
// process only.  In a parallel extension (WS2), this class would declare
// communication routines for domain decomposition:
//
//   class Communication {
//   public:
//       // Exchange ghost-cell (halo) values between neighbouring MPI ranks.
//       // Called after each SOR sweep and after flux/velocity updates.
//       static void communicate(Fields &field, const Domain &domain);
//   };
//
// Domain.hpp already defines the Domain struct (iminb/imaxb/jminb/jmaxb,
// domain_imax/jmax) and includes <mpi.h> in preparation for this extension.
// Fields.cpp and PressureSolver.cpp already include this header as
// call-site placeholders.
#pragma once
