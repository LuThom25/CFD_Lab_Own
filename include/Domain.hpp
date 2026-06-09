#pragma once
#include "Enums.hpp"
#include <mpi.h>

/**
 * @brief Data structure that holds geometrical information
 * necessary for decomposition.
 *
 */
struct Domain {

    // These indexes would only used for cell type assingment
    /// Minimum x index including ghost cells
    int iminb{-1};
    /// Maximum x index including ghost cells
    int imaxb{-1};

    /// Minimum y index including ghost cells
    int jminb{-1};
    /// Maximum y index including ghost cells
    int jmaxb{-1};

    /// Cell length
    double dx{-1.0};
    /// Cell height
    double dy{-1.0};

    /// Number of cells in x direction (without ghost), used in grid with +2 to create matrix of cells
    int size_x{-1};
    /// Number of cells in y direction
    int size_y{-1};

    /// Number of cells in x direction, not-decomposed (without ghost)
    int domain_imax{-1};
    /// Number of cells in y direction, not-decomposed
    int domain_jmax{-1};

    // WS3: neighbour MPI ranks (MPI_PROC_NULL if no neighbour on that side)
    int rank_left{MPI_PROC_NULL};
    int rank_right{MPI_PROC_NULL};
    int rank_top{MPI_PROC_NULL};
    int rank_bottom{MPI_PROC_NULL};

    // WS3: true when this subdomain borders a physical domain boundary
    bool left_physical{true};
    bool right_physical{true};
    bool top_physical{true};
    bool bottom_physical{true};
};
