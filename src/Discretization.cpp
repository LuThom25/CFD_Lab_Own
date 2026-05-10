#include "Discretization.hpp"

#include <cassert>
#include <cmath>

double Discretization::_dx    = 0.0;
double Discretization::_dy    = 0.0;
double Discretization::_gamma = 0.0;
// P3: static precomputed reciprocals initialised alongside _dx/_dy.
double Discretization::_dx2_inv = 0.0;
double Discretization::_dy2_inv = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx    = dx;
    _dy    = dy;
    _gamma = gamma;
    // P3: compute once here; every laplacian/sor_helper call reuses these.
    _dx2_inv = 1.0 / (dx * dx);
    _dy2_inv = 1.0 / (dy * dy);
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // R4: debug-mode guard — indices must be interior (not ghost cells).
    assert(i >= 1 && i <= U.num_cols() - 2);
    assert(j >= 1 && j <= U.num_rows() - 2);

    // Implementing Eq. 4: d(u^2)/dx + d(uv)/dy using the mixed central/donor-cell scheme
    // with upwinding coefficient gamma. U(i,j) sits at the right face of cell (i,j).

    // --- d(u^2)/dx ---
    // P2: fast() bypasses bounds check — indices proven valid by the assert above.
    double u_e = (U.fast(i, j) + U.fast(i + 1, j)) / 2.0;
    double u_w = (U.fast(i - 1, j) + U.fast(i, j)) / 2.0;
    double du2_dx = (u_e * u_e - u_w * u_w) / _dx
                  + _gamma / _dx * (
                        std::abs(U.fast(i, j) + U.fast(i + 1, j)) * (U.fast(i, j) - U.fast(i + 1, j)) / 4.0
                      - std::abs(U.fast(i - 1, j) + U.fast(i, j)) * (U.fast(i - 1, j) - U.fast(i, j)) / 4.0);

    // --- d(uv)/dy ---
    double v_n  = (V.fast(i, j) + V.fast(i + 1, j)) / 2.0;
    double v_s  = (V.fast(i, j - 1) + V.fast(i + 1, j - 1)) / 2.0;
    double u_n  = (U.fast(i, j) + U.fast(i, j + 1)) / 2.0;
    double u_s  = (U.fast(i, j - 1) + U.fast(i, j)) / 2.0;
    double duv_dy = (v_n * u_n - v_s * u_s) / _dy
                  + _gamma / _dy * (
                        std::abs(V.fast(i, j) + V.fast(i + 1, j))         * (U.fast(i, j) - U.fast(i, j + 1))     / 4.0
                      - std::abs(V.fast(i, j - 1) + V.fast(i + 1, j - 1)) * (U.fast(i, j - 1) - U.fast(i, j))     / 4.0);

    return du2_dx + duv_dy;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // R4: debug-mode guard — indices must be interior (not ghost cells).
    assert(i >= 1 && i <= V.num_cols() - 2);
    assert(j >= 1 && j <= V.num_rows() - 2);

    // Implementing Eq. 5: d(uv)/dx + d(v^2)/dy using the mixed central/donor-cell scheme
    // with upwinding coefficient gamma. V(i,j) sits at the top face of cell (i,j).

    // --- d(uv)/dx ---
    // P2: fast() bypasses bounds check — indices proven valid by the assert above.
    double u_e  = (U.fast(i, j) + U.fast(i, j + 1)) / 2.0;
    double u_w  = (U.fast(i - 1, j) + U.fast(i - 1, j + 1)) / 2.0;
    double v_e  = (V.fast(i, j) + V.fast(i + 1, j)) / 2.0;
    double v_w  = (V.fast(i - 1, j) + V.fast(i, j)) / 2.0;
    double duv_dx = (u_e * v_e - u_w * v_w) / _dx
                  + _gamma / _dx * (
                        std::abs(U.fast(i, j) + U.fast(i, j + 1))         * (V.fast(i, j) - V.fast(i + 1, j))   / 4.0
                      - std::abs(U.fast(i - 1, j) + U.fast(i - 1, j + 1)) * (V.fast(i - 1, j) - V.fast(i, j))   / 4.0);

    // --- d(v^2)/dy ---
    double v_n  = (V.fast(i, j) + V.fast(i, j + 1)) / 2.0;
    double v_s  = (V.fast(i, j - 1) + V.fast(i, j)) / 2.0;
    double dv2_dy = (v_n * v_n - v_s * v_s) / _dy
                  + _gamma / _dy * (
                        std::abs(V.fast(i, j) + V.fast(i, j + 1))   * (V.fast(i, j) - V.fast(i, j + 1))   / 4.0
                      - std::abs(V.fast(i, j - 1) + V.fast(i, j))   * (V.fast(i, j - 1) - V.fast(i, j))   / 4.0);

    return duv_dx + dv2_dy;
}

double Discretization::convection_T(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i, int j){
    double uT_dx;
    double vT_dy;

    uT_dx = 0.5 * (U(i,j) * (T(i,j) + T(i+1,j)) - U(i-1,j) * (T(i-1,j) + T(i,j))) / _dx 
            + 0.5 * _gamma * (std::abs(U(i,j)) * (T(i,j) - T(i+1,j)) - std::abs(U(i-1,j)) * (T(i-1,j) - T(i,j))) / _dx;
    vT_dy = 0.5 * (V(i,j) * (T(i,j) + T(i,j+1)) - V(i,j-1) * (T(i,j-1) + T(i,j))) / _dy 
            + 0.5 * _gamma * (std::abs(V(i,j)) * (T(i,j) - T(i,j+1)) - std::abs(V(i,j-1)) * (T(i,j-1) - T(i,j))) / _dy;

    return uT_dx + vT_dy;
}

// We can reuse the laplacian discretization inserting T in A
double Discretization::laplacian(const Matrix<double> &A, int i, int j) {
    // Central difference discretization of the Laplace operator: d^2A/dx^2 + d^2A/dy^2
    // P2: fast() for unchecked access; P3: multiply by precomputed reciprocal.
    return (A.fast(i + 1, j) - 2.0 * A.fast(i, j) + A.fast(i - 1, j)) * _dx2_inv
         + (A.fast(i, j + 1) - 2.0 * A.fast(i, j) + A.fast(i, j - 1)) * _dy2_inv;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    // Off-diagonal part of the Laplacian stencil (excludes the -2*P(i,j) centre terms).
    // P2: fast() for unchecked access; P3: multiply by precomputed reciprocal.
    return (P.fast(i + 1, j) + P.fast(i - 1, j)) * _dx2_inv
         + (P.fast(i, j + 1) + P.fast(i, j - 1)) * _dy2_inv;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    // Linear interpolation: arithmetic mean between A(i,j) and A(i+i_offset, j+j_offset)
    return (A(i, j) + A(i + i_offset, j + j_offset)) / 2.0;
}
