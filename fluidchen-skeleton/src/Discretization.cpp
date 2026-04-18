#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // Implementing Eq. 4: d(u^2)/dx + d(uv)/dy using the mixed central/donor-cell scheme
    // with upwinding coefficient gamma. U(i,j) sits at the right face of cell (i,j).

    // --- d(u^2)/dx ---
    // Interpolated u at east and west faces of the u-control volume
    double u_e = (U(i, j) + U(i + 1, j)) / 2.0;
    double u_w = (U(i - 1, j) + U(i, j)) / 2.0;
    double du2_dx = (u_e * u_e - u_w * u_w) / _dx
                  + _gamma / _dx * (
                        std::abs(U(i, j) + U(i + 1, j)) * (U(i, j) - U(i + 1, j)) / 4.0
                      - std::abs(U(i - 1, j) + U(i, j)) * (U(i - 1, j) - U(i, j)) / 4.0);

    // --- d(uv)/dy ---
    // Interpolated v and u at north and south faces of the u-control volume
    double v_n  = (V(i, j) + V(i + 1, j)) / 2.0;
    double v_s  = (V(i, j - 1) + V(i + 1, j - 1)) / 2.0;
    double u_n  = (U(i, j) + U(i, j + 1)) / 2.0;
    double u_s  = (U(i, j - 1) + U(i, j)) / 2.0;
    double duv_dy = (v_n * u_n - v_s * u_s) / _dy
                  + _gamma / _dy * (
                        std::abs(V(i, j) + V(i + 1, j))     * (U(i, j) - U(i, j + 1))     / 4.0
                      - std::abs(V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) - U(i, j)) / 4.0);

    return du2_dx + duv_dy;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // Implementing Eq. 5: d(uv)/dx + d(v^2)/dy using the mixed central/donor-cell scheme
    // with upwinding coefficient gamma. V(i,j) sits at the top face of cell (i,j).

    // --- d(uv)/dx ---
    // Interpolated u and v at east and west faces of the v-control volume
    double u_e  = (U(i, j) + U(i, j + 1)) / 2.0;
    double u_w  = (U(i - 1, j) + U(i - 1, j + 1)) / 2.0;
    double v_e  = (V(i, j) + V(i + 1, j)) / 2.0;
    double v_w  = (V(i - 1, j) + V(i, j)) / 2.0;
    double duv_dx = (u_e * v_e - u_w * v_w) / _dx
                  + _gamma / _dx * (
                        std::abs(U(i, j) + U(i, j + 1))     * (V(i, j) - V(i + 1, j))   / 4.0
                      - std::abs(U(i - 1, j) + U(i - 1, j + 1)) * (V(i - 1, j) - V(i, j)) / 4.0);

    // --- d(v^2)/dy ---
    // Interpolated v at north and south faces of the v-control volume
    double v_n  = (V(i, j) + V(i, j + 1)) / 2.0;
    double v_s  = (V(i, j - 1) + V(i, j)) / 2.0;
    double dv2_dy = (v_n * v_n - v_s * v_s) / _dy
                  + _gamma / _dy * (
                        std::abs(V(i, j) + V(i, j + 1))   * (V(i, j) - V(i, j + 1))   / 4.0
                      - std::abs(V(i, j - 1) + V(i, j))   * (V(i, j - 1) - V(i, j))   / 4.0);

    return duv_dx + dv2_dy;
}

double Discretization::laplacian(const Matrix<double> &A, int i, int j) {
    // Central difference discretization of the Laplace operator: d^2A/dx^2 + d^2A/dy^2
    return (A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j)) / (_dx * _dx)
         + (A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1)) / (_dy * _dy);
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    // Linear interpolation: arithmetic mean between A(i,j) and A(i+i_offset, j+j_offset)
    return (A(i, j) + A(i + i_offset, j + j_offset)) / 2.0;
}