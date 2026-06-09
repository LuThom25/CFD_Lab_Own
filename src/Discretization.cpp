#include "Discretization.hpp"
#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;
// Additional reciprocals initialization
double Discretization::_dx_inv = 0.0;
double Discretization::_dy_inv = 0.0;
double Discretization::_dx2_inv = 0.0;
double Discretization::_dy2_inv = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
    // Pre-computed values
    _dx_inv = 1.0 / dx;
    _dy_inv = 1.0 / dy;
    _dx2_inv = 1.0 / (dx * dx);
    _dy2_inv = 1.0 / (dy * dy);
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // d(u^2)/dx terms
    double u_e = interpolate(U, i, j, 1, 0);     // u at east face
    double u_w = interpolate(U, i - 1, j, 1, 0); // u at west face
    double du2_dx = (u_e * u_e - u_w * u_w) * _dx_inv +
                    _gamma * _dx_inv *
                        (std::abs(u_e) * (U(i, j) - U(i + 1, j)) / 2.0 - std::abs(u_w) * (U(i - 1, j) - U(i, j)) / 2.0);

    // d(uv)/dy terms
    double v_n = interpolate(V, i, j, 1, 0);     // v at north face
    double v_s = interpolate(V, i, j - 1, 1, 0); // v at south face
    double u_n = interpolate(U, i, j, 0, 1);     // u at north face
    double u_s = interpolate(U, i, j - 1, 0, 1); // u at south face
    double duv_dy = (v_n * u_n - v_s * u_s) * _dy_inv +
                    _gamma * _dy_inv *
                        (std::abs(v_n) * (U(i, j) - U(i, j + 1)) / 2.0 - std::abs(v_s) * (U(i, j - 1) - U(i, j)) / 2.0);

    return du2_dx + duv_dy;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // d(uv)/dx terms
    double u_e = interpolate(U, i, j, 0, 1);     // u at east face
    double u_w = interpolate(U, i - 1, j, 0, 1); // u at west face
    double v_e = interpolate(V, i, j, 1, 0);     // v at east face
    double v_w = interpolate(V, i - 1, j, 1, 0); // v at west face
    double duv_dx = (u_e * v_e - u_w * v_w) * _dx_inv +
                    _gamma * _dx_inv *
                        (std::abs(u_e) * (V(i, j) - V(i + 1, j)) / 2.0 - std::abs(u_w) * (V(i - 1, j) - V(i, j)) / 2.0);

    // d(v^2)/dy terms
    double v_n = interpolate(V, i, j, 0, 1);     // v at north face
    double v_s = interpolate(V, i, j - 1, 0, 1); // v at south face
    double dv2_dy = (v_n * v_n - v_s * v_s) * _dy_inv +
                    _gamma * _dy_inv *
                        (std::abs(v_n) * (V(i, j) - V(i, j + 1)) / 2.0 - std::abs(v_s) * (V(i, j - 1) - V(i, j)) / 2.0);

    return duv_dx + dv2_dy;
}

double Discretization::convection_T(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i,
                                    int j) {
    // d(uT)/dx terms
    double uT_dx =
        (U(i, j) * interpolate(T, i, j, 1, 0) - U(i - 1, j) * interpolate(T, i - 1, j, 1, 0)) * _dx_inv +
        _gamma * _dx_inv *
            (std::abs(U(i, j)) * (T(i, j) - T(i + 1, j)) / 2.0 - std::abs(U(i - 1, j)) * (T(i - 1, j) - T(i, j)) / 2.0);

    // d(vT)/dy terms
    double vT_dy =
        (V(i, j) * interpolate(T, i, j, 0, 1) - V(i, j - 1) * interpolate(T, i, j - 1, 0, 1)) * _dy_inv +
        _gamma * _dy_inv *
            (std::abs(V(i, j)) * (T(i, j) - T(i, j + 1)) / 2.0 - std::abs(V(i, j - 1)) * (T(i, j - 1) - T(i, j)) / 2.0);

    return uT_dx + vT_dy;
}

// We can reuse the laplacian discretization inserting T in A
double Discretization::laplacian(const Matrix<double> &A, int i, int j) {
    return (A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j)) * _dx2_inv +
           (A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1)) * _dy2_inv;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    // Off-diagonal part of the Laplacian stencil (excludes the -2*P(i,j) centre terms).
    return (P(i + 1, j) + P(i - 1, j)) * _dx2_inv + (P(i, j + 1) + P(i, j - 1)) * _dy2_inv;
}
double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    return (A(i, j) + A(i + i_offset, j + j_offset)) / 2.0;
}
