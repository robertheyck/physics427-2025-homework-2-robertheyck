#pragma once
#include <tuple>
#include <cmath>
#include <limits>
#include <type_traits>

// Root-finding utilities for Physics 427 HW2
// All functions solve f(x) = 0 to a tolerance meaning |f(x_n)| < epsilon.
// They return {converged, approximate_root}. If inputs are invalid or
// a numerical issue occurs (NaN/Inf, derivative ~ 0, etc.), they return {false, last_x}.

namespace hw2 {

// A small helper to check finiteness of a double-like value.
inline bool is_finite(double x) {
    return std::isfinite(x);
}

// Fixed-point iteration for solving f(x) = 0 using a damped update:
//   x_{n+1} = x_n - lambda * f(x_n)
// with simple adaptive damping to encourage decrease in |f|.
// Note: This is a "fixed-point" style iteration for the root equation itself.
// It is intentionally simple and robust for homework purposes.
template <typename F>
std::tuple<bool, double> fixed_point_iteration(const F& f, double x0, double tolerance) {
    if (!(tolerance > 0) || !is_finite(x0)) {
        return {false, x0};
    }

    constexpr int kMaxIters = 20000;
    constexpr double kMinLambda = 1e-12;
    double x = x0;
    double fx = f(x);
    if (!is_finite(fx)) return {false, x};

    // Start with a moderate damping; adapt it if the residual grows.
    double lambda = 1.0;
    double prev_abs_fx = std::abs(fx);

    for (int iter = 0; iter < kMaxIters; ++iter) {
        if (std::abs(fx) < tolerance) {
            return {true, x};
        }
        double step = -lambda * fx;
        double x_new = x + step;
        if (!is_finite(x_new)) return {false, x};

        double fx_new = f(x_new);
        if (!is_finite(fx_new)) return {false, x_new};

        // If residual decreased, accept and slightly increase lambda (up to 1).
        if (std::abs(fx_new) <= prev_abs_fx) {
            x = x_new;
            fx = fx_new;
            prev_abs_fx = std::abs(fx);
            lambda = std::min(1.0, lambda * 1.5);
        } else {
            // Backtrack / damp more.
            lambda *= 0.5;
            if (lambda < kMinLambda) {
                return {false, x};
            }
            // Do not update x; try again with smaller lambda on next loop.
        }
    }
    // Did not meet tolerance within max iters.
    return {false, x};
}

// Bisection method on [a,b] with f(a) * f(b) < 0 required.
// Terminates when |f(m)| < tolerance (m is midpoint).
template <typename F>
std::tuple<bool, double> bisection(const F& f, double a, double b, double tolerance) {
    if (!(tolerance > 0) || !is_finite(a) || !is_finite(b) || !(b > a)) {
        return {false, std::numeric_limits<double>::quiet_NaN()};
    }
    constexpr int kMaxIters = 2000;

    double fa = f(a);
    double fb = f(b);
    if (!is_finite(fa) || !is_finite(fb)) return {false, std::numeric_limits<double>::quiet_NaN()};
    if (fa == 0.0) return {true, a};
    if (fb == 0.0) return {true, b};
    if (fa * fb >= 0.0) {
        // Invalid bracket
        return {false, (a + b) * 0.5};
    }

    double left = a, right = b;
    double mid = 0.5 * (left + right);
    double fmid = f(mid);
    if (!is_finite(fmid)) return {false, mid};

    for (int iter = 0; iter < kMaxIters; ++iter) {
        mid = 0.5 * (left + right);
        fmid = f(mid);
        if (!is_finite(fmid)) return {false, mid};

        if (std::abs(fmid) < tolerance) {
            return {true, mid};
        }

        // Decide which subinterval contains a sign change
        if (fa * fmid < 0.0) {
            right = mid;
            fb = fmid;
        } else {
            left = mid;
            fa = fmid;
        }
    }
    return {false, 0.5 * (left + right)};
}

// Newton-Raphson: x_{n+1} = x_n - f(x_n)/f'(x_n)
// Terminates when |f(x_n)| < tolerance. Returns false on derivative ~ 0, NaN/Inf, etc.
template <typename F, typename Fprime>
std::tuple<bool, double> newton_raphson(const F& f, const Fprime& fprime, double x0, double tolerance) {
    if (!(tolerance > 0) || !is_finite(x0)) {
        return {false, x0};
    }
    constexpr int kMaxIters = 200;
    constexpr double kTiny = 1e-14;

    double x = x0;
    for (int iter = 0; iter < kMaxIters; ++iter) {
        double fx = f(x);
        if (!is_finite(fx)) return {false, x};
        if (std::abs(fx) < tolerance) {
            return {true, x};
        }
        double dfx = fprime(x);
        if (!is_finite(dfx) || std::abs(dfx) < kTiny) {
            return {false, x};
        }
        double x_new = x - fx / dfx;
        if (!is_finite(x_new)) return {false, x};
        x = x_new;
    }
    return {false, x};
}

} // namespace hw2
