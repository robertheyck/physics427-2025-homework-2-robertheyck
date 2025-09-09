#include <print>
#include <tuple>
#include <cmath>
#include <limits>

template <typename F>
// Looking for f(x) = g(x) + x (want zero of f), so x = -g(x) --> the update rule is x = x - f(x)
std::tuple<bool, double> fixed_point_iteration(const F& f, double x0, double tolerance, std::size_t max_iters = 10000) {
    double x = x0;
    bool success = false;
    for(int i = 0; i < max_iters; i++) {
        // Check if x is NaN or infinite
        if(!std::isfinite(x)) {
            std::println("Failed because iteration led to x = {}", x);   // NaN or inf
            std::tuple<bool, double> sol_bad = std::make_tuple(success, x);
            return sol_bad;
        }

        if (std::abs(f(x)) < tolerance) {
            success = true;
            std::tuple<bool, double> sol_good = std::make_tuple(success, x);
            return sol_good;
        }
        x = x - f(x);
    }  
    std::println("Method did not converge.");
    std::tuple<bool, double> sol_bad = std::make_tuple(success, x);
    return sol_bad;
}

template <typename F>
std::tuple<bool, double> bisection(const F& f, double a, double b, double tolerance, std::size_t max_iters = 10000) {
    // First check whether f(a) and f(b) have the same sign
    double x1 = a;
    double x2 = b;
    bool success = false;

    // First, check whether the endpoints of the interval are valid (f must change sign)
    if(f(x1) * f(x2) < 0.0) {          
        // Regular method 
        for(int i = 0; i < max_iters; i++) {
            // Check if x1 and/or x2 is NaN or infinite
            if(!std::isfinite(x1)) {
                std::println("Failed because iteration led to x1 = {}", x1);   // NaN or inf
                std::tuple<bool, double> sol_bad = std::make_tuple(success, x1);
                return sol_bad;
            } else if (!std::isfinite(x2)) {
                std::println("Failed because iteration led to x2 = {}", x2);   // NaN or inf
                std::tuple<bool, double> sol_bad = std::make_tuple(success, x2);
                return sol_bad;
            }
            // Evaluate beforehand to speed up
            double f1 = f(x1);
            double f2 = f(x2);

            // Check if either x1 or x2 is a zero
            if(std::abs(f1) < tolerance) {
                success = true;
                std::tuple<bool, double> sol_good = std::make_tuple(success, x1);
                return sol_good;
            } else if (std::abs(f2) < tolerance) {
                success = true;
                std::tuple<bool, double> sol_good = std::make_tuple(success, x2);
                return sol_good;
            }

            // Update rule
            double x_mid = 0.5 * (x1 + x2);
            double f_mid = f(x_mid);
            // Check if x_mid works before updating
            if(std::abs(f_mid) < tolerance) {
                success = true;
                std::tuple<bool, double> sol_good = std::make_tuple(success, x_mid);
                return sol_good;
            }
            // Replace one of the endpoints by x_mid
            if(f1 * f_mid < 0.0) {       
                x2 = x_mid;
            } else if(f_mid * f2 < 0.0) {     
                x1 = x_mid;
            }
        }
        // If no solution is found after maximum number of iterations
        std::println("Method did not converge.");
        std::tuple<bool, double> sol_bad = std::make_tuple(success, 0.5 * (x1 + x2));   // Could also use x1 or x2 instead, arbitrary
        return sol_bad;
    }
    // This code only executes if endpoints are invalid
    std::println("Initial points are invalid because f does not change sign.");
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::tuple<bool, double> sol_bad = std::make_tuple(success, nan);
    return sol_bad;
}

template <typename F, typename Fprime>
std::tuple<bool, double> newton_raphson(const F& f, const Fprime& fprime, double x0, double tolerance, std::size_t max_iters = 10000) {
    double x = x0;
    bool success = false;
    for(int i = 0; i < max_iters; i++) {
        // Check if x is NaN or infinite
        if(!std::isfinite(x)) {
            std::println("Failed because iteration led to x = {}", x);   // NaN or inf
            std::tuple<bool, double> sol_bad = std::make_tuple(success, x);
            return sol_bad;
        }
        // Check derivative not zero 
        float eps_deriv = std::numeric_limits<float>::epsilon();
        if(std::abs(fprime(x)) < eps_deriv) {
            std::println("Failed because the derivative of f is zero at x = {}.", x);
            double nan = std::numeric_limits<double>::quiet_NaN();
            std::tuple<bool, double> sol_bad = std::make_tuple(success, nan);
            return sol_bad;
        }
        // Check if satisfied
        if(std::abs(f(x)) < tolerance) {
            success = true;
            std::tuple<bool, double> sol_good = std::make_tuple(success, x);
            return sol_good;
        }
        // Update rule
        x -= f(x) / fprime(x);
    }
    // If fails to converge
    std::println("Failed to converge.");
    std::tuple<bool, double> sol_bad = std::make_tuple(success, x);
    return sol_bad;
}
