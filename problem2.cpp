#include "problem1.h"
#include <iomanip>
#include <print>
#include <cmath>

// Using bisection method for all three Lagrange points because I don't want to take derivatives

// Function for finding an interval where the sign changes (instead of just guessing)
template <typename F>
std::tuple<bool, double, double> find_interval(const F& f, double initial_guess, double increment, std::size_t max_iters = 20) {
    double f0 = f(initial_guess);
    double a = initial_guess;
    double b = initial_guess;
    bool success = false;

    double var_incr1 = increment;
    double var_incr2 = increment;
    double scale = 10.0;
    for (int i = 1; i < max_iters; i++) {
        
        // To account for huge order of magnitude differences, the increment increases exponentially,
        // then stabilizes or decreases as needed to prevent the lower bound from going to zero.
        if (a - var_incr1 > 0.0 && a - var_incr1 * scale > 0.0) {
            a -= var_incr1;
            var_incr1 *= scale;
        } else if (a - var_incr1 > 0.0 && a - var_incr1 * scale <= 0.0) {
            a -= var_incr1;
        } else {
            var_incr1 *= 1.0 / scale;
        }
        
        // The upper bound can grow without issue, but need a reasonable initial_guess and should keep max_iters small
        // to prevent finding extraneous solutions.
        b += var_incr2;
        var_incr2 *= scale;

        double fa = f(a);
        double fb = f(b);

        if (fa * f0 < 0.0) {
            success = true;
            std::tuple<bool, double, double> interval = std::make_tuple(success, a, initial_guess);
            return interval;
        } else if (f0 * fb < 0.0) {
            success = true;
            std::tuple<bool, double, double> interval = std::make_tuple(success, initial_guess, b);
            return interval;
        }
    }
    std::println("No valid interval found.");
    double nan1 = std::numeric_limits<double>::quiet_NaN();
    double nan2 = std::numeric_limits<double>::quiet_NaN();
    std::tuple<bool, double, double> bad_interval = std::make_tuple(success, nan1, nan2);
    return bad_interval;
}

int main() {
    // Define constants in SI units (M1 is the sun, M2 is earth)
    const double G = 6.6743e-11;
    const double M1 = 1.988416e30;
    const double M2 = 5.9722e24;
    const double R = 1.495979e11; 

    const double mu1 = G * M1;
    const double mu2 = G * M2;
    const double omega_sq = G * (M1 + M2) / (R * R * R);
    const double CoM_to_earth = R * M1 / (M1 + M2);
    const double CoM_to_sun = R * M2 / (M1 + M2);

    auto f_L1 = [=](double r1) { 
        double d = R - r1;
        double f = omega_sq * (CoM_to_earth - r1) + mu2 / (r1 * r1) - mu1 / (d * d);
        return f; 
    };
    
    auto f_L2 = [=](double r2) { 
        double d = R + r2;
        double f = omega_sq * (CoM_to_earth + r2) - mu2 / (r2 * r2) - mu1 / (d * d);
        return f; 
    };

    auto f_L3 = [=](double r3) { 
        double d1 = R - r3;
        double d2 = d1 + R;
        double f = omega_sq * (CoM_to_sun + d1) - mu2 / (d2 * d2) - mu1 / (d1 * d1);
        return f; 
    };

    // For L1 and L2, using same initial guess and increment size 
    double initial_guess1 = 0.5 * R;
    double incr1 = CoM_to_sun;

    auto interval1 = find_interval(f_L1, initial_guess1, incr1);
    double guess1_r1 = std::get<1>(interval1);
    double guess2_r1 = std::get<2>(interval1);
    auto L1 = bisection(f_L1, guess1_r1, guess2_r1, 1e-10);
    bool success1 = std::get<0>(L1);
    double r1 = std::get<1>(L1);

    auto interval2 = find_interval(f_L2, initial_guess1, incr1);
    double guess1_r2 = std::get<1>(interval2);
    double guess2_r2 = std::get<2>(interval2);
    auto L2 = bisection(f_L2, guess1_r2, guess2_r2, 1e-10);
    bool success2 = std::get<0>(L2);
    double r2 = std::get<1>(L2);

    // Turns out L3 is many magnitudes smaller than L1 and L2, so better to use different guess and increment size
    double initial_guess2 = CoM_to_sun;
    double incr2 = 0.1 * CoM_to_sun;

    auto interval3 = find_interval(f_L3, initial_guess2, incr2);
    double guess1_r3 = std::get<1>(interval3);
    double guess2_r3 = std::get<2>(interval3);

    auto L3 = bisection(f_L3, guess1_r3, guess2_r3, 1e-10);
    bool success3 = std::get<0>(L3);
    double r3 = std::get<1>(L3);

    std::println("r1: {:.6g}m", r1);
    std::println("r2: {:.6g}m", r2);
    std::println("r3: {:.6g}m", r3);

    return 0;
}
