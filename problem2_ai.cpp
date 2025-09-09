#include <bits/stdc++.h>
#include <print>
#include "problem1_ai.h"

// Problem 2: Solve for r1, r2, r3 (Earth-Sun collinear Lagrange points) using the
// root-finding routines from Problem 1.
//
// Equations (set F(r) = 0):
//  r1 in (0, R):
//    GM1/(R - r1)^2 - GM2/r1^2 - ( (M1/(M1+M2))*R - r1 ) * G(M1+M2)/R^3 = 0
//
//  r2 > 0:
//    GM1/(R + r2)^2 + GM2/r2^2 - ( (M1/(M1+M2))*R + r2 ) * G(M1+M2)/R^3 = 0
//
//  r3 small compared to R (Sun-side distance is (R - r3)):
//    GM1/(R - r3)^2 + GM2/(2R - r3)^2 - ( (M2/(M1+M2))*R + R - r3 ) * G(M1+M2)/R^3 = 0
//
// We use Newton-Raphson with a numerical derivative for robustness and speed.

int main() {
    using namespace hw2;

    // Physical constants (SI)
    // G: CODATA 2018 value
    const double G   = 6.67430e-11;           // m^3 kg^-1 s^-2
    const double M1  = 1.98847e30;            // kg (mass of the Sun)
    const double M2  = 5.9722e24;             // kg (mass of the Earth)
    const double R   = 1.495978707e11;        // m (mean Earth-Sun distance, 1 AU)

    const double Msum = M1 + M2;
    const double omega2 = G * Msum / (R*R*R); // from Kepler's third law

    // Helper: central-difference numerical derivative
    auto num_deriv = [](auto&& F, double x) {
        // Choose h relative to x to balance truncation/roundoff
        double h = std::max(1e-6, std::abs(x) * 1e-6);
        double f1 = F(x + h);
        double f2 = F(x - h);
        return (f1 - f2) / (2.0 * h);
    };

    // --- r1 equation ---
    auto F1 = [&](double r1) -> double {
        // guard against invalid ranges
        if (r1 <= 0.0 || r1 >= R) return std::numeric_limits<double>::infinity();
        double termG = G * (M1 / ((R - r1)*(R - r1)) - M2 / (r1*r1));
        double r_orbit = (M1 / Msum) * R - r1;
        double termC = r_orbit * (G * Msum) / (R*R*R);
        return termG - termC;
    };
    auto F1p = [&](double r1) -> double { return num_deriv(F1, r1); };

    // Initial guess near the known scale (~1.5e9 m)
    double r1_guess = 1.5e9;

    auto [ok1, r1] = newton_raphson(F1, F1p, r1_guess, 1e-12);
    if (!ok1) {
        // Fallback: try bisection in (0, R). We need a bracket that changes sign.
        // We'll expand from a small positive to something before R.
        double a = 1.0;              // avoid exactly 0
        double b = R - 1.0;          // avoid exactly R
        // Search for a sign change by shrinking/expanding around Earth side.
        // We expect r1 to be much smaller than R, so narrow the bracket.
        a = 1e6; b = 5e9; // 1000 km to 5e9 m
        auto [okb, r1b] = bisection(F1, a, b, 1e-12);
        ok1 = okb; r1 = r1b;
    }

    // --- r2 equation ---
    auto F2 = [&](double r2) -> double {
        if (r2 <= 0.0) return std::numeric_limits<double>::infinity();
        double termG = G * (M1 / ((R + r2)*(R + r2)) + M2 / (r2*r2));
        double r_orbit = (M1 / Msum) * R + r2;
        double termC = r_orbit * (G * Msum) / (R*R*R);
        return termG - termC;
    };
    auto F2p = [&](double r2) -> double { return num_deriv(F2, r2); };

    double r2_guess = 1.5e9;
    auto [ok2, r2] = newton_raphson(F2, F2p, r2_guess, 1e-12);
    if (!ok2) {
        double a = 1.0;      // just beyond zero
        double b = 6e9;      // a broad upper bound
        auto [okb, r2b] = bisection(F2, a, b, 1e-12);
        ok2 = okb; r2 = r2b;
    }

    // --- r3 equation ---
    auto F3 = [&](double r3) -> double {
        // r3 can be small but must be < R to keep denominators positive
        if (r3 <= 0.0 || r3 >= R) return std::numeric_limits<double>::infinity();
        double termG = G * ( M1 / ((R - r3)*(R - r3)) + M2 / ((2.0*R - r3)*(2.0*R - r3)) );
        double r_orbit = (M2 / Msum) * R + R - r3;
        double termC = r_orbit * (G * Msum) / (R*R*R);
        return termG - termC;
    };
    auto F3p = [&](double r3) -> double { return num_deriv(F3, r3); };

    // r3 is also on the order of 1e9 m
    double r3_guess = 1.0e9;
    auto [ok3, r3] = newton_raphson(F3, F3p, r3_guess, 1e-12);
    if (!ok3) {
        double a = 1.0e6; // 1000 km
        double b = 1.0e10; // wide
        auto [okb, r3b] = bisection(F3, a, b, 1e-12);
        ok3 = okb; r3 = r3b;
    }

    // Print results with 6 significant figures.
    // Even if an individual solve failed, still print the last estimate per instructions.
    std::println("r1: {:.6g}m", r1);
    std::println("r2: {:.6g}m", r2);
    std::println("r3: {:.6g}m", r3);

    // Exit code 0 regardless; the tuple's boolean is how we signal success in library functions.
    return 0;
}