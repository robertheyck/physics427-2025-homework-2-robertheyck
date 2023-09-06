#include "problem1.h"
#include <iomanip>
#include <iostream>

int main() {
  // Test function for fixed point iteration
  auto f1 = [](double x) { return -2.0 + x + std::exp(-x); };

  // Run fixed point iteration
  auto result1 = fixed_point_iteration(f1, 1.0, 1e-6);

  // Check if the result is correct
  double root1 = 1.841406;
  if (std::get<0>(result1) == true &&
      std::abs(std::get<1>(result1) - root1) < 1e-6) {
    std::cout << "Fixed point: PASS\n";
  } else {
    std::cout << "Fixed point: FAIL, root is " << std::setprecision(10)
              << std::get<1>(result1) << ", expected to be around " << root1
              << std::endl;
  }

  // Test function for bisection
  auto f2 = [](double x) { return std::sin(x) - 0.5 * x - 0.1; };

  // Run bisection
  auto result2 = bisection(f2, 0.0, 1.0, 1e-6);

  // Check if the result is correct
  double root2 = 0.202774;
  if (std::get<0>(result2) == true &&
      std::abs(std::get<1>(result2) - root2) < 1e-6) {
    std::cout << "Bisection: PASS\n";
  } else {
    std::cout << "Bisection: FAIL, root is " << std::setprecision(10)
              << std::get<1>(result2) << ", expected to be around " << root2
              << std::endl;
  }

  // Test function and derivative for Newton-Raphson
  auto f3 = [](double x) { return 4.0 + 8.0 * x * x - x * x * x * x; };
  auto f3prime = [](double x) { return 16.0 * x - 4.0 * x * x * x; };

  // Run Newton-Raphson
  auto result3 = newton_raphson(f3, f3prime, 3.0, 1e-6);

  // Check if the result is correct
  double root3 = 2.910693;
  if (std::get<0>(result3) == true &&
      std::abs(std::get<1>(result3) - root3) < 1e-6) {
    std::cout << "Newton: PASS\n";
  } else {
    std::cout << "Newton: FAIL, root is " << std::setprecision(10)
              << std::get<1>(result3) << ", expected to be around " << root3
              << std::endl;
  }

  return 0;
}
