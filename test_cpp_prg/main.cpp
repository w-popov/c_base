/**
 * Решение СЛАУ 3х3 по формулам Крамера.
 */

#include "kramer.hpp"
#include <iostream>

int main() {
    std::vector<std::vector<double>> coefficients = {
        {1., 2., 3},
        {1., -3., 2.},
        {1., 1., 1.}
    };
    std::vector<double> constants = {7., 5., 3.};
    
    Kramer linear{coefficients, constants};
    auto solution = linear.solve();

    std::cout << "\nSolution: ";
    for (const auto& value : solution) {
        std::cout << "\033[32m" << value << " " << "\033[0m";
    }

    std::cout << std::endl << "Determinant: " << "\033[34m" << 
    linear.computeDeterminant(coefficients) << "\033[0m" << "\n\n";

    return 0;
}