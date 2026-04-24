/**
 * Решение СЛАУ 3х3 по формулам Крамера.
 */

#include "kramer.hpp"
#include <locale.h>

int main()
{
    std::setlocale(LC_ALL, "Russian");
    Kramer::matrix3x3 coefficients = {
        3., -2., 4.,
        3., 4., -2.,
        2., -1., -1.
    };
    Kramer::vector3 constants = {21., 9., 10.};
    
    Kramer::Kramer linear{coefficients, constants};
    auto solution = linear.solve();

    std::cout << "\n(Решение) Solution: ";
    for (const auto& value : solution)
    {
        std::cout << "\033[32m" << value << " " << "\033[0m";
    }

    std::cout << std::endl << "(Определитель) Determinant: " << "\033[34m" << 
    linear.computeDeterminant(coefficients) << "\033[0m" << "\n\n";

    return 0;
}