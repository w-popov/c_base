#include "kramer.hpp"

namespace Kramer
{

// Вычисление определителя матрицы (вспомогательный метод)
double Kramer::computeDeterminant(matrix3x3 matrix) const
{
    if (matrix.empty())
    {
        return 0;
    }

    size_t n = matrix.size();
    double det = 1;
    int swaps = 0;

    // Приведение матрицы к треугольному виду методом Гаусса
    for (size_t i = 0; i < n; ++i)
    {
        // Поиск ведущего элемента
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k)
        {
            if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i]))
            {
                maxRow = k;
            }
        }

        // Если ведущий элемент близок к нулю, матрица вырождена
        if (std::abs(matrix[maxRow][i]) < 1e-10)
        {
            return 0;
        }

        // Если нужно было переставить строки
        if (maxRow != i)
        {
            std::swap(matrix[i], matrix[maxRow]);
            swaps++;
        }

        // Обновить определитель
        det *= matrix[i][i];

        // Исключение переменной из остальных строк
        for (size_t k = i + 1; k < n; ++k)
        {
            double factor = matrix[k][i] / matrix[i][i];
            for (size_t j = i; j < n; ++j)
            {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }

    // Учитывать количество перестановок строк
    if (swaps % 2 == 1)
    {
        det = -det;
    }

    return det;
}

// Решение системы линейных уравнений методом Крамера
vector3 Kramer::solve() const
{
    size_t n = coefficients_.size();
    vector3 solution = {0};

    double det = computeDeterminant(coefficients_);
    if (std::abs(det) < 1e-10)
    {
        throw std::runtime_error("Система не имеет единственного решения (определитель равен нулю).");
    }

    for (size_t i = 0; i < n; ++i)
    {
        matrix3x3 modifiedMatrix = coefficients_;
        for (size_t j = 0; j < n; ++j)
        {
            modifiedMatrix[j][i] = static_cast<double>(constants_[j]);
        }
        solution[i] = computeDeterminant(modifiedMatrix) / det;
    }

    return solution;
}

}

