#ifndef KRAMER_HPP
#define KRAMER_HPP

#include <vector>
#include <iostream>

template<typename T = double, typename U = double>
class Kramer {
public:
    // Конструктор по умолчанию
    Kramer() = delete;

    // Конструктор инициализации
    Kramer(const std::vector<std::vector<T>>& coefficients, const std::vector<U>& constants) :
        coefficients_(coefficients), constants_(constants) {}

    // Конструктор копирования
    Kramer(const Kramer& other) :
        coefficients_(other.coefficients_), constants_(other.constants_) {}

    // Конструктор перемещения
    Kramer(Kramer&& other) noexcept :
        coefficients_(std::move(other.coefficients_)), constants_(std::move(other.constants_)) {}

    // Оператор присваивания копирования
    Kramer& operator=(const Kramer& other) {
        if (this != &other) {
            coefficients_ = other.coefficients_;
            constants_ = other.constants_;
        }
        return *this;
    }

    // Оператор присваивания перемещения
    Kramer& operator=(Kramer&& other) noexcept {
        if (this != &other) {
            coefficients_ = std::move(other.coefficients_);
            constants_ = std::move(other.constants_);
        }
        return *this;
    }

    // Вычисление определителя матрицы (вспомогательный метод)
    T computeDeterminant(std::vector<std::vector<T>> matrix) const {
        if (matrix.empty()) {
            return 0;
        }

        size_t n = matrix.size();
        T det = 1;
        int swaps = 0;

        // Приведение матрицы к треугольному виду методом Гаусса
        for (size_t i = 0; i < n; ++i) {
            // Поиск ведущего элемента
            size_t maxRow = i;
            for (size_t k = i + 1; k < n; ++k) {
                if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Если ведущий элемент близок к нулю, матрица вырождена
            if (std::abs(matrix[maxRow][i]) < 1e-10) {
                return 0;
            }

            // Если нужно было переставить строки
            if (maxRow != i) {
                std::swap(matrix[i], matrix[maxRow]);
                swaps++;
            }

            // Обновляем определитель
            det *= matrix[i][i];

            // Исключение переменной из остальных строк
            for (size_t k = i + 1; k < n; ++k) {
                T factor = matrix[k][i] / matrix[i][i];
                for (size_t j = i; j < n; ++j) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }

        // Учитываем количество перестановок строк
        if (swaps % 2 == 1) {
            det = -det;
        }

        return det;
    }

    // Решение системы линейных уравнений методом Крамера
    std::vector<T> solve() const {
        size_t n = coefficients_.size();
        std::vector<T> solution(n);

        T det = computeDeterminant(coefficients_);
        if (std::abs(det) < 1e-10) {
            throw std::runtime_error("Система не имеет единственного решения (определитель равен нулю).");
        }

        for (size_t i = 0; i < n; ++i) {
            std::vector<std::vector<T>> modifiedMatrix = coefficients_;
            for (size_t j = 0; j < n; ++j) {
                modifiedMatrix[j][i] = static_cast<T>(constants_[j]);
            }
            solution[i] = computeDeterminant(modifiedMatrix) / det;
        }

        return solution;
    }

private:
    std::vector<std::vector<T>> coefficients_;
    std::vector<U> constants_;

};

#endif // KRAMER_HPP    