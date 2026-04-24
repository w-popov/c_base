#ifndef KRAMER_HPP
#define KRAMER_HPP

#include <array>
#include <iostream>

namespace Kramer 
{

constexpr size_t SIZEM = 3;
using matrix3x3 = std::array<std::array<double, SIZEM>, SIZEM>;
using vector3 = std::array<double, SIZEM>;

class Kramer
{
public:
    // Конструктор по умолчанию удалён, так как класс требует инициализации данными
    Kramer() = delete;

    // Конструктор инициализации
    Kramer(const matrix3x3& coefficients, const vector3& constants) :
        coefficients_(coefficients), constants_(constants) {}

    // Конструктор копирования
    Kramer(const Kramer& other) :
        coefficients_(other.coefficients_), constants_(other.constants_) {}

    // Конструктор перемещения
    Kramer(Kramer&& other) noexcept :
        coefficients_(std::move(other.coefficients_)), constants_(std::move(other.constants_)) {}

    // Оператор присваивания копирования
    Kramer& operator=(const Kramer& other)
    {
        if (this != &other)
        {
            coefficients_ = other.coefficients_;
            constants_ = other.constants_;
        }
        return *this;
    }

    // Оператор присваивания перемещения
    Kramer& operator=(Kramer&& other) noexcept
    {
        if (this != &other)
        {
            coefficients_ = std::move(other.coefficients_);
            constants_ = std::move(other.constants_);
        }
        return *this;
    }

    // Вычисление определителя матрицы методом Гаусса
    double computeDeterminant(matrix3x3 matrix) const;

    // Решение системы линейных уравнений методом Крамера
    vector3 solve() const;

private:
    matrix3x3 coefficients_;
    vector3 constants_;
};

}

#endif // KRAMER_HPP    