#include "catch.hpp"
#include "iostream"
#include <array>
#include <algorithm>
#include "TestScanf.hpp"
#include "func_tools.hpp"
#include "HW8.h"

TEST_CASE("TEST E1: среднее арифметическое всех элементов массива")
{
    constexpr int size = 5;
    SECTION("Секция #1")
    {
        std::array<float, size> arr{4., 15., 3., 10., 14.};
        REQUIRE_THAT(average_array(arr.data(), size),
                     Catch::Matchers::WithinAbs(9.200, 0.001));
    }
    SECTION("Секция #2")
    {
        std::array<float, size> arr{0.02f, 0.8f, 5.5665f, 12.1f, 1};
        REQUIRE_THAT(average_array(arr.data(), size),
                     Catch::Matchers::WithinAbs(3.897, 0.001));
    }
}

TEST_CASE("TEST E3: максимальный и минимальный элементы массива и их номера")
{
    constexpr int SIZE_ARR = 10;
    std::array<int, SIZE_ARR> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("4 -5 3 10 -4 -6 8 -10 1 0");
        struct MinMax mm = min_max_from_array(
                            fill_array(array), 
                            SIZE_ARR);

        /* показать после падения */
        CAPTURE(mm.max_index, mm.max, mm.min_index, mm.min);
        
        REQUIRE((mm.max_index == 4 && mm.max == 10 && mm.min_index == 8 && mm.min == -10));
    }
    
}

TEST_CASE("TEST E4: Найти два максимальных элемента в массиве и напечатать их сумму")
{
    constexpr int size = 10;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("4 -5 3 10 -4 -6 8 -10 1 0");
        REQUIRE(sum_two_max_from_array(fill_array(array), size) == 18);
        array.fill(0);
    }
}

TEST_CASE("TEST E5: посчитать сумму положительных элементов массива")
{
    constexpr int size = 10;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("4 -5 3 10 -4 -6 8 -10 1 0");
        REQUIRE(sum_positive_array(fill_array(array), size) == 26);
        array.fill(0);
    }
}

TEST_CASE("TEST E7: выполнить инверсию отдельно для 1-ой и 2-ой половин массива")
{
    constexpr int size = 10;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("4 -5 3 10 -4 -6 8 -10 1 0");
        half_reverse_array(fill_array(array), size);
        std::array<int, size> res{-4, 10, 3, -5, 4, 0, 1, -10, 8, -6};
        REQUIRE( array == res );
        array.fill(0);
    }
}

TEST_CASE("TEST E8: выполнить инверсию отдельно для каждой трети массива")
{
    constexpr int size = 12;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("1 2 3 4 5 6 7 8 9 10 11 12");
        third_reverse_array(fill_array(array), size);
        std::array<int, size> res{4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9};
        REQUIRE( array == res );
        array.fill(0);
    }
}

TEST_CASE("TEST E9: Считать массив из 10 элементов и выполнить циклический сдвиг ВПРАВО на 1")
{
    constexpr int size = 10;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("1 2 3 4 5 6 7 8 9 10");
        rshift_array(fill_array(array), 1, size);
        std::array<int, size> res{10, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        REQUIRE( array == res );
        array.fill(0);
    }
}

TEST_CASE("TEST E11: Считать массив из 10 элементов и отсортировать его по последней цифре (по возрастанию)")
{
    constexpr int size = 10;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("14  25  13  30  76  58  32  11  41  97");
        sort_X_array(fill_array(array), size, compate_last_digit);
        std::array<int, size> res{30, 11, 41, 32, 13, 14, 25, 76, 97, 58};
        REQUIRE( array == res );
        array.fill(0);
    }
}

TEST_CASE("TEST E13: массив из 10 элементов, вторая с конца цифра(число десятков) – ноль")
{
    constexpr int size = 10;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("40 105 203 1 14 1000 22 33 44 55");
        std::array<int, size> second_array{0};
        std::vector<int> result_array;
        std::vector<int> res{105, 203, 1, 1000};
        round_number_array(fill_array(array), second_array.data(), size);
        std::for_each(second_array.begin(), second_array.end(), 
                    [&result_array](int value){
                        if (value != INT32_MIN)
                            result_array.push_back(value);
                    }
        );
        REQUIRE( result_array ==  res);
        array.fill(0);
    }
}

TEST_CASE("TEST E16: Определить, какое число в массиве встречается чаще всего")
{
    constexpr int size = 10;
    std::array<int, size> array{0};
    SECTION("Секция #1")
    {
        TestScanf ts("4  1  2  1  11  2  34  11  0  11");
        int often = more_often (fill_array(array), size);
        REQUIRE(often == 11);
        array.fill(0);
    }
    SECTION("Секция #1")
    {
        TestScanf ts("1  1  1  1  1  1  1  1  1  1");
        int often = more_often (fill_array(array), size);
        REQUIRE(often == 1);
        array.fill(0);
    }
}

TEST_CASE("TEST E20: Переставить цифры в числе так, что бы получилось максимальное число")
{
    SECTION("Секция #1")
    {
        REQUIRE(the_greatest_number(1229) == 9221);
    }
    SECTION("Секция #2")
    {
        REQUIRE(the_greatest_number(1) == 1);
    }
    SECTION("Секция #3")
    {
        REQUIRE(the_greatest_number(0) == 0);
    }
    SECTION("Секция #4")
    {
        REQUIRE(the_greatest_number(123456) == 654321);
    }
}
