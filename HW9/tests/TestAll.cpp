#include "catch.hpp"
#include <array>
#include <algorithm>
#include "TestScanf.hpp"
#include "func_tools.hpp"
#include "HW9.h"

TEST_CASE("TEST F1: сортирует массив по возрастанию")
{
    SECTION("Секция #1")
    {
        std::vector<int> base{20, 19, 4,  3,  2,  1, 18, 17, 13, 12,
                              11, 16, 15, 14, 10, 9, 8,  7,  6,  5};
        std::vector<int> expected{1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
        sort_array(static_cast<int>(base.size()), base.data());
        REQUIRE(base == expected);
    }
}

TEST_CASE(
    "TEST F2: начало массива все четные элементы, а в конец – все нечетные")
{
    SECTION("Секция #1")
    {
        std::vector<int> base{20, 19, 18, 17, 16, 15, 14, 13, 12, 11,
                              10, 9,  8,  7,  6,  5,  4,  3,  2,  1};
        std::vector<int> expected{20, 18, 16, 14, 12, 10, 8, 6, 4, 2,
                                  19, 17, 15, 13, 11, 9,  7, 5, 3, 1};
        sort_even_odd(static_cast<int>(base.size()), base.data());
        REQUIRE(base == expected);
    }
}

TEST_CASE("TEST F5: максимальный элемент в массиве")
{
    SECTION("Секция #1")
    {
        std::vector<int> base{
            773, 307, 371, 548, 531, 765, 402, 27,  573, 591, 217, 859, 662,
            493, 173, 174, 125, 591, 324, 231, 130, 394, 573, 65,  570, 258,
            343, 3,   586, 14,  785, 296, 140, 726, 598, 262, 807, 794, 510,
            465, 66,  895, 182, 218, 302, 34,  205, 252, 687, 660, 952, 737,
            2,   32,  310, 680, 36,  139, 346, 139, 489, 217, 767, 544, 158,
            774, 883, 154, 802, 136, 569, 377, 867, 423, 224, 176, 118, 660,
            513, 734, 45,  978, 983, 749, 909, 601, 270, 147, 433, 737, 789,
            304, 842, 769, 815, 503, 190, 399, 3};
        int expected = 983;
        int result = find_max_array(static_cast<int>(base.size()), base.data());
        REQUIRE(result == expected);
    }
    SECTION("Секция #2")
    {
        std::vector<int> base{
            -227, -693, -629, -452, -469, -235, -598, -973, -427, -409,
            -783, -141, -338, -507, -827, -826, -875, -409, -676, -769,
            -870, -606, -427, -935, -430, -742, -657, -997, -414, -986,
            -215, -704, -860, -274, -402, -738, -193, -206, -490, -535,
            -934, -105, -818, -782, -698, -966, -795, -748, -313, -340,
            -48,  -263, -998, -968, -690, -320, -964, -861, -654, -861,
            -511, -783, -233, -456, -842, -226, -117, -846, -198, -864,
            -431, -623, -133, -577, -776, -824, -882, -340, -487, -266,
            -955, -22,  -17,  -251, -91,  -399, -730, -853, -567, -263,
            -211, -696, -158, -231, -185, -497, -810, -601, -997, -900};
        int expected = -17;
        int result = find_max_array(static_cast<int>(base.size()), base.data());
        REQUIRE(result == expected);
    }
    SECTION("Секция #3")
    {
        std::vector<int> base{5, 5, 5, 5, 5};
        int expected = 5;
        int result = find_max_array(static_cast<int>(base.size()), base.data());
        REQUIRE(result == expected);
    }
}

TEST_CASE("TEST F6: верно ли, что среди элементов массива есть два одинаковых.")
{
    SECTION("Секция #1")
    {
        std::vector<int> base{
            1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
            35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
            52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
            69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85,
            86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
        int expected = 0;
        int result = is_two_same(static_cast<int>(base.size()), base.data());
        REQUIRE(result == expected);
    }
    SECTION("Секция #2")
    {
        std::vector<int> base{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 3};
        int expected = 1;
        int result = is_two_same(static_cast<int>(base.size()), base.data());
        REQUIRE(result == expected);
    }
}

TEST_CASE("TEST F7: Функция, которая сжимает серии массива")
{
    SECTION("Секция #1")
    {
        std::vector<int> base{1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0,
                              1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1,
                              1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0,
                              1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1,
                              0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1,
                              0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0};

        std::vector<int> expected{0, 2, 2, 3, 1, 2, 2, 1, 1, 1, 2, 1, 4, 4,
                                1, 1, 4, 3, 1, 2, 1, 1, 3, 2, 2, 1, 3, 1,
                                2, 3, 1, 1, 6, 1, 1, 1, 2, 2, 1, 1, 4, 1,
                                2, 2, 1, 1, 3, 1, 1, 6, 1, 1, 2};

        std::vector<int> result(expected.size());
                                
        int size_expected = 53;
        int size = compression(base.data(), result.data(), static_cast<int>(base.size()));
        REQUIRE(result == expected);
        REQUIRE(size == size_expected);
    }
}

TEST_CASE("TEST F8: Определить пропущенное число")
{
    SECTION("Секция #1.1: Тест ф-ции заполнения массива")
    {
        constexpr int size = 9;
        constexpr int maxsize = 1000;
        TestScanf ts("2   3   1   4   7   6   9   8  10   0");
        std::array<int, size> arr{};
        std::array<int, size> expected{2, 3, 1, 4, 7, 6, 9, 8, 10};
        int sz = fill_array_F8(arr.data(), maxsize);
        REQUIRE(size == sz);
        REQUIRE(arr == expected);
    }
    SECTION("Секция #1.2: Тест ф-ции заполнения массива")
    {
        constexpr int size = 0;
        constexpr int maxsize = 1000;
        TestScanf ts("");
        std::array<int, size> arr{};
        std::array<int, size> expected{};
        int sz = fill_array_F8(arr.data(), maxsize);
        REQUIRE(size == sz);
        REQUIRE(arr == expected);
    }
}

TEST_CASE("TEST F9: функция в массиве находит максимальный из отрицательных элементов и меняет его местами с последним элементом")
{
    SECTION("Секция #1")
    {
        std::vector<int> base{1, -2, -3, -4, 5, 6, 7, 8, 9, 10};
        std::vector<int> expected{1, 10, -3, -4, 5, 6, 7, 8, 9, -2};
        swap_negmax_last(static_cast<int>(base.size()), base.data());
        REQUIRE(base == expected);
    }
}
