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
        std::vector<int> base{20, 19, 4, 3, 2, 1, 18, 17, 13, 12, 11, 16, 15, 14, 10, 9, 8, 7, 6, 5};
        std::vector<int> expected{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
        sort_array(static_cast<int>(base.size()), base.data());
        REQUIRE(base == expected); 
    }
}

TEST_CASE("TEST F2: начало массива все четные элементы, а в конец – все нечетные")
{
    SECTION("Секция #1")
    {
        std::vector<int> base{20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
        std::vector<int> expected{20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1};
        sort_even_odd(static_cast<int>(base.size()), base.data());
        REQUIRE(base == expected); 
    }
}
