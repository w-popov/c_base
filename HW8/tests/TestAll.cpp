#include "catch.hpp"
#include <array>
#include "TestScanf.hpp"
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
