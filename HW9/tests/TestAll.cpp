#include "catch.hpp"
#include <array>
#include <algorithm>
#include "TestScanf.hpp"
#include "func_tools.hpp"
#include "HW9.h"

TEST_CASE("TEST F1: test")
{
    SECTION("Секция #1")
    {
        REQUIRE(test_add(1, 2) == 3);
    }
}


