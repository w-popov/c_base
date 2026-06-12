#include "catch.hpp"
#include "TestScanf.hpp"
#include "HW8.h"

TEST_CASE("TEST E: test")
{
    SECTION("Секция #1")
    {
        TestScanf tscf("12 10");
        REQUIRE(test_sum() == 22);
    }
}

