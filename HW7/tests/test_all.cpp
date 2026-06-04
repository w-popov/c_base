#include "catch.hpp"
#include "HW7.h"


TEST_CASE( "TEST D2: рекурсивно определяет сумму всех чисел от 1 до N" )
{
    SECTION("Секция #1")
    {
        REQUIRE( sum_range_numbers_rec(0) == 0 );
        REQUIRE( sum_range_numbers_rec(1) == 1 );
        REQUIRE( sum_range_numbers_rec(3) == 6 );
        REQUIRE( sum_range_numbers_rec(5) == 15 );
    }
}

