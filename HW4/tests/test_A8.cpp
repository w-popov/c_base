
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

extern "C"
{
    int find_max(int a, int b, int c);
}

TEST_CASE( "find_max returns the maximum of three integers" )
{
    REQUIRE( find_max(1, 2, 3) == 3 );
    REQUIRE( find_max(3, 2, 1) == 3 );
    REQUIRE( find_max(1, 3, 2) == 3 );
    REQUIRE( find_max(-1, -2, -3) == -1 );
    REQUIRE( find_max(-1000, 0, 3135) == 3135 );
    
    SECTION( "find_max handles equal values correctly" )
    {
        REQUIRE( find_max(5, 5, 3) == 5 );
        REQUIRE( find_max(5, 3, 5) == 5 );
        REQUIRE( find_max(3, 5, 5) == 5 );
        REQUIRE( find_max(5, 5, 5) == 5 );
        REQUIRE( find_max(0, 0, 0) == 0 );
        REQUIRE( find_max(-3, -3, -3) == -3 );
        REQUIRE( find_max(3, -3, -3) == 3 );
        REQUIRE( find_max(-3, 3, -3) == 3 );
        REQUIRE( find_max(-3, -3, 3) == 3 );
    }
}
