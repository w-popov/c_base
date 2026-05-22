#define CATCH_CONFIG_MAIN
#include "catch.hpp"

extern "C"
{
    int find_max(int a, int b, int c);
    int min_from_5_numbers(int, int, int, int, int);
    int max_from_5_numbers(int, int, int, int, int);
    int sum_min_max_from_5_numbers(int, int, int, int, int);
}

TEST_CASE( "TEST A8: find_max returns the maximum of three integers" )
{
    REQUIRE( find_max(1, 2, 3) == 3 );
    REQUIRE( find_max(3, 2, 1) == 3 );
    REQUIRE( find_max(1, 3, 2) == 3 );
    REQUIRE( find_max(-1, -2, -3) == -1 );
    REQUIRE( find_max(-1000, 0, 3135) == 3135 );
    
    SECTION( "TEST A8: find_max handles equal values correctly" )
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
TEST_CASE( "TEST A9: max_from_5_numbers returns the maximum of five integers" )
{    
    SECTION( "TEST A9: find_max handles equal values correctly" )
    {
        REQUIRE( max_from_5_numbers(5, 5, 3, 0, 789643) == 789643 );
        REQUIRE( max_from_5_numbers(5, 5, 5, 5, 5) == 5 );
        REQUIRE( max_from_5_numbers(-5, -5, -5, -5, -5) == -5 );
        REQUIRE( max_from_5_numbers(0, 0, 0, 0, 0) == 0 );
        REQUIRE( max_from_5_numbers(50, 500, 500, 50, 500) == 500 );
    }
}
TEST_CASE( "TEST A10: min_from_5_numbers returns the minimum of five integers" )
{    
    SECTION( "TEST A10: min_from_5_numbers handles equal values correctly" )
    {
        REQUIRE( min_from_5_numbers(5, 5, 3, 0, 789643) == 0 );
        REQUIRE( min_from_5_numbers(5, 5, 5, 5, 5) == 5 );
        REQUIRE( min_from_5_numbers(-5, -5, -5, -5, -5) == -5 );
        REQUIRE( min_from_5_numbers(0, 0, 0, 0, 0) == 0 );
        REQUIRE( min_from_5_numbers(50, 500, -500, 50, 500) == -500 );
    }
}
TEST_CASE( "TEST A11: returns the summ maximum and minimum of five integers" )
{    
    SECTION( "min + max handles equal values correctly" )
    {
        REQUIRE( sum_min_max_from_5_numbers(5, 5, 3, 0, 789643) == 789643 );
        REQUIRE( sum_min_max_from_5_numbers(5, 5, 5, 5, 5) == 10 );
        REQUIRE( sum_min_max_from_5_numbers(-5, -5, -5, -5, -5) == -10 );
        REQUIRE( sum_min_max_from_5_numbers(0, 0, 0, 0, 0) == 0 );
        REQUIRE( sum_min_max_from_5_numbers(50, 500, -500, 50, 500) == 0 );
        
    }
}