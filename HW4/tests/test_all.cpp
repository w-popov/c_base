#define CATCH_CONFIG_MAIN
#include "catch.hpp"

extern "C"
{
    int find_max(int a, int b, int c);
    int min_from_5_numbers(int, int, int, int, int);
    int max_from_5_numbers(int, int, int, int, int);
    int sum_min_max_from_5_numbers(int, int, int, int, int);
    unsigned max_digit_from_number(unsigned);
    struct EquationLine { float k; float b; };
    struct EquationLine equation_line(int, int, int, int);
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
TEST_CASE( "TEST A14: returns the maximum digit of input unsigned int number" )
{    
    SECTION( "TEST A14: handles equal values correctly" )
    {
        REQUIRE( max_digit_from_number(345) == 5 );
        REQUIRE( max_digit_from_number(435) == 5 );
        REQUIRE( max_digit_from_number(543) == 5 );
        REQUIRE( max_digit_from_number(453) == 5 );
        REQUIRE( max_digit_from_number(111) == 1 );
        REQUIRE( max_digit_from_number(000) == 0 );
        REQUIRE( max_digit_from_number(122) == 2 );
        REQUIRE( max_digit_from_number(221) == 2 );
        REQUIRE( max_digit_from_number(121) == 2 );
        REQUIRE( max_digit_from_number(112) == 2 );
        REQUIRE( max_digit_from_number(211) == 2 );

    }
}
TEST_CASE( "TEST A15: Determine the equation of a line given the coordinates of two points" )
{   
    SECTION( "TEST A15: #1" )
    {
        struct EquationLine eq = {.k = 0.0, .b = 0.0};
        eq = equation_line(6, 9, -1, 3);
        REQUIRE_THAT(floor(eq.k * 1000) / 1000, Catch::Matchers::WithinAbs(0.86, 0.01f));  
        REQUIRE_THAT(floor(eq.b * 1000) / 1000, Catch::Matchers::WithinAbs(3.86, 0.01f));  

    }
    SECTION( "TEST A15: #2" )
    {
        struct EquationLine eq = {.k = 0.0, .b = 0.0};
        eq = equation_line(1, 2, 3, 4);
        REQUIRE_THAT(floor(eq.k * 1000) / 1000, Catch::Matchers::WithinAbs(1.00, 0.01f));  
        REQUIRE_THAT(floor(eq.b * 1000) / 1000, Catch::Matchers::WithinAbs(1.00, 0.01f));  

    }
}
