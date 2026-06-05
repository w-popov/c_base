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

TEST_CASE( "TEST D9: рекурсивно вычислить сумму цифр N" )
{
    SECTION("Секция #1")
    {
        REQUIRE( sum_digits(0) == 0 );
        REQUIRE( sum_digits(1) == 1 );
        REQUIRE( sum_digits(3) == 3 );
        REQUIRE( sum_digits(15) == 6 );
        REQUIRE( sum_digits(10) == 1 );
        REQUIRE( sum_digits(1111) == 4 );
        REQUIRE( sum_digits(10000000) == 1 );
    }
}

TEST_CASE( "TEST D10: рекурсивно проверить число на простоту" )
{
    SECTION("Секция #1")
    {
        REQUIRE( is_prime(0, 2) == 0 );
        REQUIRE( is_prime(1, 2) == 0 );
        REQUIRE( is_prime(2, 2) == 1 );
        REQUIRE( is_prime(3, 2) == 1 );
        REQUIRE( is_prime(11, 2) == 1 );
        REQUIRE( is_prime(12, 2) == 0 );
        REQUIRE( is_prime(73, 2) == 1 );        
    }
}

