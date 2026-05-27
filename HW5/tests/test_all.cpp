#include "catch.hpp"
#include <ranges>

extern "C"
{
    int square_power(int number);   // B2
}


TEST_CASE( "TEST B2: Квадрат числа " )
{    
    int min = 1, max = 3;
    std::vector<int> x;
    SECTION( "B2 Секция #1: 1 <=> 3" )
    {
        std::vector<int> v {1, 4, 9};
        x.reserve(v.size());
        for (int i : std::views::iota(min, max + 1))
        {
            x.push_back(square_power(i));
        }
        REQUIRE_THAT( x, Catch::Matchers::Equals(v) );   
    }

    SECTION( "B2 Секция #2: 1 <=> 5" )
    {
        x.clear();
        min = 1; 
        max = 5;
        std::vector<int> v {1, 4, 9, 16, 25};
        x.reserve(v.size());
        for (int i : std::views::iota(min, max + 1))
        {
            x.push_back(square_power(i));
        }
        REQUIRE_THAT( x, Catch::Matchers::Equals(v) );
    }
}