#include "catch.hpp"
#include <ranges>

extern "C"
{
    int square_power(int number);                           // B2
    int summ_squares_range(int a, int b);                   // B3
    const char* is_it_number_have_3_digit(unsigned number); // B4
    unsigned summ_digits_of_number(unsigned number);        // B5
    const char* is_both_equal_digits_in_number(int number); // B6
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

TEST_CASE( "TEST B3: Сумма квадратов чисел от a до b " )
{
    SECTION("B3 Секция #1")
    {
        REQUIRE( summ_squares_range(4, 10) == 371 );
        REQUIRE( summ_squares_range(1, 5) == 55 );
        REQUIRE( summ_squares_range(1, 1) == 1 );
        REQUIRE( summ_squares_range(5, 5) == 25 );
    }
}

TEST_CASE( "TEST B4: В числе ровно 3 цифры " )
{
    std::string Y = "YES", N = "NO";
    auto match = [](unsigned num)->std::string {
        return static_cast<std::string>(is_it_number_have_3_digit(num));
    };
    SECTION("B4 Секция #1")
    {
        REQUIRE( match(100) == Y );
        REQUIRE( match(123) == Y );
        REQUIRE( match(10) == N );
        REQUIRE( match(1) == N );
        REQUIRE( match(1000) == N );
        REQUIRE( match(1998121) == N );
    }
}

TEST_CASE( "TEST B5: Сумма цифр числа " )
{
    SECTION("B5 Секция #1")
    {
        REQUIRE( summ_digits_of_number(1234) == 10 );
        REQUIRE( summ_digits_of_number(123) == 6 );
        REQUIRE( summ_digits_of_number(12) == 3 );
        REQUIRE( summ_digits_of_number(1) == 1 );
        REQUIRE( summ_digits_of_number(0) == 0 );
    }
}

TEST_CASE( "TEST B6: верно ли, что в записи числа есть две одинаковые цифры, стоящие рядом " )
{
    std::string Y = "YES", N = "NO";
    auto match = [](int num)->std::string {
        return static_cast<std::string>(is_both_equal_digits_in_number(num));
    };
    SECTION("B6 Секция #1")
    {
        REQUIRE( match(100) == Y );
        REQUIRE( match(1223) == Y );
        REQUIRE( match(10) == N );
        REQUIRE( match(1) == N );
        REQUIRE( match(1000) == Y );
        REQUIRE( match(199811) == Y );
        REQUIRE( match(0) == N );
        REQUIRE( match(11) == Y );
        REQUIRE( match(11111) == Y );
    }
}