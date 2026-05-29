#include "catch.hpp"
#include <ranges>

extern "C"
{
    int square_power(int number);                           // B2
    int summ_squares_range(int a, int b);                   // B3
    const char* is_it_number_have_3_digit(unsigned number); // B4
    unsigned summ_digits_of_number(unsigned number);        // B5
    const char* is_both_equal_digits_in_number(int number); // B6
    const char* is_equal_digits_in_number(int number);      // B7
    const char* is_nine_digit_in_number(int number);        // B8
    const char* is_even_all_digits_in_number(int number);   // B9
    const char* is_in_ascending_order(int number);          // B10
    unsigned reverse_number(unsigned number);               // B11
    int calculate_min_max_digits_number(unsigned* min, 
                                        unsigned* max, 
                                        unsigned _number);  // B12


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

TEST_CASE( "TEST B7: верно ли, что в записи числа есть две одинаковые цифры, НЕ обязательно стоящие рядом" )
{
    std::string Y = "YES", N = "NO";
    auto match = [](int num)->std::string {
        return static_cast<std::string>(is_equal_digits_in_number(num));
    };
    SECTION("B7 Секция #1")
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
        REQUIRE( match(11511) == Y );
        REQUIRE( match(200002) == Y );
        REQUIRE( match(1494) == Y );
    }
}

TEST_CASE( "TEST B8: верно ли, что в записи числа ровно одна цифра «9»" )
{
    std::string Y = "YES", N = "NO";
    auto match = [](int num)->std::string {
        return static_cast<std::string>(is_nine_digit_in_number(num));
    };
    SECTION("B8 Секция #1")
    {
        REQUIRE( match(100) == N );
        REQUIRE( match(9) == Y );
        REQUIRE( match(1234) == N );
        REQUIRE( match(1009) == Y );
        REQUIRE( match(10099) == N );
        REQUIRE( match(9009) == N );
        REQUIRE( match(19999) == N );
        REQUIRE( match(99999) == N );
    }
}

TEST_CASE( "TEST B9: верно ли, что в записи числа все его цифры четные" )
{
    std::string Y = "YES", N = "NO";
    auto match = [](int num)->std::string {
        return static_cast<std::string>(is_even_all_digits_in_number(num));
    };
    SECTION("B9 Секция #1")
    {
        REQUIRE( match(100) == N );
        REQUIRE( match(9) == N );
        REQUIRE( match(1234) == N );
        REQUIRE( match(1009) == N );
        REQUIRE( match(0) == Y );
        REQUIRE( match(2468) == Y );
        REQUIRE( match(228645) == N );
        REQUIRE( match(6) == Y );        
    }
}

TEST_CASE( "TEST B10: верно ли, что все цифры числа расположены в порядке возрастания" )
{
    std::string Y = "YES", N = "NO";
    auto match = [](int num)->std::string {
        return static_cast<std::string>(is_in_ascending_order(num));
    };
    SECTION("B10 Секция #1")
    {
        REQUIRE( match(100) == N );
        REQUIRE( match(9) == Y );
        REQUIRE( match(1234) == Y );
        REQUIRE( match(1009) == N );
        REQUIRE( match(0) == Y );
        REQUIRE( match(2468) == Y );
        REQUIRE( match(228645) == N );
        REQUIRE( match(68) == Y );
        REQUIRE( match(11) == N );        

    }
}

TEST_CASE( "TEST B11: Реверс цифр числа" )
{
    SECTION("B11 Секция #1")
    {
        REQUIRE( reverse_number(1234) == 4321 );
        REQUIRE( reverse_number(782) == 287 );
        REQUIRE( reverse_number(5) == 5 );
        REQUIRE( reverse_number(468756) == 657864 );
        REQUIRE( reverse_number(23569) == 96532 );
        REQUIRE( reverse_number(1000) == 1 );
        REQUIRE( reverse_number(11111) == 11111 );        
    }
}

SCENARIO("TEST B12: минимальную и максимальную цифры числа number ")
{
    WHEN("#12-1: number=2457")
    {
        unsigned number = 2457, min = UINT32_MAX, max = 0;
        THEN("#12-1: min=2, max=7")
        {
            REQUIRE( calculate_min_max_digits_number(&min, &max, number) == 0 );
            REQUIRE( min == 2 );
            REQUIRE( max == 7 );
        }
    }
    WHEN("#12-2: number=0")
    {
        unsigned number = 0, min = UINT32_MAX, max = 0;
        THEN("#12-2: min=0, max=0")
        {
            REQUIRE( calculate_min_max_digits_number(&min, &max, number) == 0 );
            REQUIRE( min == 0 );
            REQUIRE( max == 0 );
        }
    }
    WHEN("#12-3: number=10")
    {
        unsigned number = 10, min = UINT32_MAX, max = 0;
        THEN("#12-3: min=0, max=1")
        {
            REQUIRE( calculate_min_max_digits_number(&min, &max, number) == 0 );
            REQUIRE( min == 0 );
            REQUIRE( max == 1 );
        }
    }
    WHEN("#12-3: number=92037")
    {
        unsigned number = 92037, min = UINT32_MAX, max = 0;
        THEN("#12-3: min=0, max=9")
        {
            REQUIRE( calculate_min_max_digits_number(&min, &max, number) == 0 );
            REQUIRE( min == 0 );
            REQUIRE( max == 9 );
        }
    }
    GIVEN("#12 Нулевой указатель(и)")
    {
        WHEN("#12-1: number=123, max=NULL")
        {
            unsigned number = 123, min = UINT32_MAX, *max = NULL;
            THEN("#12-1: min=UINT32_MAX, max=NULL")
            {
                REQUIRE( calculate_min_max_digits_number(&min, max, number) == 1 );
                REQUIRE( min == UINT32_MAX );
            }
        }
        WHEN("#12-2: number=123, min=NULL")
        {
            unsigned number = 123, *min = NULL, max = 0;
            THEN("#12-2: min=NULL, max=0")
            {
                REQUIRE( calculate_min_max_digits_number(min, &max, number) == 1 );
                REQUIRE( max == 0 );
            }
        }
        WHEN("#12-3: number=123, max=NULL, min=NULL")
        {
            unsigned number = 123, *min = NULL, *max = NULL;
            THEN("#12-3: min=NULL, max=0")
            {
                REQUIRE( calculate_min_max_digits_number(min, max, number) == 1 );
            }
        }
    }
}



