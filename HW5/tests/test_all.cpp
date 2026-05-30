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

    int calculate_even_odd_digits_number(unsigned* even, 
                                        unsigned* odd, 
                                        unsigned _number);  // B13

    unsigned number_of_even_numbers (char*, const int);     // B15
    
    struct UniquePtr_u { unsigned* u_ptr; };
    void auto_free(struct UniquePtr_u*);
    #define UNIQE_PTR_U __attribute__((cleanup(auto_free))) struct UniquePtr_u
    struct UniquePtr_u happy_numbers (unsigned);            // B17

    const char* to_lowercase (char* text, char* lower);     // B21
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

SCENARIO("TEST B13: Посчитать количество четных и нечетных цифр числа ")
{
    WHEN("#13-1: number=1234")
    {
        unsigned number = 1234, even = 0, odd = 0;
        THEN("#12-1: even=2, odd=2")
        {
            REQUIRE( calculate_even_odd_digits_number(&even, &odd, number) == 0 );
            REQUIRE( even == 2 );
            REQUIRE( odd == 2 );
        }
    }
    WHEN("#13-2: number=0")
    {
        unsigned number = 0, even = 0, odd = 0;
        THEN("#12-2: even=1, odd=0")
        {
            REQUIRE( calculate_even_odd_digits_number(&even, &odd, number) == 0 );
            REQUIRE( even == 1 );
            REQUIRE( odd == 0 );
        }
    }
    WHEN("#13-3: number=1")
    {
        unsigned number = 1, even = 0, odd = 0;
        THEN("#12-3: even=0, odd=1")
        {
            REQUIRE( calculate_even_odd_digits_number(&even, &odd, number) == 0 );
            REQUIRE( even == 0 );
            REQUIRE( odd == 1 );
        }
    }
    WHEN("#13-4: number=92037165")
    {
        unsigned number = 92037165, even = 0, odd = 0;
        THEN("#12-3: even=3, odd=5")
        {
            REQUIRE( calculate_even_odd_digits_number(&even, &odd, number) == 0 );
            REQUIRE( even == 3 );
            REQUIRE( odd == 5 );
        }
    }
     WHEN("#13-5: number=11")
    {
        unsigned number = 11, even = 0, odd = 0;
        THEN("#12-5: even=0, odd=2")
        {
            REQUIRE( calculate_even_odd_digits_number(&even, &odd, number) == 0 );
            REQUIRE( even == 0 );
            REQUIRE( odd == 2 );
        }
    }
     WHEN("#13-6: number=82")
    {
        unsigned number = 82, even = 0, odd = 0;
        THEN("#13-6: even=2, odd=0")
        {
            REQUIRE( calculate_even_odd_digits_number(&even, &odd, number) == 0 );
            REQUIRE( even == 2 );
            REQUIRE( odd == 0 );
        }
    }
    GIVEN("#13 Нулевой указатель(и)") // *************
    {
        WHEN("#13-1: odd=NULL")
        {
            unsigned number = 123, even = 0, *odd = NULL;
            THEN("#13-1: even=0 ")
            {
                REQUIRE( calculate_even_odd_digits_number(&even, odd, number) == 1 );
                REQUIRE( even == 0 );
            }
        }
        WHEN("#13-2: even=NULL")
        {
            unsigned number = 123, *even = NULL, odd = 0;
            THEN("#13-2: odd=0")
            {
                REQUIRE( calculate_even_odd_digits_number(even, &odd, number) == 1 );
                REQUIRE( odd == 0 );
            }
        }
        WHEN("#13-3: number=123, even=NULL, odd=NULL")
        {
            unsigned number = 123, *even = NULL, *odd = NULL;
            THEN("#13-3: return error")
            {
                REQUIRE( calculate_even_odd_digits_number(even, odd, number) == 1 );
            }
        }
    }
}

TEST_CASE( "TEST B15: Посчитать количество четных чисел" )
{
    constexpr int size = 128;
    char s1[size] = {"0\n"};
    char s2[size] = {"1 0\n"};
    char s3[size] = {" 2 0\n"};
    char s4[size] = {"1 2 3 4 0\n"};
    char s5[size] = {"2 4 6 5 0\n"};
    char s6[size] = {"1  1 1 1  9 0 \n"};

    SECTION("B15 Секция #1")
    {
        REQUIRE( number_of_even_numbers(s1, size) == 0 );
        REQUIRE( number_of_even_numbers(s2, size) == 0 );
        REQUIRE( number_of_even_numbers(s3, size) == 1 );
        REQUIRE( number_of_even_numbers(s4, size) == 2 );
        REQUIRE( number_of_even_numbers(s5, size) == 3 );
        REQUIRE( number_of_even_numbers(s6, size) == 0 );
    }
}

TEST_CASE( "TEST B17: все числа от 10 до введенного числа - у которых сумма цифр равна произведению цифр " )
{    
    unsigned number = 0;
    std::vector<unsigned> x;
    
    SECTION( "B17 Секция #1:" )
    {
        number = 200;
        std::vector<unsigned> v {22, 123, 132};
        x.reserve(v.size());
        UNIQE_PTR_U uptr = happy_numbers(number);
        for (uint i = 0; i < v.size(); ++i) {
            if (uptr.u_ptr)
                x.push_back(uptr.u_ptr[i]);
        }
        REQUIRE_THAT( x, Catch::Matchers::Equals(v) );   
    }
    SECTION( "B17 Секция #2:" )
    {
        x.clear();
        number = 1000;
        std::vector<unsigned> v {22, 123, 132, 213, 231, 312, 321};
        x.reserve(v.size());
        UNIQE_PTR_U uptr = happy_numbers(number);
        for (uint i = 0; i < v.size(); ++i) {
            if (uptr.u_ptr)
                x.push_back(uptr.u_ptr[i]);
        } 
        REQUIRE_THAT( x, Catch::Matchers::Equals(v) );
    }
}

TEST_CASE( "TEST B21: Перевести все заглавные английские буквы в строчные" )
{
    constexpr int size = 512;
    char u1[size] = "HELlo WoRlD";
    char u2[size] = "  AGd * dO-R+W  ";
    char u3[size] = "small letters";
    char l1[size];
    char l2[size];
    char l3[size];
    
    SECTION("B21 Секция #1")
    {
        REQUIRE( std::string(to_lowercase(u1, l1)) ==  std::string("hello world"));
        REQUIRE( std::string(to_lowercase(u2, l2)) ==  std::string("  agd * do-r+w  "));
        REQUIRE( std::string(to_lowercase(u3, l3)) ==  std::string("small letters"));
    }
}

