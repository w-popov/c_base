#include <cstdio>
#include <fstream>
#include <string>
#include <chrono>
#include <stdexcept>
#include "catch.hpp"
#include "HW7.h"

#ifdef _WIN32
#include <io.h>
#define TEST_DUP _dup
#define TEST_DUP2 _dup2
#define TEST_FILENO _fileno
#define TEST_CLOSE _close
#else
#include <unistd.h>
#define TEST_DUP dup
#define TEST_DUP2 dup2
#define TEST_FILENO fileno
#define TEST_CLOSE close
#endif

/**
 * Класс для тестирования си функций вызывающих scanf();
 * Принимает строку как будто через входной поток.
 * Пример #1:
 *   SECTION("D19: Секция #1")
 *   {
 *       TestScanf tscf("abcd a.");
 *       REQUIRE(acounter() == 2);
 *   }
 * Пример #2:
 *   int sum_two_numbers(void)
 *   {
 *       int a, b;
 *       scanf("%d %d", &a, &b);
 *       return a + b;
 *   }
 *   SECTION("Сумма 2 чисел")
 *   {
 *       TestScanf tscf("10 20");
 *       REQUIRE(sum_two_numbers() == 30);
 *   }
 */
class TestScanf
{
  private:
    int _stdin_backup_fd = -1;
    std::string _temp_fname = "";
    std::string _check_text = "";

  public:
    ~TestScanf()
    {
        if (_stdin_backup_fd != -1)
        {
            TEST_DUP2(_stdin_backup_fd, TEST_FILENO(stdin));
            TEST_CLOSE(_stdin_backup_fd);
        }
        if (!_temp_fname.empty())
        {
            std::remove(_temp_fname.c_str());
        }
    }

    TestScanf() = delete;
    TestScanf(const TestScanf &) = delete;
    TestScanf(TestScanf &&) noexcept = delete;
    TestScanf &operator=(const TestScanf &) = delete;
    TestScanf &operator=(TestScanf &&) noexcept = delete;

    explicit TestScanf(std::string text) : _check_text{std::move(text)}
    {
        _temp_fname = _gen_temp_fname();
        std::ofstream test_scanf_file(_temp_fname);
        if (!test_scanf_file)
        {
            throw std::runtime_error("Временный файл " + _temp_fname +
                                     " не создан!\n");
        }
        test_scanf_file << _check_text;
        test_scanf_file.close();

        _stdin_backup_fd = TEST_DUP(TEST_FILENO(stdin));
        if (_stdin_backup_fd == -1)
        {
            throw std::runtime_error("Ошибка сохранения stdin");
        }
        if (!freopen(_temp_fname.c_str(), "r", stdin))
        {
            TEST_CLOSE(_stdin_backup_fd);
            _stdin_backup_fd = -1;
            throw std::runtime_error("Ошибка перенаправления stdin");
        }
    }

  private:
    std::string _gen_temp_fname () const
    {
        static unsigned counter = 0;
        auto now = std::chrono::steady_clock::now().time_since_epoch().count();
        return "test_scanf_" + std::to_string(now) + "_"
               + std::to_string(counter++) + ".txt";
    }
};

TEST_CASE("TEST D2: рекурсивно определяет сумму всех чисел от 1 до N")
{
    SECTION("Секция #1")
    {
        REQUIRE(sum_range_numbers_rec(0) == 0);
        REQUIRE(sum_range_numbers_rec(1) == 1);
        REQUIRE(sum_range_numbers_rec(3) == 6);
        REQUIRE(sum_range_numbers_rec(5) == 15);
    }
}

TEST_CASE("TEST D9: рекурсивно вычислить сумму цифр N")
{
    SECTION("Секция #1")
    {
        REQUIRE(sum_digits(0) == 0);
        REQUIRE(sum_digits(1) == 1);
        REQUIRE(sum_digits(3) == 3);
        REQUIRE(sum_digits(15) == 6);
        REQUIRE(sum_digits(10) == 1);
        REQUIRE(sum_digits(1111) == 4);
        REQUIRE(sum_digits(10000000) == 1);
    }
}

TEST_CASE("TEST D10: рекурсивно проверить число на простоту")
{
    SECTION("Секция #1")
    {
        REQUIRE(is_prime(0, 2) == 0);
        REQUIRE(is_prime(1, 2) == 0);
        REQUIRE(is_prime(2, 2) == 1);
        REQUIRE(is_prime(3, 2) == 1);
        REQUIRE(is_prime(11, 2) == 1);
        REQUIRE(is_prime(12, 2) == 0);
        REQUIRE(is_prime(73, 2) == 1);
    }
}

TEST_CASE("TEST D16: рекурсивно определяет является ли введенное натуральное "
          "число точной степенью двойки")
{
    SECTION("Секция #1")
    {
        REQUIRE(is2pow(0) == 0);
        REQUIRE(is2pow(1) == 1);
        REQUIRE(is2pow(2) == 1);
        REQUIRE(is2pow(3) == 0);
        REQUIRE(is2pow(4) == 1);
        REQUIRE(is2pow(100) == 0);
        REQUIRE(is2pow(32768) == 1);
    }
}

TEST_CASE("TEST D17: рекурсивно проверить функцию Аккермана")
{
    /*
       m\n |  0  |  1  |  2  |  3  |  4  |
        0  |   1 |   2 |   3 |   4 |   5 |
        1  |   2 |   3 |   4 |   5 |   6 |
        2  |   3 |   5 |   7 |   9 |  11 |
        3  |   5 |  13 |  29 |  61 | 125 |
    */
    SECTION("Секция #1")
    {
        REQUIRE(akkerman(0, 0) == 1);
        REQUIRE(akkerman(1, 1) == 3);
        REQUIRE(akkerman(2, 1) == 5);
        REQUIRE(akkerman(2, 4) == 11);
        REQUIRE(akkerman(3, 4) == 125);
    }
}

TEST_CASE("TEST D19: считывает данную строку со стандартного потока ввода и \
    возвращает целое число - количество символов 'a' в данной строке")
{
    SECTION("Секция #1")
    {
        TestScanf tscf("abcd a.");
        REQUIRE(acounter() == 2);
    }
    SECTION("Секция #2")
    {
        TestScanf tscf("bcd .");
        REQUIRE(acounter() == 0);
    }
    SECTION("Секция #3")
    {
        TestScanf tscf("aaaaaaaaaa.");
        REQUIRE(acounter() == 10);
    }
    SECTION("Секция #4")
    {
        TestScanf tscf(".");
        REQUIRE(acounter() == 0);
    }
}

TEST_CASE(
    "TEST D20: рекурсивно определяет возведения целого числа n в степень p")
{
    SECTION("Секция #1")
    {
        REQUIRE(recurs_power(1, 0) == 1);
        REQUIRE(recurs_power(0, 0) == 1);
        REQUIRE(recurs_power(-1, 1) == -1);
        REQUIRE(recurs_power(-2, 2) == 4);
        REQUIRE(recurs_power(3, 3) == 27);
        REQUIRE(recurs_power(-3, 3) == -27);
    }
}
