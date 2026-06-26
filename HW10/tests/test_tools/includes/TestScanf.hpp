#pragma once

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
    int _stdin_backup_fd{-1};
    std::string _temp_fname{""};
    std::string _check_text{""};

  public:
    ~TestScanf();

    TestScanf() = delete;
    TestScanf(const TestScanf &) = delete;
    TestScanf(TestScanf &&) noexcept = delete;
    TestScanf &operator=(const TestScanf &) = delete;
    TestScanf &operator=(TestScanf &&) noexcept = delete;

    explicit TestScanf(std::string text);

  private:
    std::string _gen_temp_fname () const;
};