#include <cstdio>
#include <fstream>
#include <string>
#include <chrono>
#include <stdexcept>
#include "TestScanf.hpp"

TestScanf::~TestScanf()
{
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
}

TestScanf::TestScanf(std::string text) : _check_text{std::move(text)}
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

std::string TestScanf::_gen_temp_fname () const
{
    static unsigned counter = 0;
    auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    return "test_scanf_" + std::to_string(now) + "_"
            + std::to_string(counter++) + ".txt";
}

