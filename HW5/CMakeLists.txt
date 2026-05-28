#[[
    ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
]]

cmake_minimum_required(VERSION 3.2)

project(HW5 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Включить тесты
enable_testing()

set(C_SRC_PATH ${CMAKE_SOURCE_DIR}/src)
set(TESTS_SRC ${CMAKE_SOURCE_DIR}/tests)

add_executable(B1 ${C_SRC_PATH}/B1.c)
target_compile_options(B1 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B2 ${C_SRC_PATH}/B2.c)
target_compile_options(B2 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B3 ${C_SRC_PATH}/B3.c)
target_compile_options(B3 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B4 ${C_SRC_PATH}/B4.c)
target_compile_options(B4 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B5 ${C_SRC_PATH}/B5.c)
target_compile_options(B5 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B6 ${C_SRC_PATH}/B6.c)
target_compile_options(B6 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B7 ${C_SRC_PATH}/B7.c)
target_compile_options(B7 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B8 ${C_SRC_PATH}/B8.c)
target_compile_options(B8 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B9 ${C_SRC_PATH}/B9.c)
target_compile_options(B9 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B10 ${C_SRC_PATH}/B10.c)
target_compile_options(B10 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B11 ${C_SRC_PATH}/B11.c)
target_compile_options(B11 PRIVATE -g -Wall -Wextra -Wpedantic)

# add_executable(B12 ${C_SRC_PATH}/B12.c)
# add_executable(B13 ${C_SRC_PATH}/B13.c)
# add_executable(B14 ${C_SRC_PATH}/B14.c)
# add_executable(B15 ${C_SRC_PATH}/B15.c)
# add_executable(B16 ${C_SRC_PATH}/B16.c)
# add_executable(B17 ${C_SRC_PATH}/B17.c)
# add_executable(B18 ${C_SRC_PATH}/B18.c)
# add_executable(B19 ${C_SRC_PATH}/B19.c)
# add_executable(B20 ${C_SRC_PATH}/B20.c)
# add_executable(B21 ${C_SRC_PATH}/B21.c)


# Создать библиотеку из catch.hpp
add_library(catchlib STATIC tests/catch_include.cpp)
# Где искать catch.hpp
target_include_directories(catchlib PUBLIC tests/includes/)

# Создать библиотеку для тестов из исходников Си
# main() в исходниках обрезана препроцессором
add_library(test_hw STATIC 
            ${C_SRC_PATH}/B1.c
            ${C_SRC_PATH}/B2.c
            ${C_SRC_PATH}/B3.c
            ${C_SRC_PATH}/B4.c
            ${C_SRC_PATH}/B5.c
            ${C_SRC_PATH}/B6.c
            ${C_SRC_PATH}/B7.c
            ${C_SRC_PATH}/B8.c
            ${C_SRC_PATH}/B9.c
            ${C_SRC_PATH}/B10.c
            ${C_SRC_PATH}/B11.c
            # ${C_SRC_PATH}/B12.c
            # ${C_SRC_PATH}/B13.c
            # ${C_SRC_PATH}/B14.c
            # ${C_SRC_PATH}/B15.c
            # ${C_SRC_PATH}/B16.c
            # ${C_SRC_PATH}/B17.c
            # ${C_SRC_PATH}/B18.c
            # ${C_SRC_PATH}/B19.c
            # ${C_SRC_PATH}/B20.c
            # ${C_SRC_PATH}/B21.c
            )

# Обрезать main() в исходниках Си
target_compile_definitions(test_hw PRIVATE TEST_DEF_HW5)

# Создать исполняемый файл тестов catch
add_executable(hwtests ${TESTS_SRC}/test_all.cpp)
target_compile_options(hwtests PRIVATE -g -Wall -Wextra -Wpedantic)
# Прилинковать к исполняемому файлу тестов собранную мной библиотеку catch
target_link_libraries(hwtests PRIVATE catchlib test_hw)

add_test(NAME AllHWTests COMMAND hwtests)
