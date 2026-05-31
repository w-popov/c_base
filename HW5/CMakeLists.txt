#[[
    ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
]]

cmake_minimum_required(VERSION 3.28)

project(HW5 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Включить тесты
enable_testing()

set(C_SRC_PATH ${CMAKE_SOURCE_DIR}/src)
set(TESTS_SRC ${CMAKE_SOURCE_DIR}/tests)

# Исходники с Д.З .c
add_subdirectory(${C_SRC_PATH})

# Создать библиотеку из catch.hpp
add_library(catchlib STATIC tests/catch_include.cpp)
# Где искать catch.hpp
target_include_directories(catchlib PUBLIC ${TESTS_SRC}/includes/)

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
            ${C_SRC_PATH}/B12.c
            ${C_SRC_PATH}/B13.c
            ${C_SRC_PATH}/B14.c
            ${C_SRC_PATH}/B15.c
            ${C_SRC_PATH}/B16.c
            ${C_SRC_PATH}/B17.c
            ${C_SRC_PATH}/B18.c
            ${C_SRC_PATH}/B19.c
            ${C_SRC_PATH}/B20.c
            ${C_SRC_PATH}/B21.c
            )

# Обрезать main() в исходниках Си
target_compile_definitions(test_hw PRIVATE TEST_DEF_HW5)

# Где искать HW5.h
target_include_directories(test_hw PUBLIC ${C_SRC_PATH}/includes/)

# Создать исполняемый файл тестов catch
add_executable(hwtests ${TESTS_SRC}/test_all.cpp)
target_compile_options(hwtests PRIVATE -g -Wall -Wextra -Wpedantic)

# Создать директорию /runtests если её нет
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/runtests)

# Установить директорию /runtests для исполняемого файла тестов
set_target_properties(hwtests PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/runtests
)

# Прилинковать к исполняемому файлу тестов собранную мной библиотеку catch
target_link_libraries(hwtests PRIVATE catchlib test_hw)

add_test(NAME AllHWTests COMMAND hwtests)
