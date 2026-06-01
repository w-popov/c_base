#[[
    ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
]]

cmake_minimum_required(VERSION 3.28)

project(HW5 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Линковка исполняемых файлов с библиотеками для переносимости.
# РАЗМЕР СГЕНЕРИРОВАННЫХ ИСПОЛНЯЕМЫХ ФАЙЛОВ УВЕЛИЧИТСЯ! 
if(WIN32)
    # Проверка, что компилятор — GCC или Clang (а не MSVC)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # Флаги, специфичные для GCC-подобных в Windows (MinGW/Clang)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    else()
        # Блок для MSVC (cl.exe) в Ninja
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
endif()

# Включить тесты
enable_testing()

set(C_SRC_PATH ${CMAKE_SOURCE_DIR}/src)
set(TESTS_SRC ${CMAKE_SOURCE_DIR}/tests)

set(SOURCES_C ${C_SRC_PATH}/C1.c
              ${C_SRC_PATH}/C2.c
              ${C_SRC_PATH}/C3.c
              ${C_SRC_PATH}/C4.c
              ${C_SRC_PATH}/C5.c
              ${C_SRC_PATH}/C6.c
              ${C_SRC_PATH}/C7.c
              ${C_SRC_PATH}/C8.c
              ${C_SRC_PATH}/C9.c
              ${C_SRC_PATH}/C10.c
              ${C_SRC_PATH}/C11.c
              ${C_SRC_PATH}/C12.c
              ${C_SRC_PATH}/C13.c
              ${C_SRC_PATH}/C14.c
              ${C_SRC_PATH}/C15.c
              ${C_SRC_PATH}/C16.c
              ${C_SRC_PATH}/C17.c
              ${C_SRC_PATH}/C18.c
              ${C_SRC_PATH}/C19.c
              ${C_SRC_PATH}/C20.c
              )

# Исходники с Д.З .c
add_subdirectory(${C_SRC_PATH})

# Создать библиотеку из catch.hpp
add_library(catchlib STATIC tests/catch_include.cpp)
# Где искать catch.hpp
target_include_directories(catchlib PUBLIC ${TESTS_SRC}/includes/)

# Создать библиотеку для тестов из исходников Си
# main() в исходниках обрезана препроцессором
add_library(test_hw STATIC ${SOURCES_C})

# Обрезать main() в исходниках Си
target_compile_definitions(test_hw PRIVATE TEST_DEF_HW6)

# Где искать HW6.h
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

# Правильное имя исполняемого файла для разных ОС
if(WIN32)
    set(TEST_EXECUTABLE ${CMAKE_SOURCE_DIR}/runtests/hwtests.exe)
else()
    set(TEST_EXECUTABLE ${CMAKE_SOURCE_DIR}/runtests/hwtests)
endif()

# Добавить тест
add_test(NAME AllHWTests COMMAND ${TEST_EXECUTABLE})
