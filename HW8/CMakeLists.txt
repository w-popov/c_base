#[[
    ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
]]

cmake_minimum_required(VERSION 3.28)

project(HW8 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
set(TESTS_SRC_PATH ${CMAKE_SOURCE_DIR}/tests)

set(SOURCES_C ${C_SRC_PATH}/E1.c
              ${C_SRC_PATH}/E2.c
              ${C_SRC_PATH}/E3.c
              ${C_SRC_PATH}/E4.c
              ${C_SRC_PATH}/E5.c
              ${C_SRC_PATH}/E6.c
              ${C_SRC_PATH}/E7.c
              ${C_SRC_PATH}/E8.c
              ${C_SRC_PATH}/E9.c
              ${C_SRC_PATH}/E10.c
              ${C_SRC_PATH}/E11.c
              ${C_SRC_PATH}/E12.c
              ${C_SRC_PATH}/E13.c
              ${C_SRC_PATH}/E14.c
              ${C_SRC_PATH}/E15.c
              ${C_SRC_PATH}/E16.c
              ${C_SRC_PATH}/E17.c
              ${C_SRC_PATH}/E18.c
              ${C_SRC_PATH}/E19.c
              ${C_SRC_PATH}/E20.c
              )

add_subdirectory(${C_SRC_PATH})
add_subdirectory(${TESTS_SRC_PATH})