#[[
    Основной CMakeLists.txt выше.
]]

SET(C_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# include "HW5.h"
include_directories(${C_SRC_PATH}/includes)

# Задать директорию для исполняемых файлов в корневой bin/
# CMAKE_RUNTIME_OUTPUT_DIRECTORY - переменная CMake
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Создать директорию bin если её нет
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_executable(B1 ${C_SRC_PATH}/B1.c)
target_compile_options(B1 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B2 ${C_SRC_PATH}/B2.c)
target_compile_options(B2 PRIVATE -g -Wall -Wextra -Wpedantic)
# target_include_directories(B2 PUBLIC ${C_SRC_PATH}/includes)

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

add_executable(B12 ${C_SRC_PATH}/B12.c)
target_compile_options(B12 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B13 ${C_SRC_PATH}/B13.c)
target_compile_options(B13 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B14 ${C_SRC_PATH}/B14.c)
target_compile_options(B14 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B15 ${C_SRC_PATH}/B15.c)
target_compile_options(B15 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B16 ${C_SRC_PATH}/B16.c)
target_compile_options(B16 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B17 ${C_SRC_PATH}/B17.c)
target_compile_options(B17 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B18 ${C_SRC_PATH}/B18.c)
target_compile_options(B18 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B19 ${C_SRC_PATH}/B19.c)
target_compile_options(B19 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B20 ${C_SRC_PATH}/B20.c)
target_compile_options(B20 PRIVATE -g -Wall -Wextra -Wpedantic)

add_executable(B21 ${C_SRC_PATH}/B21.c)
target_compile_options(B21 PRIVATE -g -Wall -Wextra -Wpedantic)