#[[
    Основной CMakeLists.txt выше.
]]

SET(C_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# include "HW6.h"
include_directories(${C_SRC_PATH}/includes)

# Задать директорию для исполняемых файлов в корневой bin/
# CMAKE_RUNTIME_OUTPUT_DIRECTORY - переменная CMake
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Создать директорию bin если её нет
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

foreach(C_FILE ${SOURCES_C})
    get_filename_component(TARGET_NAME ${C_FILE} NAME_WE)
    add_executable(${TARGET_NAME} ${C_FILE})
    target_compile_options(${TARGET_NAME} PRIVATE -g -Wall -Wextra -Wpedantic)
endforeach()
