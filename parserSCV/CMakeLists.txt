project(parcerCSV LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(C_SRC_PATH ${CMAKE_SOURCE_DIR}/src)
set(SOURCES_C ${C_SRC_PATH}/main.c ${C_SRC_PATH}/parcer_csv.c)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
set(TARGET_NAME main)

include_directories(${C_SRC_PATH}/headers)
add_executable(${TARGET_NAME} ${C_SRC_PATH}/main.c ${C_SRC_PATH}/parcer_csv.c)
target_compile_options(${TARGET_NAME} PRIVATE -g -Wall -Wextra -Wpedantic)
target_link_libraries(${TARGET_NAME} PRIVATE m)
