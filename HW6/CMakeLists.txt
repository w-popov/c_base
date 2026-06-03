#[[
    ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
]]

cmake_minimum_required(VERSION 3.28)

project(HW6 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(C_SRC_PATH ${CMAKE_SOURCE_DIR}/src)

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

