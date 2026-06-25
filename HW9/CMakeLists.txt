#[[
    ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
]]

cmake_minimum_required(VERSION 3.28)

project(HW9 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#[[
  Запуск с Debug:
  
    cmake -B build -DCMAKE_BUILD_TYPE=Debug
    cmake --build build
]]

# Принудительный Release (если не указан явно)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Поддержка Debug конфигурации
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Available build types" FORCE)

# Настройки для Debug конфигурации
set(CMAKE_C_FLAGS_DEBUG "-g -O0" CACHE STRING "C Debug flags" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "C++ Debug flags" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "" CACHE STRING "Debug linker flags" FORCE)

# Загрука FLTK для GUI
# Для поддержки русского языка в linux нужно установить:
# sudo apt install libxft-dev libfontconfig1-dev libx11-dev libxext-dev
# ============================================================
include(FetchContent)
message(STATUS "Detected: Fetching FLTK for graphics")

FetchContent_Declare(
    fltk
    GIT_REPOSITORY https://github.com/fltk/fltk.git
    GIT_TAG release-1.4.5
    GIT_SHALLOW TRUE
)
set(OPTION_BUILD_SHARED_LIBS OFF CACHE BOOL "Static libs" FORCE) 
set(FLTK_BUILD_TEST OFF CACHE BOOL "No tests" FORCE)
set(FLTK_BUILD_EXAMPLES OFF CACHE BOOL "No examples" FORCE)
set(FLTK_BUILD_FLUID OFF CACHE BOOL "No FLUID" FORCE)

FetchContent_MakeAvailable(fltk)

message(STATUS "FLTK fetched: ${fltk_SOURCE_DIR}")
# ============================================================

# Линковка исполняемых файлов с библиотеками для переносимости
if(WIN32)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # Настройки для Debug
        set(CMAKE_C_FLAGS_DEBUG "-g -O0")
        set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
        
        # Статическая линковка (для Release)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
        
        # Оптимизация размера (для Release)
        set(CMAKE_C_FLAGS_RELEASE "-O2 -s")
        set(CMAKE_CXX_FLAGS_RELEASE "-O2 -s")
        # Флаг -flto только при линковке EXE:
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto -s")
        
        # Для Debug отключение оптимизации и LTO
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "")
    else()
        # MSVC (cl.exe)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
        
        # Release настройки для MSVC
        set(CMAKE_CXX_FLAGS_RELEASE "/O2 /GL /GS- /Gy")
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /OPT:REF /OPT:ICF")
        
        # Debug настройки для MSVC
        set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /RTC1")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/DEBUG")
    endif()
endif()

# Включить тесты
enable_testing()

set(C_SRC_PATH ${CMAKE_SOURCE_DIR}/src)
set(TESTS_SRC_PATH ${CMAKE_SOURCE_DIR}/tests)
set(GUI_SRC_PATH ${CMAKE_SOURCE_DIR}/gui)

set(SOURCES_C ${C_SRC_PATH}/F1.c
              ${C_SRC_PATH}/F2.c
              ${C_SRC_PATH}/F3.c
              ${C_SRC_PATH}/F4.c
              ${C_SRC_PATH}/F5.c
              ${C_SRC_PATH}/F6.c
              ${C_SRC_PATH}/F7.c
              ${C_SRC_PATH}/F8.c
              ${C_SRC_PATH}/F9.c
              ${C_SRC_PATH}/F10.c
              ${C_SRC_PATH}/F11.c
              ${C_SRC_PATH}/F12.c
              ${C_SRC_PATH}/F13.c
              ${C_SRC_PATH}/F14.c
              ${C_SRC_PATH}/F15.c
              ${C_SRC_PATH}/F16.c
              ${C_SRC_PATH}/F17.c
              ${C_SRC_PATH}/F18.c
              ${C_SRC_PATH}/F19.c
              ${C_SRC_PATH}/F20.c
              )

set(SOURCES_GUI ${GUI_SRC_PATH}/gF1.cpp
                ${GUI_SRC_PATH}/gF2.cpp
                # ${GUI_SRC_PATH}/gF3.cpp
                # ${GUI_SRC_PATH}/gF4.cpp
                # ${GUI_SRC_PATH}/gF5.cpp
                # ${GUI_SRC_PATH}/gF6.cpp
                # ${GUI_SRC_PATH}/gF7.cpp
                ${GUI_SRC_PATH}/gF8.cpp
                # ${GUI_SRC_PATH}/gF9.cpp
                # ${GUI_SRC_PATH}/gF10.cpp
                # ${GUI_SRC_PATH}/gF11.cpp
                # ${GUI_SRC_PATH}/gF12.cpp
                # ${GUI_SRC_PATH}/gF13.cpp
                # ${GUI_SRC_PATH}/gF14.cpp
                # ${GUI_SRC_PATH}/gF15.cpp
                # ${GUI_SRC_PATH}/gF16.cpp
                # ${GUI_SRC_PATH}/gF17.cpp
                # ${GUI_SRC_PATH}/gF18.cpp
                # ${GUI_SRC_PATH}/gF19.cpp
                # ${GUI_SRC_PATH}/gF20.cpp
              )

add_subdirectory(${C_SRC_PATH})
add_subdirectory(${TESTS_SRC_PATH})
add_subdirectory(${GUI_SRC_PATH})

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")