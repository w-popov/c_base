#[[ HW8/src: Компиляция исходников си контеста Е ]]

# Где искать .h
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/includes)

# Задать директорию для исполняемых файлов HW8/bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Создать директорию HW8/bin если её нет
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#[[ 
    Компиляция из списка исходников SOURCES_C.
    SOURCES_C определен в корневом CMakeLists.txt
]]
foreach(C_FILE ${SOURCES_C})
    get_filename_component(TARGET_NAME ${C_FILE} NAME_WE)
    add_executable(${TARGET_NAME} ${C_FILE})
    target_compile_options(${TARGET_NAME} PRIVATE -g -Wall -Wextra -Wpedantic)
    target_link_libraries(${TARGET_NAME} PRIVATE m)    
endforeach()
