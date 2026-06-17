#[[ HW9/tests: Сборка тестов ]]

#[[ 
    Создать библиотеку из catch.hpp
]]
add_library(catch_lib_test STATIC ${CMAKE_CURRENT_SOURCE_DIR}/CatchModule.cpp)

# Где искать catch.hpp
target_include_directories(catch_lib_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/includes)

#[[ 
    Создать библиотеку для тестов из исходников Си
    main() в исходниках обрезана препроцессором
]]
add_library(c_lib_test STATIC ${SOURCES_C})

# Обрезать main() в исходниках Си
target_compile_definitions(c_lib_test PRIVATE TEST_DEF_HW9)

# Где искать HW9.h
target_include_directories(c_lib_test PUBLIC ${C_SRC_PATH}/includes)


#[[
    Создать библиотеку tools_lib_test
]]
set(TOOLS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/test_tools)
add_library(tools_lib_test STATIC ${TOOLS_PATH}/TestScanf.cpp
                                  ${TOOLS_PATH}/func_tools.cpp              
                                )

target_include_directories(tools_lib_test PUBLIC ${TOOLS_PATH}/includes)


# Создать исполняемый файл тестов
add_executable(hwtests ${CMAKE_CURRENT_SOURCE_DIR}/TestAll.cpp)
target_compile_options(hwtests PRIVATE -g -Wall -Wextra -Wpedantic)

# Создать директорию HW9/runtests если её нет
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/runtests)

# Установить директорию HW9/runtests для исполняемого файла тестов
set_target_properties(hwtests PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/runtests
)

#[[
    Линковка к исполняемому файлу тестов библиотек
]]
target_link_libraries(hwtests PRIVATE catch_lib_test 
                                      c_lib_test
                                      tools_lib_test
                                      )

# Правильное имя исполняемого файла для разных ОС
if(WIN32)
    set(TEST_EXECUTABLE ${CMAKE_SOURCE_DIR}/runtests/hwtests.exe)
else()
    set(TEST_EXECUTABLE ${CMAKE_SOURCE_DIR}/runtests/hwtests)
endif()

# Добавить тест
add_test(NAME AllHWTests COMMAND ${TEST_EXECUTABLE})