#[[ Сборка ДЗ с GUI для HW9 ]]

add_library(commonc STATIC ${C_SRC_PATH}/common.c)
target_include_directories(commonc PUBLIC ${C_SRC_PATH}/includes)

foreach(CPP_GUI_FILE ${SOURCES_GUI})
    get_filename_component(TARGET_NAME ${CPP_GUI_FILE} NAME_WE)

    add_executable(${TARGET_NAME} ${CPP_GUI_FILE})
    target_include_directories(${TARGET_NAME} PRIVATE 
                                ${fltk_INCLUDE_DIR} 
                                ${CMAKE_CURRENT_SOURCE_DIR}/includes
                                ${C_SRC_PATH}/includes
    )

    # Иконка только для Windows
    if(WIN32)
        target_sources(${TARGET_NAME} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/includes/res/icon.rc
        )
    endif()

    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/gui
    )

    if(WIN32)
        set_target_properties(${TARGET_NAME} PROPERTIES
            WIN32_EXECUTABLE TRUE
        )
    endif()

    target_compile_options(${TARGET_NAME} PRIVATE -Wall -Wextra -Wpedantic)

    target_link_libraries(${TARGET_NAME} PRIVATE
        fltk::fltk
        commonc
    )

    # libm нужна только на Unix-подобных
    if(UNIX)
        target_link_libraries(${TARGET_NAME} PRIVATE m)
    endif()
endforeach()