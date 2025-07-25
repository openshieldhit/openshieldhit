file(GLOB TEST_SOURCES test_*.c)

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_include_directories(${test_name}
        PRIVATE ${PROJECT_SOURCE_DIR}/include
                ${PROJECT_SOURCE_DIR}/src/common
    )
    target_link_libraries(${test_name}
        PRIVATE osh_common
    )
    add_test(NAME ${test_name}_tests COMMAND ${test_name})
endforeach()