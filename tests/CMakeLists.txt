function(register_named_c_tests target source_file)
    file(STRINGS "${source_file}" named_test_lines REGEX "^static void test_[A-Za-z0-9_]+\\(void\\)")

    foreach(test_line IN LISTS named_test_lines)
        string(REGEX MATCH "test_([A-Za-z0-9_]+)\\(void\\)" _ "${test_line}")
        set(test_case "${CMAKE_MATCH_1}")
        add_test(NAME "${target}::${test_case}" COMMAND "${target}" "${test_case}")
    endforeach()
endfunction()

# ---- test_osh_main: verifies shared version header values ----
add_executable(test_osh_main test_osh_main.c)
target_link_libraries(test_osh_main PRIVATE openshieldhit_api)
register_named_c_tests(test_osh_main ${CMAKE_CURRENT_SOURCE_DIR}/test_osh_main.c)
add_test(NAME openshieldhit_version_cli COMMAND $<TARGET_FILE:openshieldhit> --version)

# Collect remaining test source files
file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.c")
list(REMOVE_ITEM TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/test_osh_main.c)

# Loop over each test file and create an executable
foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)

    add_executable(${test_name} ${test_src})

    # Include all relevant source subdirectories
    target_include_directories(${test_name}
        PRIVATE
            ${PROJECT_SOURCE_DIR}/src
    )

    # Link all needed internal libraries
    target_link_libraries(${test_name}
        PRIVATE
            osh_common
            osh_random
            osh_beam
            osh_particle
    )

    # Register the test
    add_test(NAME ${test_name}_tests COMMAND ${test_name})
endforeach()
