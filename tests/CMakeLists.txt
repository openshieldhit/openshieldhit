# ---- test_osh_main: verifies shared version header values ----
add_executable(test_osh_main test_osh_main.c)
target_link_libraries(test_osh_main PRIVATE osh_all)
add_test(NAME test_osh_main_tests COMMAND test_osh_main)
add_test(NAME openshieldhit_version_cli COMMAND $<TARGET_FILE:openshieldhit> --version)

# Collect remaining test source files
file(GLOB TEST_SOURCES test_*.c)
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
