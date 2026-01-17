# Collect all test source files
file(GLOB TEST_SOURCES test_*.c)

# Loop over each test file and create an executable
foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)

    add_executable(${test_name} ${test_src})

    # Include all relevant source subdirectories
    target_include_directories(${test_name}
        PRIVATE
            ${PROJECT_SOURCE_DIR}/src/common
            ${PROJECT_SOURCE_DIR}/src/random
            ${PROJECT_SOURCE_DIR}/src/beam
            ${PROJECT_SOURCE_DIR}/src/particle
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
