# ---- Shared helper: register individual test functions as named CTest cases ----
# Scans the source file for "static void test_<name>(void)" signatures and
# registers each as its own CTest entry, so failures are reported per-function.
function(register_named_c_tests target source_file)
    file(STRINGS "${source_file}" named_test_lines
        REGEX "^static void test_[A-Za-z0-9_]+\\(void\\)")
    foreach(test_line IN LISTS named_test_lines)
        string(REGEX MATCH "test_([A-Za-z0-9_]+)\\(void\\)" _ "${test_line}")
        set(test_case "${CMAKE_MATCH_1}")
        add_test(NAME "unit::${target}::${test_case}"
                 COMMAND $<TARGET_FILE:${target}> "${test_case}")
        set_tests_properties("unit::${target}::${test_case}" PROPERTIES LABELS "unit")
    endforeach()
endfunction()

# ---- test_osh_main: public API (version, context lifecycle, setters) ----
# Registered three times — once per verbosity level (0=silent, 1=info, 2=debug)
# so that logging paths are exercised under each level.
add_executable(test_osh_main test_osh_main.c)
target_link_libraries(test_osh_main PRIVATE openshieldhit_api)
target_compile_definitions(test_osh_main PRIVATE
    OSH_PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)
foreach(verbosity 0 1 2)
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/test_osh_main.c" named_test_lines
        REGEX "^static void test_[A-Za-z0-9_]+\\(void\\)")
    foreach(test_line IN LISTS named_test_lines)
        string(REGEX MATCH "test_([A-Za-z0-9_]+)\\(void\\)" _ "${test_line}")
        set(test_case "${CMAKE_MATCH_1}")
        add_test(NAME "unit::test_osh_main::v${verbosity}::${test_case}"
                 COMMAND $<TARGET_FILE:test_osh_main> "${verbosity}" "${test_case}")
        set_tests_properties("unit::test_osh_main::v${verbosity}::${test_case}"
                             PROPERTIES LABELS "unit")
    endforeach()
endforeach()

# ---- test_osh_cli: internal CLI parser ----
add_executable(test_osh_cli test_osh_cli.c)
target_link_libraries(test_osh_cli PRIVATE osh_cli)
register_named_c_tests(test_osh_cli ${CMAKE_CURRENT_SOURCE_DIR}/test_osh_cli.c)

# ---- Remaining unit tests: auto-discovered by glob ----
file(GLOB UNIT_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.c")
list(REMOVE_ITEM UNIT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/test_osh_main.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_osh_cli.c
)

foreach(test_src ${UNIT_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)

    add_executable(${test_name} ${test_src})

    target_include_directories(${test_name} PRIVATE
        ${PROJECT_SOURCE_DIR}/src
    )

    target_link_libraries(${test_name} PRIVATE
        osh_common
        osh_random
        osh_beam
        osh_material
        osh_particle
    )

    # Pass the fixtures directory so tests can locate input files without
    # relying on fragile relative paths from the build directory.
    target_compile_definitions(${test_name} PRIVATE
        OSH_TEST_FIXTURES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../fixtures"
        OSH_PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
    )

    add_test(NAME "unit::${test_name}" COMMAND $<TARGET_FILE:${test_name}>)
    set_tests_properties("unit::${test_name}" PROPERTIES LABELS "unit")
endforeach()

# ---- CLI smoke test via the built executable ----
add_test(NAME "unit::cli::version" COMMAND $<TARGET_FILE:openshieldhit> --version)
add_test(NAME "unit::cli::help"    COMMAND $<TARGET_FILE:openshieldhit> --help)
set_tests_properties("unit::cli::version" "unit::cli::help" PROPERTIES LABELS "unit")
