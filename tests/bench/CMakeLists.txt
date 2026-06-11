# Smoke tests for the benchmark tooling in tools/bench/.
#
# These do not measure anything; they protect the harness against silent
# breakage in CI:
#   - every committed benchmark case must stay parseable (--dry-run),
#   - the C8 geometry generator must keep producing parser-valid cases,
#   - the scenario registry must stay importable and consistent.

find_program(BENCH_PYTHON NAMES python3 python)

# ---- Bench case inputs: one validation test per case directory ----
file(GLOB BENCH_CASE_DIRS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/tools/bench/cases/*/")
foreach(case_dir ${BENCH_CASE_DIRS})
    if(NOT IS_DIRECTORY "${case_dir}")
        continue()
    endif()
    get_filename_component(case_name "${case_dir}" NAME)
    add_test(
        NAME "bench::validate::${case_name}"
        COMMAND $<TARGET_FILE:openshieldhit> --dry-run -v "${case_dir}"
    )
    set_tests_properties("bench::validate::${case_name}" PROPERTIES LABELS "bench")
endforeach()

if(BENCH_PYTHON)
    # ---- C8 generator: generate a lattice case and parse it ----
    add_test(
        NAME "bench::gen_gemca_stress"
        COMMAND ${CMAKE_COMMAND}
            -DBENCH_PYTHON=${BENCH_PYTHON}
            -DOSH_EXECUTABLE=$<TARGET_FILE:openshieldhit>
            -DGEN_SCRIPT=${PROJECT_SOURCE_DIR}/tools/bench/gen_gemca_stress.py
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/gen_gemca_stress
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_gen_stress.cmake
    )
    set_tests_properties("bench::gen_gemca_stress" PROPERTIES LABELS "bench")

    # ---- Scenario registry / runner: import and self-consistency ----
    add_test(
        NAME "bench::run_bench_list"
        COMMAND ${BENCH_PYTHON} ${PROJECT_SOURCE_DIR}/tools/bench/run_bench.py --list
    )
    set_tests_properties("bench::run_bench_list" PROPERTIES LABELS "bench")

    # ---- Profiling must not change physics: bit-identical outputs ----
    add_test(
        NAME "bench::profile_bit_identity"
        COMMAND ${CMAKE_COMMAND}
            -DBENCH_PYTHON=${BENCH_PYTHON}
            -DOSH_EXECUTABLE=$<TARGET_FILE:openshieldhit>
            -DCASE_DIR=${PROJECT_SOURCE_DIR}/tests/cases/06_minimal_nucre
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/bit_identity
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_bit_identity.cmake
    )
    set_tests_properties("bench::profile_bit_identity" PROPERTIES LABELS "bench")
endif()
