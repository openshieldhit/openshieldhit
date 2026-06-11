# run_gen_stress.cmake — CTest driver for the C8 geometry generator.
#
# Generates a ~200-zone sphere-lattice case with gen_gemca_stress.py and
# validates it with `openshieldhit --dry-run`, so a geometry-parser change
# that breaks generated cases fails in CI rather than at profiling time.

cmake_minimum_required(VERSION 3.14)

foreach(required_var BENCH_PYTHON OSH_EXECUTABLE GEN_SCRIPT WORK_DIR)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_gen_stress.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

execute_process(
    COMMAND "${BENCH_PYTHON}" "${GEN_SCRIPT}" --zones 200 --out "${WORK_DIR}"
    RESULT_VARIABLE gen_rc
    OUTPUT_VARIABLE gen_out
    ERROR_VARIABLE gen_err
)
if(NOT gen_rc EQUAL 0)
    message(FATAL_ERROR "gen_gemca_stress.py failed (rc=${gen_rc}):\n${gen_out}\n${gen_err}")
endif()

execute_process(
    COMMAND "${OSH_EXECUTABLE}" --dry-run -v "${WORK_DIR}"
    RESULT_VARIABLE parse_rc
    OUTPUT_VARIABLE parse_out
    ERROR_VARIABLE parse_err
)
if(NOT parse_rc EQUAL 0)
    message(FATAL_ERROR "generated geometry rejected by the parser (rc=${parse_rc}):\n${parse_out}\n${parse_err}")
endif()
