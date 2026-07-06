# run_case.cmake — CTest driver for integration cases
#
# Invoked by CMake as:
#   cmake -DOSH_EXECUTABLE=<path> -DCASE_DIR=<path> -DWORK_DIR=<path> -P run_case.cmake
#
# What it does:
#   1. Runs the openshieldhit executable with the case directory as working dir.
#   2. Checks the exit code against expected/exit_code.txt (default: 0).
#   3. Compares stdout/stderr against expected/ references if they exist.
#   4. Compares any expected/*.dat files against the corresponding output files
#      in WORK_DIR, ignoring comment lines (#).  The reference and the test run
#      must use the same RNG seed and nstat for data lines to match exactly.
#
# Adding output comparison for a case:
#   - Capture a known-good run with the desired nstat and seed, then place the
#     output .dat files in expected/.
#   - Commit the expected/ files alongside the input files.
#
# Binary output files (*.bdo) are not compared here.

cmake_minimum_required(VERSION 3.14)

# ---- Validate parameters --------------------------------------------------
foreach(required_var OSH_EXECUTABLE CASE_DIR WORK_DIR)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_case.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

find_program(PYTHON_EXECUTABLE NAMES python3 python)
if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "run_case.cmake: python3 not found — required for .dat comparison")
endif()

# ---- Prepare work directory -----------------------------------------------
file(MAKE_DIRECTORY "${WORK_DIR}")

# ---- Determine CLI arguments ----------------------------------------------
# Default: --dry-run (OPENSHIELDHIT_RUN_VALIDATE); no transport needed.
# --outdir routes any output files to the isolated work directory so nothing
# is written back into the source-controlled case directory.
# A case can override OSH_ARGS entirely via a local args.cmake.
set(OSH_ARGS "--dry-run" "-v" "--outdir" "${WORK_DIR}")
if(EXISTS "${CASE_DIR}/args.cmake")
    include("${CASE_DIR}/args.cmake")
endif()

# ---- Run ------------------------------------------------------------------
execute_process(
    COMMAND "${OSH_EXECUTABLE}" ${OSH_ARGS} "${CASE_DIR}"
    OUTPUT_FILE "${WORK_DIR}/stdout.txt"
    ERROR_FILE  "${WORK_DIR}/stderr.txt"
    RESULT_VARIABLE exit_code
)

# ---- Check exit code ------------------------------------------------------
set(expected_exit 0)
if(EXISTS "${CASE_DIR}/expected/exit_code.txt")
    file(READ "${CASE_DIR}/expected/exit_code.txt" _raw)
    string(STRIP "${_raw}" expected_exit)
endif()

if(NOT exit_code EQUAL expected_exit)
    file(READ "${WORK_DIR}/stderr.txt" _stderr)
    message(FATAL_ERROR
        "Exit code mismatch for '${CASE_DIR}':\n"
        "  expected : ${expected_exit}\n"
        "  got      : ${exit_code}\n"
        "stderr:\n${_stderr}"
    )
endif()

# ---- Compare stdout (if reference exists) ---------------------------------
if(EXISTS "${CASE_DIR}/expected/stdout.txt")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E compare_files
            "${WORK_DIR}/stdout.txt"
            "${CASE_DIR}/expected/stdout.txt"
        RESULT_VARIABLE _diff
    )
    if(NOT _diff EQUAL 0)
        file(READ "${WORK_DIR}/stdout.txt"          _actual)
        file(READ "${CASE_DIR}/expected/stdout.txt" _expected)
        message(FATAL_ERROR
            "stdout mismatch for '${CASE_DIR}':\n"
            "--- expected ---\n${_expected}"
            "--- actual ---\n${_actual}"
        )
    endif()
endif()

# ---- Compare stderr (if reference exists) ---------------------------------
if(EXISTS "${CASE_DIR}/expected/stderr.txt")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E compare_files
            "${WORK_DIR}/stderr.txt"
            "${CASE_DIR}/expected/stderr.txt"
        RESULT_VARIABLE _diff
    )
    if(NOT _diff EQUAL 0)
        file(READ "${WORK_DIR}/stderr.txt"          _actual)
        file(READ "${CASE_DIR}/expected/stderr.txt" _expected)
        message(FATAL_ERROR
            "stderr mismatch for '${CASE_DIR}':\n"
            "--- expected ---\n${_expected}"
            "--- actual ---\n${_actual}"
        )
    endif()
endif()

# ---- Compare ASCII output .dat files (if references exist) ----------------
# compare_dat.py checks each numeric column within a 2 % relative tolerance,
# skipping comment lines (#).  The reference must be generated with the same
# geometry and particle type; nstat and RNG seed may differ.
#
# The number of leading coordinate columns defaults to 3 (Mesh X/Y/Z).  A case
# whose output has a different coordinate layout — e.g. Zone output, whose single
# leading column is the zone id — can override it with expected/coord_cols.txt.
set(_coord_cols 3)
if(EXISTS "${CASE_DIR}/expected/coord_cols.txt")
    file(READ "${CASE_DIR}/expected/coord_cols.txt" _raw_cc)
    string(STRIP "${_raw_cc}" _coord_cols)
endif()
file(GLOB _expected_dats "${CASE_DIR}/expected/*.dat")
foreach(_ref_dat IN LISTS _expected_dats)
    get_filename_component(_dat_name "${_ref_dat}" NAME)
    set(_actual_dat "${WORK_DIR}/${_dat_name}")
    if(NOT EXISTS "${_actual_dat}")
        message(FATAL_ERROR
            "Expected output file was not produced: ${_dat_name}\n"
            "  looked in: ${WORK_DIR}\n"
            "  reference: ${_ref_dat}"
        )
    endif()
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}"
            "${CMAKE_CURRENT_LIST_DIR}/compare_dat.py"
            "${_actual_dat}"
            "${_ref_dat}"
            --rtol 0.02
            --coord-cols ${_coord_cols}
            --label "${_dat_name}"
        RESULT_VARIABLE _cmp_rc
    )
    if(NOT _cmp_rc EQUAL 0)
        message(FATAL_ERROR "Data file comparison failed for '${_dat_name}' in '${CASE_DIR}'")
    endif()
endforeach()
