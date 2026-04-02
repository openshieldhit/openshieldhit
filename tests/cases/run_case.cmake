# run_case.cmake — CTest driver for integration cases
#
# Invoked by CMake as:
#   cmake -DOSH_EXECUTABLE=<path> -DCASE_DIR=<path> -DWORK_DIR=<path> -P run_case.cmake
#
# What it does:
#   1. Runs the openshieldhit executable in validate mode (--dry-run) with
#      the case directory as the working directory.
#   2. Checks the exit code against expected/exit_code.txt (default: 0).
#   3. Compares stdout/stderr against expected/ references if they exist.
#
# Adding output comparison for a case:
#   - Capture a known-good run:
#       openshieldhit --dry-run <case_dir> > expected/stdout.txt 2> expected/stderr.txt
#   - Commit the expected/ files alongside the input files.
#
# Binary output files (*.bdo) are not compared here. A numeric comparison
# tool will be wired in once transport is implemented.

cmake_minimum_required(VERSION 3.14)

# ---- Validate parameters --------------------------------------------------
foreach(required_var OSH_EXECUTABLE CASE_DIR WORK_DIR)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_case.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

# ---- Prepare work directory -----------------------------------------------
file(MAKE_DIRECTORY "${WORK_DIR}")

# ---- Determine CLI arguments ----------------------------------------------
# Default: --dry-run (OPENSHIELDHIT_RUN_VALIDATE); no transport needed.
# --outdir routes any output files to the isolated work directory so nothing
# is written back into the source-controlled case directory.
# A case can override OSH_ARGS entirely via a local args.cmake.
set(OSH_ARGS "--dry-run" "--outdir" "${WORK_DIR}")
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
