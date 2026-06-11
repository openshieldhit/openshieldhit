# run_idd.cmake — CTest driver for IDD reference benchmarks
#
# Invoked by CMake as:
#   cmake -DOSH_EXECUTABLE=<path> -DCASE_DIR=<path> -DWORK_DIR=<path> -P run_idd.cmake
#
# What it does:
#   1. Runs the openshieldhit executable on the case inputs (full transport).
#   2. Compares the produced idd.dat against every reference curve in
#      <case>/reference/*.dat using compare_idd.py (absolute per-primary
#      comparison: integral, per-bin, and peak-position criteria).
#   3. If no reference curve is present the test prints "SKIPPED" and the
#      CTest SKIP_REGULAR_EXPRESSION property marks it as skipped — the
#      benchmark infrastructure can be merged before the reference data
#      (SH12A, later FLUKA) is available.
#
# Adding a reference:
#   - Run the equivalent setup in the reference code (see the case README).
#   - Export the IDD as text: either two columns (z [cm], energy per primary
#     per bin) or openshieldhit mesh format (X Y Z E ...).
#   - Save it as <case>/reference/idd_<code>.dat (e.g. idd_sh12a.dat).
#   - Optional per-case tolerances: set COMPARE_ARGS in <case>/compare_args.cmake.

cmake_minimum_required(VERSION 3.16)

foreach(required_var OSH_EXECUTABLE CASE_DIR WORK_DIR)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_idd.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

find_program(PYTHON_EXECUTABLE NAMES python3 python)
if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "run_idd.cmake: python3 not found — required for IDD comparison")
endif()

file(MAKE_DIRECTORY "${WORK_DIR}")

# ---- Reference curves (checked first so a missing-ref skip stays cheap) ----
file(GLOB _ref_curves "${CASE_DIR}/reference/*.dat")
if(NOT _ref_curves)
    message("SKIPPED: no reference data in ${CASE_DIR}/reference/")
    return()
endif()

# ---- Run the simulation -----------------------------------------------------
execute_process(
    COMMAND "${OSH_EXECUTABLE}" -v --outdir "${WORK_DIR}" "${CASE_DIR}"
    OUTPUT_FILE "${WORK_DIR}/stdout.txt"
    ERROR_FILE  "${WORK_DIR}/stderr.txt"
    RESULT_VARIABLE exit_code
)
if(NOT exit_code EQUAL 0)
    file(READ "${WORK_DIR}/stderr.txt" _stderr)
    message(FATAL_ERROR "openshieldhit failed (exit ${exit_code}) for '${CASE_DIR}':\n${_stderr}")
endif()

if(NOT EXISTS "${WORK_DIR}/idd.dat")
    message(FATAL_ERROR "expected output idd.dat was not produced in ${WORK_DIR}")
endif()

# ---- Per-case tolerance overrides -------------------------------------------
set(COMPARE_ARGS "")
if(EXISTS "${CASE_DIR}/compare_args.cmake")
    include("${CASE_DIR}/compare_args.cmake")
endif()

# ---- Compare against every reference curve ----------------------------------
foreach(_ref IN LISTS _ref_curves)
    get_filename_component(_ref_name "${_ref}" NAME_WE)
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}"
            "${CMAKE_CURRENT_LIST_DIR}/compare_idd.py"
            "${WORK_DIR}/idd.dat"
            "${_ref}"
            --label "${_ref_name}"
            ${COMPARE_ARGS}
        RESULT_VARIABLE _cmp_rc
    )
    if(NOT _cmp_rc EQUAL 0)
        message(FATAL_ERROR "IDD comparison against '${_ref_name}' failed for '${CASE_DIR}'")
    endif()
endforeach()
