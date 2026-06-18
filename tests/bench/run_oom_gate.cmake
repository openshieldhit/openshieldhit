# run_oom_gate.cmake — CTest driver for the out-of-memory gate (issue #152).
#
# Runs a deliberately enormous scoring configuration (huge_detect.dat, ~59.6 GiB
# of accumulators) against a small explicit --mem-budget.  The point is to prove
# that OpenShieldHIT refuses such a run *before* allocating anything and prints a
# meaningful message, rather than crashing or being OOM-killed.  Using an
# explicit budget makes the refusal deterministic on any runner, so this is safe
# to run on all CI OSes.

cmake_minimum_required(VERSION 3.14)

foreach(required_var OSH_EXECUTABLE CASE_DIR DETECT_FILE WORK_DIR)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_oom_gate.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

file(MAKE_DIRECTORY "${WORK_DIR}")

execute_process(
    COMMAND "${OSH_EXECUTABLE}" --workdir "${CASE_DIR}" -d "${DETECT_FILE}"
        -n 10 -o "${WORK_DIR}" --mem-budget 256MB
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
)

set(combined "${out}${err}")

# 1. The run must NOT succeed (it would otherwise try to allocate ~60 GiB).
if(rc EQUAL 0)
    message(FATAL_ERROR "OOM gate failed: run succeeded but should have been refused.\n${combined}")
endif()

# 2. The refusal must carry a meaningful memory message, not a random crash.
string(FIND "${combined}" "memory budget" budget_pos)
string(FIND "${combined}" "Aborting before allocating memory" abort_pos)
if(budget_pos EQUAL -1 OR abort_pos EQUAL -1)
    message(FATAL_ERROR "OOM gate refused (rc=${rc}) but without the expected message:\n${combined}")
endif()

message(STATUS "OOM gate refused the oversized run as expected (rc=${rc}).")
