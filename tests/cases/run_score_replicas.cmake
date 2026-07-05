# run_score_replicas.cmake — CTest driver for the --score-replicas reproducibility
# contract (issue #230).
#
# Invoked by CMake as:
#   cmake -DOSH_EXECUTABLE=<path> -DCASE_DIR=<path> -DWORK_DIR=<path>
#         -DNSTAT=<n> [-DREPLICAS=<n>] -P run_score_replicas.cmake
#
# It runs the same case four ways (real transport, same nstat + seed) and asserts
# the deterministic core of the parallel contract (issue #168), without threads:
#
#   serial (no flag)  vs  --score-replicas 1   → BYTE-IDENTICAL data
#       one range → private set → merge into an empty master is the same
#       summation order as serial.
#   serial (no flag)  vs  --score-replicas N   → within compare_dat.py tolerance
#       same per-history physics; only cross-partition FP summation order differs.
#   --score-replicas N (run twice)             → BYTE-IDENTICAL data
#       fixed partition ⇒ deterministic merge order (nondeterminism guard).
#
# Every TEXT *.dat the serial run emits is compared, so multi-output cases (mesh +
# cyl) are all covered.

cmake_minimum_required(VERSION 3.14)

foreach(required_var OSH_EXECUTABLE CASE_DIR WORK_DIR NSTAT)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_score_replicas.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

if(NOT DEFINED REPLICAS)
    set(REPLICAS 4)
endif()

find_program(PYTHON_EXECUTABLE NAMES python3 python)
if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "run_score_replicas.cmake: python3 not found — required for .dat comparison")
endif()

# ---- Helper: run the case with the given extra args into a fresh sub-dir ------
function(run_variant subdir)
    set(out_dir "${WORK_DIR}/${subdir}")
    file(REMOVE_RECURSE "${out_dir}")
    file(MAKE_DIRECTORY "${out_dir}")
    execute_process(
        COMMAND "${OSH_EXECUTABLE}" -n "${NSTAT}" ${ARGN} --outdir "${out_dir}" "${CASE_DIR}"
        OUTPUT_FILE "${out_dir}/stdout.txt"
        ERROR_FILE  "${out_dir}/stderr.txt"
        RESULT_VARIABLE _rc
    )
    if(NOT _rc EQUAL 0)
        file(READ "${out_dir}/stderr.txt" _stderr)
        message(FATAL_ERROR "run '${subdir}' failed (exit ${_rc}) for '${CASE_DIR}':\n${_stderr}")
    endif()
endfunction()

# ---- Helper: compare every serial *.dat against another run's copy -----------
function(compare_variant subdir mode)
    file(GLOB _dats "${WORK_DIR}/serial/*.dat")
    if(_dats STREQUAL "")
        message(FATAL_ERROR "no *.dat output produced by the serial run of '${CASE_DIR}'")
    endif()
    foreach(_ref IN LISTS _dats)
        get_filename_component(_name "${_ref}" NAME)
        set(_other "${WORK_DIR}/${subdir}/${_name}")
        if(NOT EXISTS "${_other}")
            message(FATAL_ERROR "'${subdir}' did not produce ${_name}")
        endif()
        if(mode STREQUAL "exact")
            set(_cmp_args --exact --label "${subdir}:${_name}")
        else()
            set(_cmp_args --rtol 0.02 --label "${subdir}:${_name}")
        endif()
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/compare_dat.py"
                    "${_other}" "${_ref}" ${_cmp_args}
            RESULT_VARIABLE _cmp_rc
        )
        if(NOT _cmp_rc EQUAL 0)
            message(FATAL_ERROR "comparison failed (${mode}) for ${_name}: serial vs ${subdir}")
        endif()
    endforeach()
endfunction()

# ---- Run the four variants ---------------------------------------------------
run_variant("serial")
run_variant("r1"   --score-replicas 1)
run_variant("rN"   --score-replicas "${REPLICAS}")
run_variant("rN_b" --score-replicas "${REPLICAS}")

# ---- Assert the contract -----------------------------------------------------
compare_variant("r1"   exact)      # N==1 reproduces serial bit-for-bit
compare_variant("rN"   tolerance)  # N>1 matches serial within FP tolerance
# The two N-runs must equal each other exactly; compare rN_b against rN's files.
file(GLOB _rn_dats "${WORK_DIR}/rN/*.dat")
foreach(_ref IN LISTS _rn_dats)
    get_filename_component(_name "${_ref}" NAME)
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/compare_dat.py"
                "${WORK_DIR}/rN_b/${_name}" "${_ref}" --exact --label "rN_b-vs-rN:${_name}"
        RESULT_VARIABLE _cmp_rc
    )
    if(NOT _cmp_rc EQUAL 0)
        message(FATAL_ERROR "score-replicas ${REPLICAS} is not reproducible across runs for ${_name}")
    endif()
endforeach()

message(STATUS "score-replicas reproducibility OK for '${CASE_DIR}' (nstat=${NSTAT}, N=${REPLICAS})")
