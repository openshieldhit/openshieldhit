# run_variance.cmake — CTest driver for the Monte-Carlo standard-error feature
# (issue #209).
#
# Invoked by CMake as:
#   cmake -DOSH_EXECUTABLE=<path> -DCASE_DIR=<path> -DWORK_DIR=<path>
#         -DNSTAT=<n> [-DREPLICAS=<n>] -P run_variance.cmake
#
# It copies the case, enables error tracking per estimator by injecting a
# "Settings / Variance On" block into the copied detect.dat and tagging every
# Quantity line with it, then runs real transport.  It then asserts:
#
#   * a plain (no-variance) run and the variance run agree on every value column
#     (batching only reorders the FP summation) and the variance run adds a sane,
#     non-negative, not-all-zero standard-error column per quantity — checked by
#     check_variance.py;
#   * the variance run is BYTE-IDENTICAL across two invocations — the derived
#     count-cadence partition is deterministic, the reproducibility guard of #168
#     (and the regression guard against leftover shared-mutable state).
#
# When -DREPLICAS=<n> is given it also runs the same variance case under
# --score-replicas <n> and checks its error columns the same way, confirming the
# replica split doubles as the batch-means source (no dump cadence required).

cmake_minimum_required(VERSION 3.14)

foreach(required_var OSH_EXECUTABLE CASE_DIR WORK_DIR NSTAT)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_variance.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

find_program(PYTHON_EXECUTABLE NAMES python3 python)
if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "run_variance.cmake: python3 not found — required for .dat comparison")
endif()

# ---- Case sources: a plain copy, and a copy with per-estimator variance on ----
set(PLAIN_SRC "${WORK_DIR}/case_plain")
set(VAR_SRC "${WORK_DIR}/case_variance")
file(REMOVE_RECURSE "${PLAIN_SRC}" "${VAR_SRC}")
file(COPY "${CASE_DIR}/" DESTINATION "${PLAIN_SRC}")
file(COPY "${CASE_DIR}/" DESTINATION "${VAR_SRC}")
if(NOT EXISTS "${VAR_SRC}/detect.dat")
    message(FATAL_ERROR "run_variance.cmake: case '${CASE_DIR}' has no detect.dat")
endif()

# Enable variance per estimator: prepend a Settings block carrying "Variance On"
# and tag every Quantity line with its name ("Quantity <q>" -> "Quantity <q> withErr"),
# so each scoring page opts in through the same mechanism as "Quantity Dose inWater".
file(READ "${VAR_SRC}/detect.dat" _deck)
string(REGEX REPLACE "([Qq]uantity[ \t]+[A-Za-z0-9_]+)" "\\1 withErr" _deck "${_deck}")
set(_deck "Settings\n    Name withErr\n    Variance On\n\n${_deck}")
file(WRITE "${VAR_SRC}/detect.dat" "${_deck}")

# ---- Helper: run a case source with extra args into a fresh sub-dir -----------
function(run_variant src subdir)
    set(out_dir "${WORK_DIR}/${subdir}")
    file(REMOVE_RECURSE "${out_dir}")
    file(MAKE_DIRECTORY "${out_dir}")
    execute_process(
        COMMAND "${OSH_EXECUTABLE}" -n "${NSTAT}" ${ARGN} --outdir "${out_dir}" "${src}"
        OUTPUT_FILE "${out_dir}/stdout.txt"
        ERROR_FILE  "${out_dir}/stderr.txt"
        RESULT_VARIABLE _rc
    )
    if(NOT _rc EQUAL 0)
        file(READ "${out_dir}/stderr.txt" _stderr)
        message(FATAL_ERROR "run '${subdir}' failed (exit ${_rc}) for '${src}':\n${_stderr}")
    endif()
endfunction()

run_variant("${PLAIN_SRC}" "plain")
run_variant("${VAR_SRC}"   "var_a")
run_variant("${VAR_SRC}"   "var_b")

# Optional --score-replicas variant: with per-estimator variance tracking active,
# the replica split doubles as the batch-means source (issue #209) — no dump cadence needed.
# It must produce the same sane error columns as the derived-cadence run.
if(DEFINED REPLICAS)
    run_variant("${VAR_SRC}" "var_repl" --score-replicas "${REPLICAS}")
endif()

# ---- Validate the error columns against the plain run ------------------------
file(GLOB _plain_dats "${WORK_DIR}/plain/*.dat")
if(_plain_dats STREQUAL "")
    message(FATAL_ERROR "no *.dat output produced by the plain run of '${CASE_DIR}'")
endif()
foreach(_plain IN LISTS _plain_dats)
    get_filename_component(_name "${_plain}" NAME)
    set(_var "${WORK_DIR}/var_a/${_name}")
    if(NOT EXISTS "${_var}")
        message(FATAL_ERROR "VARIANCE run did not produce ${_name}")
    endif()
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/check_variance.py"
                "${_var}" "${_plain}" --label "${_name}"
        RESULT_VARIABLE _chk_rc
    )
    if(NOT _chk_rc EQUAL 0)
        message(FATAL_ERROR "standard-error validation failed for ${_name}")
    endif()
endforeach()

# ---- Validate the --score-replicas run's error columns too (if requested) ----
if(DEFINED REPLICAS)
    foreach(_plain IN LISTS _plain_dats)
        get_filename_component(_name "${_plain}" NAME)
        set(_var "${WORK_DIR}/var_repl/${_name}")
        if(NOT EXISTS "${_var}")
            message(FATAL_ERROR "--score-replicas VARIANCE run did not produce ${_name}")
        endif()
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/check_variance.py"
                    "${_var}" "${_plain}" --label "replicas:${_name}"
            RESULT_VARIABLE _chk_rc
        )
        if(NOT _chk_rc EQUAL 0)
            message(FATAL_ERROR "standard-error validation failed for ${_name} (--score-replicas)")
        endif()
    endforeach()
endif()

# ---- Reproducibility: the two VARIANCE runs must be byte-identical -----------
file(GLOB _var_dats "${WORK_DIR}/var_a/*.dat")
foreach(_ref IN LISTS _var_dats)
    get_filename_component(_name "${_ref}" NAME)
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/compare_dat.py"
                "${WORK_DIR}/var_b/${_name}" "${_ref}" --exact --label "var-repro:${_name}"
        RESULT_VARIABLE _cmp_rc
    )
    if(NOT _cmp_rc EQUAL 0)
        message(FATAL_ERROR "VARIANCE run is not reproducible across runs for ${_name}")
    endif()
endforeach()

message(STATUS "variance error-bar + reproducibility OK for '${CASE_DIR}' (nstat=${NSTAT})")
