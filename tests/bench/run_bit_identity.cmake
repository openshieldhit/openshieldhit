# run_bit_identity.cmake — CTest driver for the profiling no-effect guarantee.
#
# Issue #138 (Phase 0) requires that profiled runs are bit-identical to
# unprofiled ones: the phase timers must never touch the RNG streams or any
# physics state.  This driver runs the same case twice — once plain, once
# with --profile — and compares the scored text output exactly, ignoring
# only '#' comment lines (they embed the wall-clock timestamp, which differs
# between any two runs).

cmake_minimum_required(VERSION 3.14)

foreach(required_var BENCH_PYTHON OSH_EXECUTABLE CASE_DIR WORK_DIR)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_bit_identity.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

file(MAKE_DIRECTORY "${WORK_DIR}/plain")
file(MAKE_DIRECTORY "${WORK_DIR}/profiled")

execute_process(
    COMMAND "${OSH_EXECUTABLE}" --workdir "${CASE_DIR}" -n 500 -o "${WORK_DIR}/plain"
    RESULT_VARIABLE rc_plain
    OUTPUT_QUIET
    ERROR_VARIABLE err_plain
)
if(NOT rc_plain EQUAL 0)
    message(FATAL_ERROR "plain run failed (rc=${rc_plain}):\n${err_plain}")
endif()

execute_process(
    COMMAND "${OSH_EXECUTABLE}" --workdir "${CASE_DIR}" -n 500
        -o "${WORK_DIR}/profiled" --profile "${WORK_DIR}/profiled/profile.json"
    RESULT_VARIABLE rc_prof
    OUTPUT_QUIET
    ERROR_VARIABLE err_prof
)
if(NOT rc_prof EQUAL 0)
    message(FATAL_ERROR "profiled run failed (rc=${rc_prof}):\n${err_prof}")
endif()

execute_process(
    COMMAND "${BENCH_PYTHON}" -c
"import sys, pathlib
plain, profiled = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
names = sorted(p.name for p in plain.glob('*.dat'))
if not names:
    sys.exit('no .dat outputs found in ' + str(plain))
for name in names:
    strip = lambda p: [l for l in p.read_text().splitlines() if not l.lstrip().startswith('#')]
    if strip(plain / name) != strip(profiled / name):
        sys.exit(name + ': profiled output differs from unprofiled output')
print('bit-identical:', ', '.join(names))
"
        "${WORK_DIR}/plain" "${WORK_DIR}/profiled"
    RESULT_VARIABLE rc_cmp
    OUTPUT_VARIABLE out_cmp
    ERROR_VARIABLE err_cmp
)
if(NOT rc_cmp EQUAL 0)
    message(FATAL_ERROR "bit-identity violated:\n${out_cmp}\n${err_cmp}")
endif()
message(STATUS "${out_cmp}")
