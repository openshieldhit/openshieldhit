# run_capacity_invariance.cmake — CTest driver for the pool-capacity
# invariance guarantee (issue #148, RNG parallelization Phase 0).
#
# Per-history RNG seeding makes the random sequence each history consumes a
# function of its global index alone, never of execution order.  Consequently
# the physics outcome of a run is independent of the transport pool capacity
# (the number of histories alive simultaneously) — capacity is purely a
# cache/parallelism performance knob.
#
# This driver runs the same case at several pool capacities and compares the
# scored .dat outputs.  The case is deliberately non-nuclear (NUCRE 0) so that
# no secondaries are produced: the only remaining capacity-dependent effect is
# the floating-point summation order in the shared scoring accumulators, which
# is bounded far below any real physics difference.  We therefore compare
# numeric fields within a tight relative tolerance rather than byte-for-byte.
#
# Exact byte identity across capacities additionally requires a deterministic
# scoring reduction (per-thread tallies + fixed-order merge), tracked as the
# out-of-scope follow-up to issue #148.

cmake_minimum_required(VERSION 3.14)

foreach(required_var BENCH_PYTHON OSH_EXECUTABLE CASE_DIR WORK_DIR)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "run_capacity_invariance.cmake: required parameter -D${required_var}=... not set")
    endif()
endforeach()

# Baseline first, then the capacities compared against it.  Spanning 1 (fully
# serial: one history at a time) to a value >= nstat (all histories live at
# once) exercises the extremes of the wavefront/compaction machinery.
set(NSTAT 400)
set(BASE_CAP 1)
set(COMPARE_CAPS 8 64 4096)
set(REL_TOL 1e-9)

function(run_case cap out_dir)
    file(MAKE_DIRECTORY "${out_dir}")
    execute_process(
        COMMAND "${OSH_EXECUTABLE}" --workdir "${CASE_DIR}" -n ${NSTAT}
            -o "${out_dir}" --pool-capacity ${cap}
        RESULT_VARIABLE rc
        ERROR_VARIABLE err
        OUTPUT_QUIET
    )
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "run at pool-capacity=${cap} failed (rc=${rc}):\n${err}")
    endif()
endfunction()

run_case(${BASE_CAP} "${WORK_DIR}/cap_${BASE_CAP}")
foreach(cap ${COMPARE_CAPS})
    run_case(${cap} "${WORK_DIR}/cap_${cap}")
endforeach()

# Compare every compared capacity's .dat outputs against the baseline.
foreach(cap ${COMPARE_CAPS})
    execute_process(
        COMMAND "${BENCH_PYTHON}" -c
"import sys, pathlib
base, other, tol = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), float(sys.argv[3])
names = sorted(p.name for p in base.glob('*.dat'))
if not names:
    sys.exit('no .dat outputs found in ' + str(base))
other_names = sorted(p.name for p in other.glob('*.dat'))
if other_names != names:
    sys.exit('output set differs: baseline %r vs %r' % (names, other_names))
def rows(p):
    return [l.split() for l in p.read_text().splitlines()
            if l.strip() and not l.lstrip().startswith('#')]
worst = 0.0
for name in names:
    a, b = rows(base / name), rows(other / name)
    if len(a) != len(b):
        sys.exit(name + ': row count differs (%d vs %d)' % (len(a), len(b)))
    for i, (ra, rb) in enumerate(zip(a, b)):
        if len(ra) != len(rb):
            sys.exit('%s row %d: column count differs' % (name, i))
        for xa, xb in zip(ra, rb):
            try:
                fa, fb = float(xa), float(xb)
            except ValueError:
                if xa != xb:
                    sys.exit('%s row %d: non-numeric field differs (%r vs %r)' % (name, i, xa, xb))
                continue
            denom = abs(fa) + abs(fb)
            rel = 0.0 if denom == 0.0 else abs(fa - fb) / denom
            if rel > worst:
                worst = rel
            if rel > tol:
                sys.exit('%s row %d: %r vs %r relative diff %.3e exceeds tol %.1e'
                         % (name, i, xa, xb, rel, tol))
print('capacity-invariant within %.1e (worst observed %.3e): %s'
      % (tol, worst, ', '.join(names)))
"
            "${WORK_DIR}/cap_${BASE_CAP}" "${WORK_DIR}/cap_${cap}" "${REL_TOL}"
        RESULT_VARIABLE rc_cmp
        OUTPUT_VARIABLE out_cmp
        ERROR_VARIABLE err_cmp
    )
    if(NOT rc_cmp EQUAL 0)
        message(FATAL_ERROR "capacity invariance violated (cap=${BASE_CAP} vs cap=${cap}):\n${out_cmp}\n${err_cmp}")
    endif()
    message(STATUS "cap=${BASE_CAP} vs cap=${cap}: ${out_cmp}")
endforeach()
