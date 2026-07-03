# Wall-time budget integration case (issue #192 / #195).
#
# Requests far more primaries than can finish in the budget and uses a small
# wavefront pool so the run is forced to stop cleanly mid-injection.  The run
# must still exit 0: in-flight histories drain, all secondary families drain,
# and the partial result is saved normalised by the true completed-primary
# count.  The exact count is timing-dependent (a wall-time budget is not
# deterministic by design — see issue #195), so this case only asserts the
# end-to-end exit status; the precise, deterministic count is verified in the
# unit test test_osh_simulation::clean_stop_partial_result_is_exact.
set(OSH_ARGS
    "-v"
    "-n" "5000000"
    "--pool-capacity" "200"
    "--max-time" "0.2"
    "--outdir" "${WORK_DIR}"
)
