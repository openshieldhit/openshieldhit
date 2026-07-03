# Periodic partial-dump integration case (issue #193).
#
# Runs 200 primaries with a deterministic count cadence (--dump-every-primaries
# 40), so the run reaches five family-complete checkpoints and overwrites the
# output file with the exact partial result at each of the first four.  Because a
# dump is a non-destructive snapshot (out-of-place postprocess of a shadow copy),
# the FINAL saved result must be byte-for-byte identical to a plain no-dump run of
# the same 200 primaries at the same seed — that reference lives in expected/
# NB_msh.dat and run_case.cmake compares against it (ignoring '#' comment lines,
# so the extra "# COMPLETENESS:" header does not affect the match).
#
# The deterministic "dump at K == truncated K-run" invariant this rests on is
# proven at the unit level (test_osh_simulation::test_dump_control_is_non_destructive
# plus the checkpoint-batching invariance tests); here we exercise the full CLI
# path — flag parsing, the file sink, and the completeness stamp — end to end.
set(OSH_ARGS
    "-v"
    "-n" "200"
    "--dump-every-primaries" "40"
    "--outdir" "${WORK_DIR}"
)
