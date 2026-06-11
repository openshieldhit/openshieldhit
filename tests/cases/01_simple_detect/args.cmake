# --profile together with --dry-run exercises the "profile requested but
# transport skipped" warning path; no profile file is written.
set(OSH_ARGS
    "--dry-run"
    "-v"
    "--outdir" "${WORK_DIR}"
    "--profile" "${WORK_DIR}/profile.json"
)
