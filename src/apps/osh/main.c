#include <stdio.h>
#include <stdlib.h>

#include "apps/osh/osh_run.h"
#include "cli/osh_cli.h"
#include "common/osh_exit.h"
#include "common/osh_logger.h"
#include "openshieldhit/file.h"
#include "openshieldhit/status.h"
#include "openshieldhit/version.h"

static int exit_code_for_status(enum osh_status status);

int main(int argc, char *argv[]) {
    struct osh_cli_options opt;
    struct osh_run_options run_opt;
    char err[256];
    enum osh_status rc;

    if (osh_cli_parse(argc, argv, &opt, err, sizeof(err)) != 0) {
        fprintf(stderr, "Error: %s\n\n", err);
        osh_cli_print_help(stderr, argv[0]);
        return EX_USAGE;
    }

    if (opt.action == OSH_CLI_ACTION_HELP) {
        osh_cli_print_help(stdout, argv[0]);
        return EX_OK;
    }

    if (opt.action == OSH_CLI_ACTION_VERSION) {
        printf("OpenShieldHIT version %s\n", osh_version_string());
        return EX_OK;
    }

    /* Initialise the logger before any library calls.
     *   no -v  → OSH_LOG_WARN  (silent; only warnings and errors)
     *      -v  → OSH_LOG_INFO  (normal informational output)
     *     -vv  → OSH_LOG_DEBUG (verbose debug output)
     * Stdout is disabled by default; enable it here for the CLI. */
    {
        int log_level = (opt.verbose == 0) ? OSH_LOG_WARN : (opt.verbose == 1) ? OSH_LOG_INFO : OSH_LOG_DEBUG;
        osh_log_init(log_level, OSH_LOG_F_NONE);
        osh_log_enable_stdout(1);
    }

    /* Normalize path separators once here so all library code can assume '/'
     * without per-subsystem #ifdef _WIN32 blocks. No-op on non-Windows. */
    osh_path_normalize((char *) opt.workdir);
    osh_path_normalize((char *) opt.geo_path);
    osh_path_normalize((char *) opt.beam_path);
    osh_path_normalize((char *) opt.mat_path);
    osh_path_normalize((char *) opt.detect_path);
    osh_path_normalize((char *) opt.out_dir);

    run_opt.workdir = opt.workdir;
    run_opt.out_dir = opt.out_dir;
    run_opt.geo_path = opt.geo_path;
    run_opt.beam_path = opt.beam_path;
    run_opt.mat_path = opt.mat_path;
    run_opt.detect_path = opt.detect_path;
    run_opt.nstat = opt.nstat;
    run_opt.has_nstat = opt.has_nstat;
    run_opt.validate_only = opt.dry_run ? 1 : 0;

    rc = osh_run(&run_opt, stdout, stderr);
    osh_log_close();
    return exit_code_for_status(rc);
}

static int exit_code_for_status(enum osh_status status) {
    switch (status) {
    case OSH_OK:
        return EX_OK;
    case OSH_EINVAL:
        return EX_USAGE;
    case OSH_EIO:
        return EX_NOINPUT;
    case OSH_EPARSE:
    case OSH_EINCOMPLETE:
        return EX_CONFIG;
    default:
        return 1;
    }
}
