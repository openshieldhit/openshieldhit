#include <stdio.h>
#include <stdlib.h>

#include "apps/osh/osh_run.h"
#include "cli/osh_cli.h"
#include "common/osh_exit.h"
#include "common/osh_file.h"
#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"
#include "openshieldhit/version.h"

/**
 * @brief Map library status codes to conventional process exit codes.
 *
 * @details
 * This keeps CLI-facing shell behavior stable even if the internal status enum
 * grows. Only a small subset currently gets dedicated sysexits mappings; all
 * other failures fall back to a generic non-zero exit status.
 */
static int exit_code_for_status(enum osh_status status);
static void cli_diag_emit(void *user, int level, char const *file, int line, char const *function, char const *message);

int main(int argc, char *argv[]) {
    struct osh_cli_options opt;
    struct osh_run_options run_opt;
    struct osh_diag_sink diag;
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

    diag.emit = cli_diag_emit;
    diag.min_level = (opt.verbose == 0)   ? OSH_DIAG_LEVEL_WARN
                     : (opt.verbose == 1) ? OSH_DIAG_LEVEL_INFO
                                          : OSH_DIAG_LEVEL_DEBUG;

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
    run_opt.seed_offset = opt.seed_offset;
    run_opt.has_seed_offset = opt.has_seed_offset;
    run_opt.validate_only = opt.dry_run ? 1 : 0;
    run_opt.diag = &diag;

    rc = osh_run(&run_opt, stdout, stderr);
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

static void
cli_diag_emit(void *user, int level, char const *file, int line, char const *function, char const *message) {
    FILE *fp;
    (void) user;
    (void) file;
    (void) line;
    (void) function;

    if (!message) {
        return;
    }

    /* Keep informational chatter pipe-friendly on stdout while routing
     * warnings and errors to stderr. To add log-file output, replace the
     * NULL user pointer with a small context struct carrying a log FILE*,
     * then write to it here alongside fp — no API change required. */
    fp = (level >= OSH_DIAG_LEVEL_WARN) ? stderr : stdout;

    fprintf(fp, "[%s] %s\n", osh_diag_level_name(level), message);
    fflush(fp);
}
