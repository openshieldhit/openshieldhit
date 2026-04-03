#include <stdio.h>
#include <stdlib.h>

#include "cli/osh_cli.h"
#include "common/osh_exit.h"
#include "openshieldhit/openshieldhit.h"

/* Generic failure — not covered by a specific EX_* code. */
#define OSH_EXIT_FAIL 1

static int exit_code_for_status(enum openshieldhit_status status);

int main(int argc, char *argv[]) {
    struct osh_cli_options opt;
    openshieldhit_config_t cfg = OPENSHIELDHIT_CONFIG_INIT;
    openshieldhit_context_t *ctx;
    char err[256];
    enum openshieldhit_status rc;

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
        printf("OpenShieldHIT version %s\n", openshieldhit_version_string());
        return EX_OK;
    }

    ctx = openshieldhit_context_create();
    if (!ctx) {
        fprintf(stderr, "Error: out of memory\n");
        return OSH_EXIT_FAIL;
    }

    cfg.workdir = opt.workdir;
    cfg.out_dir = opt.out_dir;
    cfg.geo_path = opt.geo_path;
    cfg.beam_path = opt.beam_path;
    cfg.mat_path = opt.mat_path;
    cfg.detect_path = opt.detect_path;
    cfg.run_mode = opt.dry_run ? OPENSHIELDHIT_RUN_VALIDATE : OPENSHIELDHIT_RUN_NORMAL;
    cfg.log_level = opt.verbose;
    cfg.nstat = opt.nstat;
    cfg.has_nstat = opt.has_nstat;

    rc = openshieldhit_context_configure(ctx, &cfg);
    if (rc != OPENSHIELDHIT_STATUS_OK) {
        fprintf(stderr, "Error: failed to configure context\n");
        openshieldhit_context_destroy(ctx);
        return OSH_EXIT_FAIL;
    }

    rc = openshieldhit_run(ctx, stdout, stderr);
    openshieldhit_context_destroy(ctx);
    return exit_code_for_status(rc);
}

static int exit_code_for_status(enum openshieldhit_status status) {
    switch (status) {
    case OPENSHIELDHIT_STATUS_OK:
        return EX_OK;
    case OPENSHIELDHIT_STATUS_INVALID_ARGUMENT:
        return EX_USAGE;
    case OPENSHIELDHIT_STATUS_IO_ERROR:
        return EX_NOINPUT;
    case OPENSHIELDHIT_STATUS_PARSE_ERROR:
    case OPENSHIELDHIT_STATUS_INCOMPLETE:
        return EX_CONFIG;
    default:
        return OSH_EXIT_FAIL;
    }
}
