#include <stdio.h>
#include <stdlib.h>

#include "cli/osh_cli.h"
#include "openshieldhit/openshieldhit.h"

enum osh_main_exit { OSH_EXIT_OK = 0, OSH_EXIT_FAIL = 1, OSH_EXIT_USAGE = 2 };

static int exit_code_for_status(enum openshieldhit_status status);
static enum openshieldhit_status apply_options(openshieldhit_context_t *ctx, struct osh_cli_options const *opt);

int main(int argc, char *argv[]) {
    struct osh_cli_options opt;
    openshieldhit_context_t *ctx;
    char err[256];
    enum openshieldhit_status rc;

    if (osh_cli_parse(argc, argv, &opt, err, sizeof(err)) != 0) {
        fprintf(stderr, "Error: %s\n\n", err);
        osh_cli_print_help(stderr, argv[0]);
        return OSH_EXIT_USAGE;
    }

    if (opt.action == OSH_CLI_ACTION_HELP) {
        osh_cli_print_help(stdout, argv[0]);
        return OSH_EXIT_OK;
    }

    if (opt.action == OSH_CLI_ACTION_VERSION) {
        printf("OpenShieldHIT version %s\n", openshieldhit_version_string());
        return OSH_EXIT_OK;
    }

    ctx = openshieldhit_context_create();
    if (!ctx) {
        fprintf(stderr, "Error: out of memory\n");
        return OSH_EXIT_FAIL;
    }

    rc = apply_options(ctx, &opt);
    if (rc != OPENSHIELDHIT_STATUS_OK) {
        fprintf(stderr, "Error: out of memory while configuring context\n");
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
        return OSH_EXIT_OK;
    case OPENSHIELDHIT_STATUS_INVALID_ARGUMENT:
        return OSH_EXIT_USAGE;
    default:
        return OSH_EXIT_FAIL;
    }
}

static enum openshieldhit_status apply_options(openshieldhit_context_t *ctx, struct osh_cli_options const *opt) {
    enum openshieldhit_status rc;

#define TRY(call)                                                                                                      \
    do {                                                                                                               \
        rc = (call);                                                                                                   \
        if (rc != OPENSHIELDHIT_STATUS_OK) {                                                                           \
            return rc;                                                                                                 \
        }                                                                                                              \
    } while (0)

    TRY(openshieldhit_context_set_workdir(ctx, opt->workdir));
    TRY(openshieldhit_context_set_out_dir(ctx, opt->out_dir));
    TRY(openshieldhit_context_set_log_level(ctx, opt->verbose));
    TRY(openshieldhit_context_set_geo_path(ctx, opt->geo_path));
    TRY(openshieldhit_context_set_beam_path(ctx, opt->beam_path));
    TRY(openshieldhit_context_set_mat_path(ctx, opt->mat_path));
    TRY(openshieldhit_context_set_detect_path(ctx, opt->detect_path));
    if (opt->has_nstat) {
        TRY(openshieldhit_context_set_nstat(ctx, opt->nstat));
    }
    if (opt->dry_run) {
        TRY(openshieldhit_context_set_run_mode(ctx, OPENSHIELDHIT_RUN_VALIDATE));
    }

#undef TRY

    return OPENSHIELDHIT_STATUS_OK;
}
