#include <stdio.h>
#include <stdlib.h>

#include "openshieldhit/openshieldhit.h"

enum osh_main_exit_code {
    OSH_MAIN_EXIT_OK = 0,
    OSH_MAIN_EXIT_FAIL = 1,
    OSH_MAIN_EXIT_USAGE = 2
};

static size_t const OSH_MAIN_ERR_BUF_CAP = 256;
static int osh_main_exit_code_for_status(int status);

int main(int argc, char *argv[]) {
    int rc;
    char err[OSH_MAIN_ERR_BUF_CAP];
    struct openshieldhit_cli_options opt;

    rc = openshieldhit_cli_parse(argc, argv, &opt, err, OSH_MAIN_ERR_BUF_CAP);
    if (rc != 0) {
        fprintf(stderr, "Error: %s\n\n", err[0] ? err : "invalid command line");
        openshieldhit_cli_print_help(stderr, argv[0]);
        return OSH_MAIN_EXIT_USAGE;
    }

    if (opt.action == OPENSHIELDHIT_CLI_ACTION_HELP) {
        openshieldhit_cli_print_help(stdout, argv[0]);
        return OSH_MAIN_EXIT_OK;
    }

    if (opt.action == OPENSHIELDHIT_CLI_ACTION_VERSION) {
        printf("OpenShieldHIT version %s\n", openshieldhit_version_string());
        return OSH_MAIN_EXIT_OK;
    }

    return osh_main_exit_code_for_status(openshieldhit_run(&opt, stdout, stderr));
}

static int osh_main_exit_code_for_status(int status) {
    if (status == OPENSHIELDHIT_STATUS_OK) {
        return OSH_MAIN_EXIT_OK;
    }
    if (status == OPENSHIELDHIT_STATUS_INVALID_ARGUMENT) {
        return OSH_MAIN_EXIT_USAGE;
    }
    return OSH_MAIN_EXIT_FAIL;
}
