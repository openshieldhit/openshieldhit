#ifndef OPENSHIELDHIT_OPENSHIELDHIT_H
#define OPENSHIELDHIT_OPENSHIELDHIT_H

#include <stddef.h>
#include <stdio.h>

enum openshieldhit_status {
    OPENSHIELDHIT_STATUS_OK = 0,
    OPENSHIELDHIT_STATUS_INVALID_ARGUMENT,
    OPENSHIELDHIT_STATUS_NO_MEMORY,
    OPENSHIELDHIT_STATUS_IO_ERROR,
    OPENSHIELDHIT_STATUS_PARSE_ERROR,
    OPENSHIELDHIT_STATUS_INCOMPLETE,
    OPENSHIELDHIT_STATUS_NOT_SUPPORTED,
    OPENSHIELDHIT_STATUS_STATE_ERROR
};

enum openshieldhit_cli_action {
    OPENSHIELDHIT_CLI_ACTION_RUN = 0,
    OPENSHIELDHIT_CLI_ACTION_HELP = 1,
    OPENSHIELDHIT_CLI_ACTION_VERSION = 2
};

struct openshieldhit_cli_options {
    enum openshieldhit_cli_action action;
    int dry_run;
    int verbose;
    char const *workdir;
    char const *geo_path;
    char const *beam_path;
    char const *mat_path;
    char const *detect_path;
    char const *out_dir;
    unsigned long long nstat;
    int has_nstat;
};

char const *openshieldhit_version_string(void);
int openshieldhit_version_major(void);
int openshieldhit_version_minor(void);
int openshieldhit_version_patch(void);

int openshieldhit_cli_parse(
    int argc, char *argv[], struct openshieldhit_cli_options *opt, char *err, size_t err_cap);
void openshieldhit_cli_print_help(FILE *out, char const *prog);

int openshieldhit_run(struct openshieldhit_cli_options const *opt, FILE *out, FILE *err);

#endif /* OPENSHIELDHIT_OPENSHIELDHIT_H */
