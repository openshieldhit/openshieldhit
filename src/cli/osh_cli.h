#ifndef OSH_CLI_H
#define OSH_CLI_H

#include <stddef.h>
#include <stdio.h>

enum osh_cli_action {
    OSH_CLI_ACTION_RUN = 0,
    OSH_CLI_ACTION_HELP = 1,
    OSH_CLI_ACTION_VERSION = 2
};

struct osh_cli_options {
    enum osh_cli_action action;
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

int osh_cli_parse(int argc, char *argv[], struct osh_cli_options *opt, char *err, size_t err_cap);
void osh_cli_print_help(FILE *out, char const *prog);

#endif /* OSH_CLI_H */
