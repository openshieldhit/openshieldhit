#include "cli/osh_cli.h"

#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/osh_rc.h"

enum osh_cli_long_opt_id {
    OSH_CLI_OPT_DRY_RUN = 1000,
    OSH_CLI_OPT_WORKDIR,
    OSH_CLI_OPT_GEO,
    OSH_CLI_OPT_BEAM,
    OSH_CLI_OPT_MAT,
    OSH_CLI_OPT_DETECT,
    OSH_CLI_OPT_OUTPUT
};

static char const *const OSH_CLI_SHORT_OPTS = "hVvN:o:";

static int set_err(char *err, size_t err_cap, char const *fmt, char const *arg);
static int parse_u64(char const *s, unsigned long long *out);

int osh_cli_parse(int argc, char *argv[], struct osh_cli_options *opt, char *err, size_t err_cap) {
    int c;
    int long_idx = 0;
    static struct option const long_opts[] = {
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {"verbose", no_argument, NULL, 'v'},
        {"nstat", required_argument, NULL, 'N'},
        {"dry-run", no_argument, NULL, OSH_CLI_OPT_DRY_RUN},
        {"workdir", required_argument, NULL, OSH_CLI_OPT_WORKDIR},
        {"geo", required_argument, NULL, OSH_CLI_OPT_GEO},
        {"beam", required_argument, NULL, OSH_CLI_OPT_BEAM},
        {"mat", required_argument, NULL, OSH_CLI_OPT_MAT},
        {"detect", required_argument, NULL, OSH_CLI_OPT_DETECT},
        {"outdir", required_argument, NULL, 'o'},
        {"out", required_argument, NULL, 'o'},
        {"output", required_argument, NULL, OSH_CLI_OPT_OUTPUT},
        {NULL, 0, NULL, 0}
    };

    if (!opt) {
        return OSH_EINVAL;
    }

    opt->action = OSH_CLI_ACTION_RUN;
    opt->dry_run = 0;
    opt->verbose = 0;
    opt->workdir = NULL;
    opt->geo_path = NULL;
    opt->beam_path = NULL;
    opt->mat_path = NULL;
    opt->detect_path = NULL;
    opt->out_dir = NULL;
    opt->nstat = 0;
    opt->has_nstat = 0;

    if (argc <= 1) {
        opt->action = OSH_CLI_ACTION_HELP;
        return OSH_OK;
    }

    /* Reset getopt state in case parser is called multiple times in-process. */
    optind = 1;
    opterr = 0;
    optopt = 0;

    while ((c = getopt_long(argc, argv, OSH_CLI_SHORT_OPTS, long_opts, &long_idx)) != -1) {
        switch (c) {
        case 'h':
            opt->action = OSH_CLI_ACTION_HELP;
            break;
        case 'V':
            opt->action = OSH_CLI_ACTION_VERSION;
            break;
        case 'v':
            opt->verbose += 1;
            break;
        case 'N':
            if (!parse_u64(optarg, &opt->nstat)) {
                return set_err(err, err_cap, "invalid integer value for option '%s'", "-N/--nstat");
            }
            opt->has_nstat = 1;
            break;
        case 'o':
            opt->out_dir = optarg;
            break;
        case OSH_CLI_OPT_DRY_RUN:
            opt->dry_run = 1;
            break;
        case OSH_CLI_OPT_WORKDIR:
            opt->workdir = optarg;
            break;
        case OSH_CLI_OPT_GEO:
            opt->geo_path = optarg;
            break;
        case OSH_CLI_OPT_BEAM:
            opt->beam_path = optarg;
            break;
        case OSH_CLI_OPT_MAT:
            opt->mat_path = optarg;
            break;
        case OSH_CLI_OPT_DETECT:
            opt->detect_path = optarg;
            break;
        case OSH_CLI_OPT_OUTPUT:
            opt->out_dir = optarg;
            break;
        case '?':
        default:
            if ((optind > 0) && argv[optind - 1]) {
                return set_err(err, err_cap, "unknown or invalid option '%s'", argv[optind - 1]);
            }
            return set_err(err, err_cap, "unknown or invalid option '%s'", "(unknown)");
        }
    }

    while (optind < argc) {
        char const *arg = argv[optind++];
        if (!opt->workdir) {
            opt->workdir = arg;
        } else {
            return set_err(err, err_cap, "unexpected positional argument '%s'", arg);
        }
    }

    return OSH_OK;
}

void osh_cli_print_help(FILE *out, char const *prog) {
    char const *cmd = (prog && *prog) ? prog : "openshieldhit";

    fprintf(out, "Usage: %s [OPTIONS] [WORKDIR]\n", cmd);
    fprintf(out, "OpenShieldHIT - Monte Carlo Particle Transport\n\n");
    fprintf(out, "Options:\n");
    fprintf(out, "  -h, --help            Show this help message\n");
    fprintf(out, "  -V, --version         Print version information\n");
    fprintf(out, "  -v, --verbose         Increase verbosity\n");
    fprintf(out, "  -N, --nstat <n>       Number of requested primary histories\n");
    fprintf(out, "      --dry-run         Parse/load inputs only, do not run transport\n");
    fprintf(out, "      --workdir <dir>   Working directory for default input/output files\n");
    fprintf(out, "      --geo <file>      Override geometry input file\n");
    fprintf(out, "      --beam <file>     Override beam input file\n");
    fprintf(out, "      --mat <file>      Override material input file\n");
    fprintf(out, "      --detect <file>   Override scoring input file\n");
    fprintf(out, "  -o, --outdir <dir>    Override output directory\n");
    fprintf(out, "\n");
    fprintf(out, "Notes:\n");
    fprintf(out, "  WORKDIR defaults input files to WORKDIR/{geo,beam,mat,detect}.dat.\n");
    fprintf(out, "  Input overrides replace only the corresponding default file.\n");
}

static int set_err(char *err, size_t err_cap, char const *fmt, char const *arg) {
    if (err && (err_cap > 0)) {
        (void) snprintf(err, err_cap, fmt, arg);
        err[err_cap - 1] = '\0';
    }
    return OSH_EINVAL;
}

static int parse_u64(char const *s, unsigned long long *out) {
    unsigned long long v;
    char *end = NULL;

    if (!s || !*s || !out) {
        return 0;
    }

    errno = 0;
    v = strtoull(s, &end, 10);
    if ((errno != 0) || (end == s) || (end && *end != '\0')) {
        return 0;
    }

    *out = v;
    return 1;
}
