#include "cli/osh_cli.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_duration.h"

/* Keep CLI-compatible behavior with original SHIELD-HIT seed offset bounds. */
#define OSH_CLI_MAX_SEED_OFFSET 9999ull

/* This parser intentionally does not use getopt/getopt_long because MSVC
 * does not provide <getopt.h> in its standard C runtime. Keeping the option
 * handling local avoids an extra compatibility dependency on Windows. */

static int set_err(char *err, size_t err_cap, char const *fmt, char const *arg);
static int parse_u64(char const *s, unsigned long long *out);
static int parse_long_option(int argc, char *argv[], int *idx, struct osh_cli_options *opt, char *err, size_t err_cap);
static int
parse_short_options(int argc, char *argv[], int *idx, struct osh_cli_options *opt, char *err, size_t err_cap);
static int consume_option_arg(int argc, char *argv[], int *idx, char const *current, char const **value_out);

int osh_cli_parse(int argc, char *argv[], struct osh_cli_options *opt, char *err, size_t err_cap) {
    int i;
    int positional_only = 0;

    if (!opt) {
        return 1;
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
    opt->seed_offset = 0;
    opt->has_seed_offset = 0;
    opt->pool_capacity = 0;
    opt->has_pool_capacity = 0;
    opt->mem_budget = NULL;
    opt->profile_path = NULL;
    opt->max_time_s = 0.0;
    opt->has_max_time = 0;

    if (argc <= 1) {
        opt->action = OSH_CLI_ACTION_HELP;
        return 0;
    }

    for (i = 1; i < argc; ++i) {
        char const *arg = argv[i];

        if (!arg) {
            continue;
        }

        if (!positional_only && strcmp(arg, "--") == 0) {
            positional_only = 1;
            continue;
        }

        if (!positional_only && (strncmp(arg, "--", 2) == 0)) {
            if (parse_long_option(argc, argv, &i, opt, err, err_cap) != 0) {
                return 1;
            }
            continue;
        }

        if (!positional_only && (arg[0] == '-') && arg[1] != '\0') {
            if (parse_short_options(argc, argv, &i, opt, err, err_cap) != 0) {
                return 1;
            }
        } else {
            if (!opt->workdir) {
                opt->workdir = arg;
            } else {
                return set_err(err, err_cap, "unexpected positional argument '%s'", arg);
            }
        }
    }

    return 0;
}

void osh_cli_print_help(FILE *out, char const *prog) {
    char const *cmd = (prog && *prog) ? prog : "openshieldhit";

    fprintf(out, "Usage: %s [OPTIONS] [WORKDIR]\n", cmd);
    fprintf(out, "OpenShieldHIT - Monte Carlo Particle Transport\n\n");
    fprintf(out, "Options:\n");
    fprintf(out, "  -h, --help            Show this help message\n");
    fprintf(out, "  -V, --version         Print version information\n");
    fprintf(out, "      --print-resources Detect and print host CPU/RAM, then exit\n");
    fprintf(out, "  -v, --verbose         Increase verbosity\n");
    fprintf(out, "  -n, --nstat <n>       Number of requested primary histories\n");
    fprintf(out, "  -N, --seedoffset <n>  Random seed offset override (max 9999)\n");
    fprintf(out, "      --dry-run         Parse/load inputs only, do not run transport\n");
    fprintf(out, "      --workdir <dir>   Working directory for default input/output files\n");
    fprintf(out, "  -g, --geo <file>      Override geometry input file\n");
    fprintf(out, "  -b, --beam <file>     Override beam input file\n");
    fprintf(out, "  -m, --mat <file>      Override material input file\n");
    fprintf(out, "  -d, --detect <file>   Override scoring input file\n");
    fprintf(out, "  -o, --outdir <dir>    Override output directory\n");
    fprintf(out,
            "      --pool-capacity <n>  Live-history pool size (perf knob; 0 = default; physics unchanged,\n"
            "                           scored output matches up to FP rounding)\n");
    fprintf(out,
            "      --mem-budget <size>  Max memory OSH may use for scoring (e.g. 8GB, 512MiB, 80%%);\n"
            "                           a run whose scoring exceeds this is refused. Default: 80%% of\n"
            "                           available RAM\n");
    fprintf(out, "      --profile <file>  Write a one-line JSON timing/counter profile to <file>\n");
    fprintf(out,
            "      --max-time <dur>  Wall-time budget; stop cleanly at the next safe point and save the\n"
            "                        partial result (e.g. 30s, 30m, 1h, or a bare number of seconds).\n"
            "                        Overrides the beam.dat MAXTIME card.\n");
    fprintf(out, "\n");
    fprintf(out, "Notes:\n");
    fprintf(out, "  WORKDIR defaults input files to WORKDIR/{geo,beam,mat,detect}.dat.\n");
    fprintf(out, "  Input overrides replace only the corresponding default file.\n");
    fprintf(out, "  Ctrl-C (SIGINT) requests the same clean stop: in-flight histories finish,\n");
    fprintf(out, "  all secondaries drain, and the partial result is saved normalised by the\n");
    fprintf(out, "  true number of completed primaries.\n");
}

/**
 * @brief Write a formatted error message into err and return 1.
 *
 * @details
 * Convenience so callers can write `return set_err(...)` in one statement.
 * fmt must contain exactly one %s placeholder, which is filled by arg.
 * Always returns 1 regardless of whether err is NULL.
 *
 * @param[out] err      Destination buffer; silently ignored when NULL.
 * @param[in]  err_cap  Capacity of err in bytes.
 * @param[in]  fmt      printf-style format string with exactly one %s.
 * @param[in]  arg      String substituted for %s.
 *
 * @returns 1 always, so callers can `return set_err(...)`.
 */
static int set_err(char *err, size_t err_cap, char const *fmt, char const *arg) {
    if (err && (err_cap > 0)) {
        (void) snprintf(err, err_cap, fmt, arg);
        err[err_cap - 1] = '\0';
    }
    return 1;
}

/**
 * @brief Parse a decimal unsigned 64-bit integer from a string.
 *
 * @details
 * Rejects empty strings, leading whitespace, sign prefixes, and any trailing
 * non-digit characters. Uses strtoull with base 10.
 *
 * @param[in]  s    Null-terminated string to parse.
 * @param[out] out  Receives the parsed value on success.
 *
 * @returns 1 on success, 0 if the string is invalid or out-of-range.
 */
static int parse_u64(char const *s, unsigned long long *out) {
    unsigned long long v;
    char *end = NULL;

    if (!s || !*s || !out) {
        return 0;
    }
    if (s[0] < '0' || s[0] > '9') {
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

/**
 * @brief Parse a single long option from argv[*idx].
 *
 * @details
 * Handles both `--name=value` and `--name value` forms. After a successful
 * parse, *idx points at the last token consumed; the main loop will advance
 * it once more before the next iteration.
 *
 * @param[in]     argc     Total argument count.
 * @param[in]     argv     Argument vector.
 * @param[in,out] idx      Index of the current `--...` token; may be
 *                         advanced when the value is a separate token.
 * @param[in,out] opt      Options struct to update.
 * @param[out]    err      Error message buffer.
 * @param[in]     err_cap  Capacity of err in bytes.
 *
 * @returns 0 on success, 1 on error.
 */
static int parse_long_option(int argc, char *argv[], int *idx, struct osh_cli_options *opt, char *err, size_t err_cap) {
    char *eq;
    char const *arg = argv[*idx];
    char const *name = arg + 2;
    char const *value = NULL;
    size_t name_len;

    eq = strchr(name, '=');
    if (eq) {
        name_len = (size_t) (eq - name);
        value = eq + 1;
    } else {
        name_len = strlen(name);
    }

    if ((name_len == 4) && (strncmp(name, "help", name_len) == 0) && !value) {
        opt->action = OSH_CLI_ACTION_HELP;
        return 0;
    }
    if ((name_len == 7) && (strncmp(name, "version", name_len) == 0) && !value) {
        opt->action = OSH_CLI_ACTION_VERSION;
        return 0;
    }
    if ((name_len == 15) && (strncmp(name, "print-resources", name_len) == 0) && !value) {
        opt->action = OSH_CLI_ACTION_PRINT_RESOURCES;
        return 0;
    }
    if ((name_len == 7) && (strncmp(name, "verbose", name_len) == 0) && !value) {
        opt->verbose += 1;
        return 0;
    }
    if ((name_len == 7) && (strncmp(name, "dry-run", name_len) == 0) && !value) {
        opt->dry_run = 1;
        return 0;
    }

    if ((name_len == 5) && (strncmp(name, "nstat", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        if (!parse_u64(value, &opt->nstat)) {
            return set_err(err, err_cap, "invalid integer value for option '%s'", "-n/--nstat");
        }
        opt->has_nstat = 1;
        return 0;
    }
    if (((name_len == 10) && (strncmp(name, "seedoffset", name_len) == 0))
        || ((name_len == 11) && (strncmp(name, "seed-offset", name_len) == 0))) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        if (!parse_u64(value, &opt->seed_offset) || (opt->seed_offset > OSH_CLI_MAX_SEED_OFFSET)) {
            return set_err(err, err_cap, "invalid integer value for option '%s'", "-N/--seedoffset");
        }
        opt->has_seed_offset = 1;
        return 0;
    }
    if ((name_len == 13) && (strncmp(name, "pool-capacity", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        /* 0 is accepted and means "use the compiled default", matching the
         * pool_capacity semantics in the transport params and public API. */
        if (!parse_u64(value, &opt->pool_capacity)) {
            return set_err(err, err_cap, "invalid integer value for option '%s'", "--pool-capacity");
        }
        opt->has_pool_capacity = 1;
        return 0;
    }
    if ((name_len == 10) && (strncmp(name, "mem-budget", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->mem_budget = value;
        return 0;
    }
    if ((name_len == 8) && (strncmp(name, "max-time", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        if (!osh_parse_duration(value, &opt->max_time_s)) {
            return set_err(err, err_cap, "invalid duration value for option '%s' (use e.g. 30s, 30m, 1h, or 500)", arg);
        }
        opt->has_max_time = 1;
        return 0;
    }
    if ((name_len == 7) && (strncmp(name, "workdir", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->workdir = value;
        return 0;
    }
    if ((name_len == 3) && (strncmp(name, "geo", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->geo_path = value;
        return 0;
    }
    if ((name_len == 4) && (strncmp(name, "beam", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->beam_path = value;
        return 0;
    }
    if ((name_len == 3) && (strncmp(name, "mat", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->mat_path = value;
        return 0;
    }
    if ((name_len == 6) && (strncmp(name, "detect", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->detect_path = value;
        return 0;
    }
    if ((name_len == 7) && (strncmp(name, "profile", name_len) == 0)) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->profile_path = value;
        return 0;
    }
    if (((name_len == 6) && (strncmp(name, "outdir", name_len) == 0))
        || ((name_len == 3) && (strncmp(name, "out", name_len) == 0))
        || ((name_len == 6) && (strncmp(name, "output", name_len) == 0))) {
        if (!value && !consume_option_arg(argc, argv, idx, arg, &value)) {
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
        opt->out_dir = value;
        return 0;
    }

    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
}

/**
 * @brief Parse all short options packed in argv[*idx] (e.g. `-vvn 1000`).
 *
 * @details
 * Iterates over each character in the flag cluster. Flags that take a value
 * (-n, -N, -o) consume either the remainder of the current token or the next
 * argv token, then return immediately since no further flags can follow.
 *
 * @param[in]     argc     Total argument count.
 * @param[in]     argv     Argument vector.
 * @param[in,out] idx      Index of the current `-...` token; may be advanced.
 * @param[in,out] opt      Options struct to update.
 * @param[out]    err      Error message buffer.
 * @param[in]     err_cap  Capacity of err in bytes.
 *
 * @returns 0 on success, 1 on error.
 */
static int
parse_short_options(int argc, char *argv[], int *idx, struct osh_cli_options *opt, char *err, size_t err_cap) {
    char const *arg = argv[*idx];
    size_t pos;

    for (pos = 1; arg[pos] != '\0'; ++pos) {
        char const *value = NULL;
        char ch = arg[pos];

        switch (ch) {
        case 'h':
            opt->action = OSH_CLI_ACTION_HELP;
            break;
        case 'V':
            opt->action = OSH_CLI_ACTION_VERSION;
            break;
        case 'v':
            opt->verbose += 1;
            break;
        case 'n':
            value = &arg[pos + 1];
            if (*value == '\0') {
                if (!consume_option_arg(argc, argv, idx, arg, &value)) {
                    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
                }
            }
            if (!parse_u64(value, &opt->nstat)) {
                return set_err(err, err_cap, "invalid integer value for option '%s'", "-n/--nstat");
            }
            opt->has_nstat = 1;
            return 0;
        case 'N':
            value = &arg[pos + 1];
            if (*value == '\0') {
                if (!consume_option_arg(argc, argv, idx, arg, &value)) {
                    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
                }
            }
            if (!parse_u64(value, &opt->seed_offset) || (opt->seed_offset > OSH_CLI_MAX_SEED_OFFSET)) {
                return set_err(err, err_cap, "invalid integer value for option '%s'", "-N/--seedoffset");
            }
            opt->has_seed_offset = 1;
            return 0;
        case 'b':
            value = &arg[pos + 1];
            if (*value == '\0') {
                if (!consume_option_arg(argc, argv, idx, arg, &value)) {
                    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
                }
            }
            opt->beam_path = value;
            return 0;
        case 'g':
            value = &arg[pos + 1];
            if (*value == '\0') {
                if (!consume_option_arg(argc, argv, idx, arg, &value)) {
                    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
                }
            }
            opt->geo_path = value;
            return 0;
        case 'm':
            value = &arg[pos + 1];
            if (*value == '\0') {
                if (!consume_option_arg(argc, argv, idx, arg, &value)) {
                    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
                }
            }
            opt->mat_path = value;
            return 0;
        case 'd':
            value = &arg[pos + 1];
            if (*value == '\0') {
                if (!consume_option_arg(argc, argv, idx, arg, &value)) {
                    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
                }
            }
            opt->detect_path = value;
            return 0;
        case 'o':
            value = &arg[pos + 1];
            if (*value == '\0') {
                if (!consume_option_arg(argc, argv, idx, arg, &value)) {
                    return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
                }
            }
            opt->out_dir = value;
            return 0;
        default:
            return set_err(err, err_cap, "unknown or invalid option '%s'", arg);
        }
    }

    return 0;
}

/**
 * @brief Advance the argument index and return the next argv token as a value.
 *
 * @details
 * Used when an option's value is a separate token (e.g. `--geo file.dat`
 * rather than `--geo=file.dat`). Increments *idx so the main loop does not
 * consume the value token a second time.
 *
 * @param[in]     argc       Total argument count.
 * @param[in]     argv       Argument vector.
 * @param[in,out] idx        Current index; incremented to point at the value.
 * @param[in]     current    The option token requiring a value (used only to
 *                           guard against a NULL or empty current argument).
 * @param[out]    value_out  Receives a pointer into argv[*idx].
 *
 * @returns 1 if a value token was found, 0 if argv is exhausted.
 */
static int consume_option_arg(int argc, char *argv[], int *idx, char const *current, char const **value_out) {
    if ((current != NULL) && (current[0] != '\0')) {
        ++(*idx);
    }
    if ((*idx >= argc) || (argv[*idx] == NULL)) {
        return 0;
    }
    *value_out = argv[*idx];
    return 1;
}
