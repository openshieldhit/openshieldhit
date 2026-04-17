#ifndef OSH_CLI_H
#define OSH_CLI_H

#include <stddef.h>
#include <stdio.h>

/* Internal CLI parser — not part of the public library API.
 * This module is the project's cross-platform getopt/getopt_long replacement.
 * Used by src/apps/osh/main.c and the test_osh_cli test suite. */

/** Actions the CLI parser may request the caller to perform. */
enum osh_cli_action { OSH_CLI_ACTION_RUN = 0, OSH_CLI_ACTION_HELP = 1, OSH_CLI_ACTION_VERSION = 2 };

/**
 * @brief Parsed and validated command-line options.
 *
 * @details
 * Populated by osh_cli_parse(). All pointer fields alias into the original
 * argv array and are valid only for its lifetime.
 */
struct osh_cli_options {
    enum osh_cli_action action; /**< Requested action (run / help / version). */
    int dry_run;                /**< Non-zero: load inputs but skip transport. */
    int verbose;                /**< Verbosity level; incremented per -v flag. */
    char const *workdir;        /**< Working directory for default file resolution. */
    char const *geo_path;       /**< Override path for the geometry input file. */
    char const *beam_path;      /**< Override path for the beam input file. */
    char const *mat_path;       /**< Override path for the material input file. */
    char const *detect_path;    /**< Override path for the scoring input file. */
    char const *out_dir;        /**< Override path for the output directory. */
    unsigned long long nstat;   /**< Requested number of primary histories. */
    int has_nstat;              /**< Non-zero if --nstat/-n was explicitly given. */
    unsigned long long seed_offset; /**< Random-seed stream offset override. */
    int has_seed_offset;            /**< Non-zero if --seedoffset/-N was explicitly given. */
};

/**
 * @brief Parse command-line arguments into an options struct.
 *
 * @details
 * Iterates over argv[1..argc-1] and fills *opt. Supports '--' to stop option
 * processing. If argc <= 1 (no arguments given), sets action to
 * OSH_CLI_ACTION_HELP so the caller can print usage without treating it as
 * an error. Pointer fields in *opt alias into argv.
 *
 * @param[in]     argc     Argument count from main().
 * @param[in]     argv     Argument vector from main().
 * @param[out]    opt      Receives the parsed options; must not be NULL.
 * @param[out]    err      Buffer for a human-readable error message, or NULL.
 * @param[in]     err_cap  Capacity of err in bytes.
 *
 * @returns 0 on success, 1 on parse error (err is populated when non-NULL).
 */
int osh_cli_parse(int argc, char *argv[], struct osh_cli_options *opt, char *err, size_t err_cap);

/**
 * @brief Print a usage summary to the given stream.
 *
 * @param[in] out   Destination stream (e.g. stdout or stderr).
 * @param[in] prog  argv[0] used as the program name; falls back to
 *                  "openshieldhit" when NULL or empty.
 */
void osh_cli_print_help(FILE *out, char const *prog);

#endif /* OSH_CLI_H */
