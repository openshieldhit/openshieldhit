#ifndef OPENSHIELDHIT_OPENSHIELDHIT_H
#define OPENSHIELDHIT_OPENSHIELDHIT_H

#include <stddef.h>
#include <stdio.h>

/*
 * Public API for libopenshieldhit
 * ================================
 *
 * DESIGN RULE: No internal structs or types may appear in this header.
 *
 * All mutable state is held in the opaque handle `openshieldhit_context_t`.
 * Callers create a context, configure it through setter functions, run it,
 * and destroy it. They never access fields directly.
 *
 * Why opaque handles?
 * -------------------
 * Exposing a struct in a public header freezes its layout as part of the ABI.
 * Every field addition or reordering breaks binary compatibility with code
 * compiled against the old header. Opaque handles avoid this entirely: the
 * struct definition lives in openshieldhit.c and is invisible to consumers.
 *
 * This also keeps the door open for alternative configuration sources without
 * any API change — a CLI parser, a JSON string, a config file, or a pipe all
 * become different ways to call the same setter functions.
 *
 * What belongs here?
 * ------------------
 *  - The opaque typedef and its lifecycle functions (create / destroy).
 *  - Setter functions for configuration.
 *  - The run function and error retrieval.
 *  - Version query functions.
 *
 * What does NOT belong here?
 * --------------------------
 *  - Internal structs or enums (keep them in src/).
 *  - CLI-specific types (osh_cli_options lives in src/cli/osh_cli.h).
 *  - Implementation details of any kind.
 *
 * Adding new functionality?
 * -------------------------
 * Add a setter function, not a new field in a struct. If you find yourself
 * wanting to put a struct in this header, stop and reconsider.
 */

/* ---- Status codes ---- */
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

/* ---- Run mode ----
 *
 * Controls what openshieldhit_run() actually does.
 * Prefer this over boolean flags — it can be extended without breaking callers.
 */
enum openshieldhit_run_mode {
    OPENSHIELDHIT_RUN_NORMAL = 0,  /* full simulation — returns NOT_SUPPORTED until implemented */
    OPENSHIELDHIT_RUN_VALIDATE = 1 /* parse and validate inputs only; no transport */
};

/* ---- Opaque context ----
 *
 * All runtime configuration is held inside this handle. Callers never
 * touch its fields directly, which keeps the ABI stable across releases
 * and makes it straightforward to add new configuration sources later
 * (e.g. JSON string, config file, pipe) without changing the public API.
 */
typedef struct openshieldhit_context openshieldhit_context_t;

/* Lifecycle */
openshieldhit_context_t *openshieldhit_context_create(void);
void openshieldhit_context_destroy(openshieldhit_context_t *ctx);

/* Configuration setters
 *
 * All setters return OPENSHIELDHIT_STATUS_OK on success or a status code on
 * failure (typically OPENSHIELDHIT_STATUS_NO_MEMORY or
 * OPENSHIELDHIT_STATUS_INVALID_ARGUMENT).
 *
 * String arguments are deep-copied into the context. The caller does not need
 * to keep the original string alive after the setter returns. Passing NULL
 * clears the previously set value.
 *
 * Future extension points (not yet implemented):
 *   openshieldhit_context_configure_from_json(ctx, json_string)
 *   openshieldhit_context_configure_from_file(ctx, path)
 *   openshieldhit_context_configure_from_fd(ctx, fd)   -- pipe / stdin
 */
enum openshieldhit_status openshieldhit_context_set_workdir(openshieldhit_context_t *ctx, char const *path);
enum openshieldhit_status openshieldhit_context_set_out_dir(openshieldhit_context_t *ctx, char const *path);
enum openshieldhit_status openshieldhit_context_set_run_mode(openshieldhit_context_t *ctx,
                                                             enum openshieldhit_run_mode mode);
enum openshieldhit_status openshieldhit_context_set_log_level(openshieldhit_context_t *ctx, int level);
enum openshieldhit_status openshieldhit_context_set_nstat(openshieldhit_context_t *ctx, unsigned long long nstat);
enum openshieldhit_status openshieldhit_context_set_geo_path(openshieldhit_context_t *ctx, char const *path);
enum openshieldhit_status openshieldhit_context_set_beam_path(openshieldhit_context_t *ctx, char const *path);
enum openshieldhit_status openshieldhit_context_set_mat_path(openshieldhit_context_t *ctx, char const *path);
enum openshieldhit_status openshieldhit_context_set_detect_path(openshieldhit_context_t *ctx, char const *path);

/* Run
 *
 * out / err may be NULL, in which case that output channel is suppressed.
 * Errors are always stored in the context and retrievable via
 * openshieldhit_last_error() regardless of whether err is NULL.
 */
enum openshieldhit_status openshieldhit_run(openshieldhit_context_t *ctx, FILE *out, FILE *err);

/* Returns the last error message set by openshieldhit_run(), or an empty
 * string if no error has occurred yet. The pointer is valid until the
 * context is destroyed or the next call to openshieldhit_run(). */
char const *openshieldhit_last_error(openshieldhit_context_t const *ctx);

/* ---- Version ---- */
char const *openshieldhit_version_string(void);
int openshieldhit_version_major(void);
int openshieldhit_version_minor(void);
int openshieldhit_version_patch(void);

#endif /* OPENSHIELDHIT_OPENSHIELDHIT_H */
