#ifndef OPENSHIELDHIT_OPENSHIELDHIT_H
#define OPENSHIELDHIT_OPENSHIELDHIT_H

#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Public API for libopenshieldhit
 * ================================
 *
 * DESIGN RULE: No internal structs or types may appear in this header.
 *
 * All mutable state is held in the opaque handle `openshieldhit_context_t`.
 * Callers create a context, configure it through one config struct, run it,
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
 *  - A config struct plus a configuration function.
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
 * Add a field to `openshieldhit_config_t` when source compatibility is enough.
 * If true binary compatibility becomes important later, evolve the config
 * struct through its `size` field.
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

/* ---- Configuration ----
 *
 * Initialize with OPENSHIELDHIT_CONFIG_INIT, then override only the fields
 * you care about. Zero values select defaults.
 *
 * The `size` field lets newer libraries distinguish which fields are present
 * in older caller-compiled structs. Set it with OPENSHIELDHIT_CONFIG_INIT.
 */
typedef struct openshieldhit_config {
    size_t size;
    char const *workdir;
    char const *out_dir;
    char const *geo_path;
    char const *beam_path;
    char const *mat_path;
    char const *detect_path;
    enum openshieldhit_run_mode run_mode;
    int log_level;
    unsigned long long nstat;
    int has_nstat;
} openshieldhit_config_t;

#define OPENSHIELDHIT_CONFIG_INIT                                                                                      \
    {sizeof(openshieldhit_config_t), NULL, NULL, NULL, NULL, NULL, NULL, OPENSHIELDHIT_RUN_NORMAL, 0, 0ULL, 0}

/* Lifecycle */
openshieldhit_context_t *openshieldhit_context_create(void);
void openshieldhit_context_destroy(openshieldhit_context_t *ctx);

/* Configure
 *
 * String arguments are deep-copied into the context. The caller does not need
 * to keep the original strings alive after this call returns.
 *
 * Passing NULL or zero values selects defaults for the corresponding field.
 * Future extension points can add fields to openshieldhit_config_t while older
 * callers remain source-compatible through the `size` field.
 */
enum openshieldhit_status openshieldhit_context_configure(openshieldhit_context_t *ctx,
                                                          openshieldhit_config_t const *cfg);

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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* OPENSHIELDHIT_OPENSHIELDHIT_H */
