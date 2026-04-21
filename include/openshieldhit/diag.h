#ifndef OPENSHIELDHIT_DIAG_H
#define OPENSHIELDHIT_DIAG_H

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Severity levels for caller-owned diagnostics sinks.
 *
 * @details
 * The library emits diagnostics through a borrowed @ref osh_diag_sink. The
 * caller decides whether messages are printed, logged, buffered, forwarded to a
 * GUI, or ignored entirely.
 */
enum osh_diag_level {
    OSH_DIAG_LEVEL_TRACE = 0,
    OSH_DIAG_LEVEL_DEBUG = 1,
    OSH_DIAG_LEVEL_INFO = 2,
    OSH_DIAG_LEVEL_WARN = 3,
    OSH_DIAG_LEVEL_ERROR = 4,
    OSH_DIAG_LEVEL_FATAL = 5,
    OSH_DIAG_LEVEL_OFF = 6
};

/**
 * @brief Diagnostics callback signature.
 *
 * @param[in] user      Caller-owned context pointer from @ref osh_diag_sink.user.
 * @param[in] level     One of @ref osh_diag_level.
 * @param[in] file      Source file where the diagnostic originated, or NULL.
 * @param[in] line      Source line where the diagnostic originated, or 0.
 * @param[in] function  Source function where the diagnostic originated, or NULL.
 * @param[in] message   Fully formatted human-readable message.
 */
typedef void (*osh_diag_emit_fn)(
    void *user, int level, char const *file, int line, char const *function, char const *message);

/**
 * @brief Caller-owned diagnostics sink configuration.
 *
 * @details
 * The library borrows this struct and never takes ownership of @p user.
 * `emit == NULL` or `min_level == OSH_DIAG_LEVEL_OFF` silences diagnostics.
 */
struct osh_diag_sink {
    osh_diag_emit_fn emit;
    void *user;
    int min_level;
};

/**
 * @brief Return a stable human-readable name for one diagnostics level.
 *
 * @param[in] level  One of @ref osh_diag_level.
 *
 * @returns String literal such as `"INFO"` or `"ERROR"`.
 */
char const *osh_diag_level_name(int level);

/**
 * @brief Format and emit one diagnostics message through a caller-owned sink.
 *
 * @details
 * This helper is intended primarily for library-internal paths that already
 * own a borrowed @ref osh_diag_sink. It checks for NULL sinks and performs
 * level filtering before invoking @ref osh_diag_emit_fn.
 *
 * @param[in] diag      Borrowed diagnostics sink; NULL means silent.
 * @param[in] level     One of @ref osh_diag_level.
 * @param[in] file      Source file of the diagnostic, or NULL.
 * @param[in] line      Source line of the diagnostic, or 0.
 * @param[in] function  Source function of the diagnostic, or NULL.
 * @param[in] fmt       `printf`-style format string.
 * @param[in] ap        Variadic argument list for @p fmt.
 */
void osh_diag_emitfv(struct osh_diag_sink const *diag,
                     int level,
                     char const *file,
                     int line,
                     char const *function,
                     char const *fmt,
                     va_list ap);

/**
 * @brief Convenience wrapper around @ref osh_diag_emitfv.
 *
 * @details
 * Typical caller-side setup is intentionally simple:
 *
 * ```c
 * static void cli_diag(void *user,
 *                      int level,
 *                      char const *file,
 *                      int line,
 *                      char const *function,
 *                      char const *message)
 * {
 *     FILE *fp = (FILE *) user;
 *     (void) file;
 *     (void) line;
 *     (void) function;
 *     fprintf(fp, "[%s] %s\n", osh_diag_level_name(level), message);
 * }
 *
 * struct osh_diag_sink diag = {
 *     .emit = cli_diag,
 *     .user = stderr,
 *     .min_level = OSH_DIAG_LEVEL_INFO
 * };
 * ```
 *
 * The sink is borrowed by library objects such as @ref osh_simulation.
 *
 * @param[in] diag      Borrowed diagnostics sink; NULL means silent.
 * @param[in] level     One of @ref osh_diag_level.
 * @param[in] file      Source file of the diagnostic, or NULL.
 * @param[in] line      Source line of the diagnostic, or 0.
 * @param[in] function  Source function of the diagnostic, or NULL.
 * @param[in] fmt       `printf`-style format string.
 */
void osh_diag_emitf(
    struct osh_diag_sink const *diag, int level, char const *file, int line, char const *function, char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 6, 7)))
#endif
    ;

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_DIAG_H */
