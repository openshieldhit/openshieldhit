/* Public diagnostics and legacy logger API for OpenShieldHIT.
 *
 * The preferred runtime-facing abstraction is the caller-owned diagnostics
 * sink (`struct osh_diag_sink`).  It lets embedding applications decide how
 * messages are handled without forcing the library to own global logging
 * policy.
 *
 * The older process-global logger API remains available during the migration
 * and is still used by untouched modules outside the current runtime path.
 * New runtime-facing code should prefer `osh_diag_sink`.
 */

#ifndef OPENSHIELDHIT_LOGGER_H
#define OPENSHIELDHIT_LOGGER_H

#include <stdarg.h>
#include <stddef.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Severity levels for caller-owned diagnostics sinks.
 *
 * @details
 * This is the minimal diagnostics surface for the runtime path.  Callers may
 * route emitted messages to stderr, a GUI, a file, a test buffer, or ignore
 * them entirely.
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
 * This helper is intended primarily for library-internal runtime paths that
 * already own a borrowed @ref osh_diag_sink.  It checks for NULL sinks and
 * performs level filtering before invoking @ref osh_diag_emit_fn.
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

/* -------------------------------------------------------------------------
 * Legacy global logger API
 *
 * Kept during migration for modules that have not yet moved to the explicit
 * diagnostics-sink path.  New runtime-facing code should prefer
 * `struct osh_diag_sink`.
 * ------------------------------------------------------------------------- */

enum osh_log_level {
    OSH_LOG_TRACE = 0,
    OSH_LOG_DEBUG = 1,
    OSH_LOG_INFO = 2,
    OSH_LOG_WARN = 3,
    OSH_LOG_ERROR = 4,
    OSH_LOG_FATAL = 5,
    OSH_LOG_OFF = 6
};

char const *osh_log_level_name(int level);

enum osh_log_flag {
    OSH_LOG_F_NONE = 0u,
    OSH_LOG_F_TIMESTAMP = 1u << 0,
    OSH_LOG_F_THREAD_ID = 1u << 1,
    OSH_LOG_F_FILELINE = 1u << 2,
    OSH_LOG_F_FUNCTION = 1u << 3
};

struct osh_logger;

enum osh_status osh_log_init(int level, unsigned flags);
void osh_log_close(void);
enum osh_status osh_log_add_file(char const *path, int append);
enum osh_status osh_log_enable_stdout(int enable);
enum osh_status osh_log_set_level(int level);
int osh_log_get_level(void);
void osh_log_set_flags(unsigned flags);
unsigned osh_log_get_flags(void);
void osh_log_flush(void);

struct osh_logger *osh_logger_create(int level, unsigned flags);
void osh_logger_destroy(struct osh_logger *lg);
enum osh_status osh_logger_add_file(struct osh_logger *lg, char const *path, int append);
enum osh_status osh_logger_enable_stdout(struct osh_logger *lg, int enable);
void osh_logger_set_level(struct osh_logger *lg, int level);
int osh_logger_get_level(struct osh_logger const *lg);
void osh_logger_set_flags(struct osh_logger *lg, unsigned flags);
unsigned osh_logger_get_flags(struct osh_logger const *lg);
void osh_logger_flush(struct osh_logger *lg);

typedef void (*osh_log_write_cb)(void *user, char const *msg, size_t len);
enum osh_status osh_logger_set_callback(struct osh_logger *lg, osh_log_write_cb cb, void *user);

void osh_logger_logv_ex(struct osh_logger *lg,
                        int level,
                        unsigned flags_override,
                        char const *file,
                        int line,
                        char const *function,
                        char const *fmt,
                        va_list ap);

void osh_logger_log_ex(struct osh_logger *lg,
                       int level,
                       unsigned flags_override,
                       char const *file,
                       int line,
                       char const *function,
                       char const *fmt,
                       ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 7, 8)))
#endif
    ;

void osh_logger_logv(struct osh_logger *lg, int level, char const *fmt, va_list ap);
void osh_logger_log(struct osh_logger *lg, int level, char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 3, 4)))
#endif
    ;

void osh_trace(char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_debug(char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_info(char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_warn(char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_error(char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_alloc_failed(char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2))) __attribute__((noreturn))
#endif
    ;

struct osh_logger *osh_log_default(void);

#ifndef OSH_LOG_USE_FILELINE
#define OSH_LOG_USE_FILELINE 1
#endif

#if OSH_LOG_USE_FILELINE
#define OSH_LOG(lvl, fmt, ...)                                                                                         \
    osh_logger_log_ex(osh_log_default(), (lvl), 0u, __FILE__, __LINE__, __func__, (fmt), ##__VA_ARGS__)
#else
#define OSH_LOG(lvl, fmt, ...) osh_logger_log(osh_log_default(), (lvl), (fmt), ##__VA_ARGS__)
#endif

#define OSH_TRACE(fmt, ...) OSH_LOG(OSH_LOG_TRACE, (fmt), ##__VA_ARGS__)
#define OSH_DEBUG(fmt, ...) OSH_LOG(OSH_LOG_DEBUG, (fmt), ##__VA_ARGS__)
#define OSH_INFO(fmt, ...) OSH_LOG(OSH_LOG_INFO, (fmt), ##__VA_ARGS__)
#define OSH_WARN(fmt, ...) OSH_LOG(OSH_LOG_WARN, (fmt), ##__VA_ARGS__)
#define OSH_ERROR(fmt, ...) OSH_LOG(OSH_LOG_ERROR, (fmt), ##__VA_ARGS__)

#define OSH_LOG_HLINE "------------------------------------------------------------\n"
#define OSH_LOG_INDENT "    "

static char const *const osh_log_offon[2] = {"OFF", "ON"};

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_LOGGER_H */
