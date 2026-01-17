/* osh_logger.h - portable, thread-safe, re-entrant logger (C)
 *
 * Goals:
 *  - re-entrant (no recursion inside logger impl)
 *  - thread-safe
 *  - builds on POSIX + Windows (MinGW)
 *  - minimal surface area
 *  - no user-defined typedefs
 */

#ifndef OSH_LOGGER_H
#define OSH_LOGGER_H

#include <stdarg.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------
 * Levels
 * ----------------------------- */
enum osh_log_level {
    OSH_LOG_TRACE = 0,
    OSH_LOG_DEBUG = 1,
    OSH_LOG_INFO = 2,
    OSH_LOG_WARN = 3,
    OSH_LOG_ERROR = 4,
    OSH_LOG_FATAL = 5,
    OSH_LOG_OFF = 6
};

const char *osh_log_level_name(int level);

/* -----------------------------
 * Flags (bitmask)
 * ----------------------------- */
enum osh_log_flag {
    OSH_LOG_F_NONE = 0u,
    OSH_LOG_F_TIMESTAMP = 1u << 0,
    OSH_LOG_F_THREAD_ID = 1u << 1,
    OSH_LOG_F_FILELINE = 1u << 2,
    OSH_LOG_F_FUNCTION = 1u << 3
};

/* -----------------------------
 * Opaque logger handle
 * ----------------------------- */
struct osh_logger;

/* -----------------------------
 * Default global logger (keeps existing style)
 * ----------------------------- */

/* Create default logger:
 *  - always logs to stderr
 *  - level and flags set from args
 * Returns 0 on success, nonzero on error.
 */
int osh_log_init(int level, unsigned flags);

/* Close default logger and free all resources (idempotent). */
void osh_log_close(void);

/* Optional sinks on the default logger */
int osh_log_add_file(const char *path, int append);
int osh_log_enable_stdout(int enable); /* 0/1, for INFO+ for example; implementation decides */
int osh_log_set_level(int level);
int osh_log_get_level(void);
void osh_log_set_flags(unsigned flags);
unsigned osh_log_get_flags(void);
void osh_log_flush(void);

/* -----------------------------
 * Explicit logger instances (optional but useful for tests / libraries)
 * ----------------------------- */
struct osh_logger *osh_logger_create(int level, unsigned flags);
void osh_logger_destroy(struct osh_logger *lg);

int osh_logger_add_file(struct osh_logger *lg, const char *path, int append);
int osh_logger_enable_stdout(struct osh_logger *lg, int enable);
void osh_logger_set_level(struct osh_logger *lg, int level);
int osh_logger_get_level(const struct osh_logger *lg);
void osh_logger_set_flags(struct osh_logger *lg, unsigned flags);
unsigned osh_logger_get_flags(const struct osh_logger *lg);
void osh_logger_flush(struct osh_logger *lg);

/* Optional: user callback sink (no vtable, just one callback) */
typedef void (*osh_log_write_cb)(void *user, const char *msg, size_t len);
int osh_logger_set_callback(struct osh_logger *lg, osh_log_write_cb cb, void *user);

/* -----------------------------
 * Canonical logging entry points
 * ----------------------------- */

/* Extended form allows attaching file/line/function explicitly.
 * flags_override: if 0, uses logger flags. If nonzero, uses that value.
 */
void osh_logger_logv_ex(struct osh_logger *lg, int level, unsigned flags_override, const char *file, int line,
                        const char *function, const char *fmt, va_list ap);

void osh_logger_log_ex(struct osh_logger *lg, int level, unsigned flags_override, const char *file, int line,
                       const char *function, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 7, 8)))
#endif
    ;

/* Plain form without file/line/function */
void osh_logger_logv(struct osh_logger *lg, int level, const char *fmt, va_list ap);
void osh_logger_log(struct osh_logger *lg, int level, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 3, 4)))
#endif
    ;

/* Default logger convenience functions (useful when migrating existing code) */
void osh_trace(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_debug(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_info(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_warn(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;
void osh_error(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

/* Logs then exits(exit_code). noreturn attribute is best-effort. */
void osh_fatal(int exit_code, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 2, 3), noreturn))
#endif
    ;

/* -----------------------------
 * Macro layer capturing file/line/function
 * ----------------------------- */
#ifndef OSH_LOG_USE_FILELINE
#define OSH_LOG_USE_FILELINE 1
#endif

/* Accessor for default logger (may return NULL if not initialized) */
struct osh_logger *osh_log_default(void);

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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* OSH_LOGGER_H */
