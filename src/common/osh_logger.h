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

#include "osh_exit.h"

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------
 * Levels
 * ----------------------------- */

/**
 * @enum osh_log_level
 *
 * @brief Logging levels for controlling verbosity.
 */
enum osh_log_level {
    OSH_LOG_TRACE = 0,
    OSH_LOG_DEBUG = 1,
    OSH_LOG_INFO = 2,
    OSH_LOG_WARN = 3,
    OSH_LOG_ERROR = 4,
    OSH_LOG_FATAL = 5,
    OSH_LOG_OFF = 6
};

/**
 * @brief Get the name of a logging level.
 *
 * @param[in] level The logging level.
 *
 * @return The name of the logging level as a string.
 */
const char *osh_log_level_name(int level);

/* -----------------------------
 * Flags (bitmask)
 * ----------------------------- */

/**
 * @enum osh_log_flag
 *
 * @brief Flags for customizing log output.
 */
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

/**
 * @struct osh_logger
 *
 * @brief Opaque handle for a logger instance.
 */
struct osh_logger;

/* -----------------------------
 * Default global logger (keeps existing style)
 * ----------------------------- */

/**
 * @brief Initialize the default logger.
 *
 * @param[in] level The logging level.
 * @param[in] flags A bitmask of logging flags.
 *
 * @return 0 on success, nonzero on error.
 */
int osh_log_init(int level, unsigned flags);

/**
 * @brief Close the default logger and free all resources.
 * This function is idempotent.
 */
void osh_log_close(void);

/**
 * @brief Add a file sink to the default logger.
 *
 * @param[in] path The file path to log to.
 * @param[in] append If nonzero, append to the file; otherwise, overwrite.
 *
 * @return 0 on success, nonzero on error.
 */
int osh_log_add_file(const char *path, int append);

/**
 * @brief Enable or disable stdout logging for the default logger.
 *
 * @param[in] enable 0 to disable, 1 to enable.
 *
 * @return 0 on success, nonzero on error.
 */
int osh_log_enable_stdout(int enable);

/**
 * @brief Set the logging level for the default logger.
 *
 * @param[in] level The new logging level.
 *
 * @return 0 on success, nonzero on error.
 */
int osh_log_set_level(int level);

/**
 * @brief Get the current logging level of the default logger.
 *
 * @return The current logging level.
 */
int osh_log_get_level(void);

/**
 * @brief Set the logging flags for the default logger.
 *
 * @param[in] flags A bitmask of logging flags.
 */
void osh_log_set_flags(unsigned flags);

/**
 * @brief Get the current logging flags of the default logger.
 *
 * @return The current logging flags.
 */
unsigned osh_log_get_flags(void);

/**
 * @brief Flush all pending log messages in the default logger.
 */
void osh_log_flush(void);

/* -----------------------------
 * Explicit logger instances (optional but useful for tests / libraries)
 * ----------------------------- */

/**
 * @brief Create a new logger instance.
 *
 * @param[in] level The logging level.
 * @param[in] flags A bitmask of logging flags.
 *
 * @return A pointer to the new logger instance, or NULL on error.
 */
struct osh_logger *osh_logger_create(int level, unsigned flags);

/**
 * @brief Destroy a logger instance and free its resources.
 *
 * @param[in] lg The logger instance to destroy.
 */
void osh_logger_destroy(struct osh_logger *lg);

/**
 * @brief Add a file sink to a logger instance.
 *
 * @param[in] lg The logger instance.
 * @param[in] path The file path to log to.
 * @param[in] append If nonzero, append to the file; otherwise, overwrite.
 *
 * @return 0 on success, nonzero on error.
 */
int osh_logger_add_file(struct osh_logger *lg, const char *path, int append);

/**
 * @brief Enable or disable stdout logging for a logger instance.
 *
 * @param[in] lg The logger instance.
 * @param[in] enable 0 to disable, 1 to enable.
 *
 * @return 0 on success, nonzero on error.
 */
int osh_logger_enable_stdout(struct osh_logger *lg, int enable);

/**
 * @brief Set the logging level for a logger instance.
 *
 * @param[in] lg The logger instance.
 * @param[in] level The new logging level.
 */
void osh_logger_set_level(struct osh_logger *lg, int level);

/**
 * @brief Get the current logging level of a logger instance.
 *
 * @param[in] lg The logger instance.
 * @return The current logging level.
 */
int osh_logger_get_level(const struct osh_logger *lg);

/**
 * @brief Set the logging flags for a logger instance.
 *
 * @param[in] lg The logger instance.
 * @param[in] flags A bitmask of logging flags.
 */
void osh_logger_set_flags(struct osh_logger *lg, unsigned flags);

/**
 * @brief Get the current logging flags of a logger instance.
 *
 * @param[in] lg The logger instance.
 * @return The current logging flags.
 */
unsigned osh_logger_get_flags(const struct osh_logger *lg);

/**
 * @brief Flush all pending log messages in a logger instance.
 *
 * @param[in] lg The logger instance.
 */
void osh_logger_flush(struct osh_logger *lg);

/**
 * @typedef osh_log_write_cb
 * @brief Callback type for custom log sinks.
 *
 * @param[in] user User-defined context.
 * @param[in] msg The log message.
 * @param[in] len The length of the log message.
 */
typedef void (*osh_log_write_cb)(void *user, const char *msg, size_t len);

/**
 * @brief Set a custom callback sink for a logger instance.
 *
 * @param[in] lg The logger instance.
 * @param[in] cb The callback function.
 * @param[in] user User-defined context for the callback.
 *
 * @return 0 on success, nonzero on error.
 */
int osh_logger_set_callback(struct osh_logger *lg, osh_log_write_cb cb, void *user);

/* ------------------------------
 * Canonical logging entry points
 * ------------------------------
 */

/**
 * @brief Log a message with extended options.
 *
 * @param[in] lg The logger instance.
 * @param[in] level The logging level.
 * @param[in] flags_override If 0, use logger flags; otherwise, use this value.
 * @param[in] file The source file name.
 * @param[in] line The line number in the source file.
 * @param[in] function The function name.
 * @param[in] fmt The format string.
 * @param[in] ap The variable argument list.
 */
void osh_logger_logv_ex(struct osh_logger *lg,
                        int level,
                        unsigned flags_override,
                        const char *file,
                        int line,
                        const char *function,
                        const char *fmt,
                        va_list ap);

/**
 * @brief Log a message with extended options (variadic version).
 *
 * @param[in] lg The logger instance.
 * @param[in] level The logging level.
 * @param[in] flags_override If 0, use logger flags; otherwise, use this value.
 * @param[in] file The source file name.
 * @param[in] line The line number in the source file.
 * @param[in] function The function name.
 * @param[in] fmt The format string.
 * @param[in] ... The variable arguments.
 */
void osh_logger_log_ex(struct osh_logger *lg,
                       int level,
                       unsigned flags_override,
                       const char *file,
                       int line,
                       const char *function,
                       const char *fmt,
                       ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 7, 8)))
#endif
    ;

/**
 * @brief Log a message.
 *
 * @param[in] lg The logger instance.
 * @param[in] level The logging level.
 * @param[in] fmt The format string.
 * @param[in] ap The variable argument list.
 */
void osh_logger_logv(struct osh_logger *lg, int level, const char *fmt, va_list ap);

/**
 * @brief Log a message (variadic version).
 *
 * @param[in] lg The logger instance.
 * @param[in] level The logging level.
 * @param[in] fmt The format string.
 * @param[in] ... The variable arguments.
 */
void osh_logger_log(struct osh_logger *lg, int level, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 3, 4)))
#endif
    ;

/**
 * @brief Log a trace-level message using the default logger.
 *
 * @param[in] fmt The format string.
 * @param[in] ... The variable arguments.
 */
void osh_trace(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

/**
 * @brief Log a debug-level message using the default logger.
 *
 * @param[in] fmt The format string.
 * @param[in] ... The variable arguments.
 */
void osh_debug(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

/**
 * @brief Log an informational message using the default logger.
 *
 * @param[in] fmt The format string.
 * @param[in] ... The variable arguments.
 */
void osh_info(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

/**
 * @brief Log a warning message using the default logger.
 *
 * @param[in] fmt The format string.
 * @param[in] ... The variable arguments.
 */
void osh_warn(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

/**
 * @brief Log an error message and exit using the default logger.
 *
 * @param[in] exit_code The exit code as in osh_exit.h.
 * @param[in] fmt The format string.
 * @param[in] ... The variable arguments.
 */
void osh_error(int exit_code, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 2, 3), noreturn))
#endif
    ;

/**
 * @brief Log a fatal memory allocation failure and terminate the program.
 *
 * This function logs a formatted fatal error message and then calls exit().
 * It does not return.
 *
 * @param[in] fmt printf-style format string describing the failure.
 * @param[in] ... Optional format arguments.
 */
void osh_alloc_failed(const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2))) __attribute__((noreturn))
#endif
    ;

/* -----------------------------
 * Macro layer capturing file/line/function
 * ----------------------------- */

/**
 * @brief Get the default logger instance.
 *
 * @return A pointer to the default logger, or NULL if not initialized.
 */
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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* OSH_LOGGER_H */
