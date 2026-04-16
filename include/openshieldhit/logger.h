/* Public logger API for OpenShieldHIT.
 *
 * This is the stable logging surface used by core, frontends, and embedding
 * applications.  Implementation details remain in src/common/.
 */

#ifndef OPENSHIELDHIT_LOGGER_H
#define OPENSHIELDHIT_LOGGER_H

#include <stdarg.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

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

int osh_log_init(int level, unsigned flags);
void osh_log_close(void);
int osh_log_add_file(char const *path, int append);
int osh_log_enable_stdout(int enable);
int osh_log_set_level(int level);
int osh_log_get_level(void);
void osh_log_set_flags(unsigned flags);
unsigned osh_log_get_flags(void);
void osh_log_flush(void);

struct osh_logger *osh_logger_create(int level, unsigned flags);
void osh_logger_destroy(struct osh_logger *lg);
int osh_logger_add_file(struct osh_logger *lg, char const *path, int append);
int osh_logger_enable_stdout(struct osh_logger *lg, int enable);
void osh_logger_set_level(struct osh_logger *lg, int level);
int osh_logger_get_level(struct osh_logger const *lg);
void osh_logger_set_flags(struct osh_logger *lg, unsigned flags);
unsigned osh_logger_get_flags(struct osh_logger const *lg);
void osh_logger_flush(struct osh_logger *lg);

typedef void (*osh_log_write_cb)(void *user, char const *msg, size_t len);
int osh_logger_set_callback(struct osh_logger *lg, osh_log_write_cb cb, void *user);

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
