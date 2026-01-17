#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "osh_logger.h"

#if defined(_WIN32)
#include <windows.h>
#else
#include <pthread.h>
#endif

struct osh_mutex {
#if defined(_WIN32)
    CRITICAL_SECTION cs;
#else
    pthread_mutex_t m;
#endif
};

/* -----------------------------
 * forward declarations
 * ----------------------------- */

/* mutex helpers (file-local) */
static void _mutex_init(struct osh_mutex *);
static void _mutex_lock(struct osh_mutex *);
static void _mutex_unlock(struct osh_mutex *);
static void _mutex_destroy(struct osh_mutex *);

/* logger lifecycle (public; declared in .h, but fine to list here too) */
struct osh_logger *osh_logger_create(int level, unsigned flags);
void osh_logger_destroy(struct osh_logger *lg);

/* logger config (public; declared in .h, but fine to list here too) */
void osh_logger_set_level(struct osh_logger *lg, int level);
int osh_logger_get_level(struct osh_logger const *lg);
void osh_logger_set_flags(struct osh_logger *lg, unsigned flags);
unsigned int osh_logger_get_flags(struct osh_logger const *lg);
void osh_logger_flush(struct osh_logger *lg);

/* default logger (public; declared in .h, but fine to list here too) */
struct osh_logger *osh_log_default(void);
int osh_log_init(int level, unsigned flags);
void osh_log_close(void);
int osh_log_set_level(int level);
int osh_log_get_level(void);
void osh_log_set_flags(unsigned flags);
unsigned osh_log_get_flags(void);
void osh_log_flush(void);

/* sinks / callbacks (public) */
int osh_log_add_file(const char *path, int append);
int osh_log_enable_stdout(int enable);

int osh_logger_add_file(struct osh_logger *lg, const char *path, int append);
int osh_logger_enable_stdout(struct osh_logger *lg, int enable);
int osh_logger_set_callback(struct osh_logger *lg, osh_log_write_cb cb, void *user);

/* wrappers / convenience (public) */
const char *osh_log_level_name(int level);

void osh_logger_log_ex(struct osh_logger *lg, int level, unsigned flags_override, const char *file, int line,
                       const char *function, const char *fmt, ...);
void osh_logger_logv(struct osh_logger *lg, int level, const char *fmt, va_list ap);
void osh_logger_log(struct osh_logger *lg, int level, const char *fmt, ...);

void osh_trace(const char *fmt, ...);
void osh_debug(const char *fmt, ...);
void osh_info(const char *fmt, ...);
void osh_warn(const char *fmt, ...);
void osh_error(int exit_code, const char *fmt, ...);
void osh_alloc_failed(size_t size);

/* formatting / output helpers (file-local) */
static unsigned _effective_flags(struct osh_logger *lg, unsigned flags_override);
static const char *_level_name(int level);
static int _is_errorish(int level);

static size_t _vsnprintf0(char *dst, size_t cap, const char *fmt, va_list ap);
static size_t _snprintf0(char *dst, size_t cap, const char *fmt, ...);

static size_t _append(char *dst, size_t cap, size_t off, const char *s);
static size_t _append_n(char *dst, size_t cap, size_t off, const char *s, size_t n);
static size_t _append_char(char *dst, size_t cap, size_t off, char c);
static size_t _append_int(char *dst, size_t cap, size_t off, int v);
static size_t _append_timestamp(char *dst, size_t cap, size_t off);

static size_t _format_prefix(char *dst, size_t cap, int level, unsigned flags, const char *file, int line,
                             const char *function);

static FILE *_open_log_file(const char *path, int append);

static void _write_record_unlocked(struct osh_logger *lg, int level, const char *prefix, size_t prefix_len,
                                   const char *fmt, va_list ap);

/* -----------------------------
 * types / state
 * ----------------------------- */

struct osh_logger {
    /* configuration */
    int level;
    unsigned int flags;

    /* outputs */
    FILE *fp_file;  /* optional logfile */
    int use_stdout; /* 0/1 */

    /* callback sink */
    void (*cb)(void *user, const char *msg, size_t len);
    void *cb_user;

    /* synchronization */
    struct osh_mutex lock;

    /* lifecycle */
    int closed;
};

static struct osh_logger *g_default_logger = NULL;

/* -----------------------------
 * mutex helpers
 * ----------------------------- */

static void _mutex_init(struct osh_mutex *m) {
#if defined(_WIN32)
    InitializeCriticalSection(&m->cs);
#else
    /* default mutex is fine; we don't need recursive mutex if we avoid re-entry */
    (void)pthread_mutex_init(&m->m, NULL);
#endif
}

static void _mutex_lock(struct osh_mutex *m) {
#if defined(_WIN32)
    EnterCriticalSection(&m->cs);
#else
    (void)pthread_mutex_lock(&m->m);
#endif
}

static void _mutex_unlock(struct osh_mutex *m) {
#if defined(_WIN32)
    LeaveCriticalSection(&m->cs);
#else
    (void)pthread_mutex_unlock(&m->m);
#endif
}

static void _mutex_destroy(struct osh_mutex *m) {
#if defined(_WIN32)
    DeleteCriticalSection(&m->cs);
#else
    (void)pthread_mutex_destroy(&m->m);
#endif
}

/* -----------------------------
 * explicit logger instances
 * ----------------------------- */

struct osh_logger *osh_logger_create(int level, unsigned flags) {
    struct osh_logger *lg = (struct osh_logger *)calloc(1, sizeof(*lg));
    if (!lg)
        return NULL;

    lg->level = level;
    lg->flags = flags;
    lg->fp_file = NULL;
    lg->use_stdout = 0;
    lg->cb = NULL;
    lg->cb_user = NULL;
    lg->closed = 0;

    _mutex_init(&lg->lock);
    return lg;
}

void osh_logger_destroy(struct osh_logger *lg) {
    if (!lg)
        return;

    _mutex_lock(&lg->lock);

    if (!lg->closed) {
        if (lg->fp_file) {
            fclose(lg->fp_file);
            lg->fp_file = NULL;
        }
        lg->closed = 1;
    }

    _mutex_unlock(&lg->lock);
    _mutex_destroy(&lg->lock);
    free(lg);
}

void osh_logger_set_level(struct osh_logger *lg, int level) {
    if (!lg)
        return;

    _mutex_lock(&lg->lock);
    lg->level = level;
    _mutex_unlock(&lg->lock);
}

int osh_logger_get_level(struct osh_logger const *lg) {
    int level = OSH_LOG_OFF;
    if (!lg)
        return level;

    /* Cast away const only to use the same mutex; or store level atomically later. */
    _mutex_lock((struct osh_mutex *)&lg->lock);
    level = lg->level;
    _mutex_unlock((struct osh_mutex *)&lg->lock);
    return level;
}

void osh_logger_set_flags(struct osh_logger *lg, unsigned flags) {
    if (!lg)
        return;

    _mutex_lock(&lg->lock);
    lg->flags = flags;
    _mutex_unlock(&lg->lock);
}

unsigned int osh_logger_get_flags(struct osh_logger const *lg) {
    unsigned flags = 0u;
    if (!lg)
        return flags;

    _mutex_lock((struct osh_mutex *)&lg->lock);
    flags = lg->flags;
    _mutex_unlock((struct osh_mutex *)&lg->lock);
    return flags;
}

void osh_logger_flush(struct osh_logger *lg) {
    if (!lg)
        return;

    _mutex_lock(&lg->lock);
    if (lg->fp_file)
        fflush(lg->fp_file);
    fflush(stderr);
    if (lg->use_stdout)
        fflush(stdout);
    _mutex_unlock(&lg->lock);
}

/* -----------------------------
 * default global logger
 * ----------------------------- */

struct osh_logger *osh_log_default(void) {
    return g_default_logger;
}

int osh_log_init(int level, unsigned flags) {
    if (g_default_logger) {
        /* already initialized: just update config */
        osh_logger_set_level(g_default_logger, level);
        osh_logger_set_flags(g_default_logger, flags);
        return 0;
    }

    g_default_logger = osh_logger_create(level, flags);
    return g_default_logger ? 0 : -1;
}

void osh_log_close(void) {
    if (!g_default_logger)
        return;

    osh_logger_destroy(g_default_logger);
    g_default_logger = NULL;
}

int osh_log_set_level(int level) {
    if (!g_default_logger)
        return -1;

    osh_logger_set_level(g_default_logger, level);
    return 0;
}

int osh_log_get_level(void) {
    if (!g_default_logger)
        return OSH_LOG_OFF;

    return osh_logger_get_level(g_default_logger);
}

void osh_log_set_flags(unsigned flags) {
    if (!g_default_logger)
        return;

    osh_logger_set_flags(g_default_logger, flags);
}

unsigned osh_log_get_flags(void) {
    if (!g_default_logger)
        return 0u;

    return osh_logger_get_flags(g_default_logger);
}

void osh_log_flush(void) {
    if (!g_default_logger)
        return;

    osh_logger_flush(g_default_logger);
}

/* ------------------------------------------------------------
 * helpers
 * ------------------------------------------------------------ */

static unsigned _effective_flags(struct osh_logger *lg, unsigned flags_override) {
    return flags_override ? flags_override : lg->flags;
}

static const char *_level_name(int level) {
    switch (level) {
    case OSH_LOG_TRACE:
        return "TRACE";
    case OSH_LOG_DEBUG:
        return "DEBUG";
    case OSH_LOG_INFO:
        return "INFO";
    case OSH_LOG_WARN:
        return "WARN";
    case OSH_LOG_ERROR:
        return "ERROR";
    case OSH_LOG_FATAL:
        return "FATAL";
    case OSH_LOG_OFF:
        return "OFF";
    default:
        return "LOG";
    }
}

static int _is_errorish(int level) {
    return (level >= OSH_LOG_WARN);
}

static size_t _vsnprintf0(char *dst, size_t cap, const char *fmt, va_list ap) {
    int n = vsnprintf(dst, cap, fmt, ap);
    if (n < 0) {
        if (cap)
            dst[0] = '\0';
        return 0u;
    }
    return (size_t)n; /* would-have length (excluding NUL) */
}

static size_t _snprintf0(char *dst, size_t cap, const char *fmt, ...) {
    va_list ap;
    size_t n;

    va_start(ap, fmt);
    n = _vsnprintf0(dst, cap, fmt, ap);
    va_end(ap);

    return n;
}

static size_t _append_n(char *dst, size_t cap, size_t off, const char *s, size_t n) {
    if (!dst || cap == 0)
        return off;
    if (off >= cap)
        return off;

    size_t space = cap - off - 1; /* keep room for NUL */
    size_t tocpy = (n < space) ? n : space;

    if (tocpy)
        memcpy(dst + off, s, tocpy);
    off += tocpy;
    dst[off] = '\0';
    return off;
}

static size_t _append(char *dst, size_t cap, size_t off, const char *s) {
    if (!s)
        return off;
    return _append_n(dst, cap, off, s, strlen(s));
}

static size_t _append_char(char *dst, size_t cap, size_t off, char c) {
    if (!dst || cap == 0)
        return off;
    if (off + 1 >= cap)
        return off;
    dst[off++] = c;
    dst[off] = '\0';
    return off;
}

static size_t _append_int(char *dst, size_t cap, size_t off, int v) {
    char tmp[32];
    (void)_snprintf0(tmp, sizeof(tmp), "%d", v);
    return _append(dst, cap, off, tmp);
}

static size_t _append_timestamp(char *dst, size_t cap, size_t off) {
    /* Format: YYYY-MM-DD HH:MM:SS.mmm */
    time_t t = time(NULL);

#if defined(_WIN32)
    struct tm tmv;
    /* localtime_s returns 0 on success */
    if (localtime_s(&tmv, &t) != 0) {
        return off;
    }
#else
    struct tm tmv;
    if (localtime_r(&t, &tmv) == NULL) {
        return off;
    }
#endif

    char buf[64];
    size_t n = strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tmv);
    if (n == 0) {
        return off;
    }

    /* milliseconds */
    unsigned ms = 0;

#if defined(_WIN32)
    {
        /* GetSystemTimeAsFileTime: 100-ns intervals since Jan 1, 1601 */
        FILETIME ft;
        ULARGE_INTEGER ui;

        GetSystemTimeAsFileTime(&ft);
        ui.LowPart = ft.dwLowDateTime;
        ui.HighPart = ft.dwHighDateTime;

        /* 100-ns -> ms */
        ms = (unsigned)((ui.QuadPart / 10000ULL) % 1000ULL);
    }
#else
    {
#if defined(CLOCK_REALTIME)
        struct timespec ts;
        if (clock_gettime(CLOCK_REALTIME, &ts) == 0) {
            ms = (unsigned)(ts.tv_nsec / 1000000L);
        }
#endif
    }
#endif

    off = _append_n(dst, cap, off, buf, n);
    off = _append_char(dst, cap, off, '.');

    /* print ms as 3 digits */
    {
        char msbuf[8];
        (void)_snprintf0(msbuf, sizeof(msbuf), "%03u", ms);
        off = _append(dst, cap, off, msbuf);
    }

    off = _append_char(dst, cap, off, ' ');
    return off;
}

static size_t _format_prefix(char *dst, size_t cap, int level, unsigned flags, const char *file, int line,
                             const char *function) {
    size_t off = 0;

    if (flags & OSH_LOG_F_TIMESTAMP) {
        off = _append_timestamp(dst, cap, off);
    }

    off = _append_char(dst, cap, off, '[');
    off = _append(dst, cap, off, _level_name(level));
    off = _append(dst, cap, off, "] ");

    if ((flags & OSH_LOG_F_FILELINE) && file && *file) {
        off = _append(dst, cap, off, file);
        off = _append_char(dst, cap, off, ':');
        off = _append_int(dst, cap, off, line);
        off = _append(dst, cap, off, " ");
    }

    if ((flags & OSH_LOG_F_FUNCTION) && function && *function) {
        off = _append_char(dst, cap, off, '(');
        off = _append(dst, cap, off, function);
        off = _append(dst, cap, off, "): ");
    }

    return off;
}

static void _write_record_unlocked(struct osh_logger *lg, int level, const char *prefix, size_t prefix_len,
                                   const char *fmt, va_list ap) {
    FILE *fp_primary = NULL;

    if (_is_errorish(level)) {
        fp_primary = stderr;
    } else {
        fp_primary = lg->use_stdout ? stdout : stderr;
    }

    /* Make copies BEFORE any consumption */
    va_list ap_primary;
    va_list ap_file;

    va_copy(ap_primary, ap);
    va_copy(ap_file, ap);

    if (fp_primary) {
        if (prefix && prefix_len)
            (void)fwrite(prefix, 1, prefix_len, fp_primary);

        (void)vfprintf(fp_primary, fmt, ap_primary);
        (void)fputc('\n', fp_primary);
    }

    if (lg->fp_file) {
        if (prefix && prefix_len)
            (void)fwrite(prefix, 1, prefix_len, lg->fp_file);

        (void)vfprintf(lg->fp_file, fmt, ap_file);
        (void)fputc('\n', lg->fp_file);
    }

    va_end(ap_primary);
    va_end(ap_file);

    /* Callback — optional, bounded */
    if (lg->cb) {
        if (prefix && prefix_len)
            lg->cb(lg->cb_user, prefix, prefix_len);

        char tmp[512];
        va_list ap_cb;
        va_copy(ap_cb, ap);
        int n = vsnprintf(tmp, sizeof(tmp), fmt, ap_cb);
        va_end(ap_cb);

        if (n > 0) {
            lg->cb(lg->cb_user, tmp, (size_t)((n < (int)sizeof(tmp)) ? n : (int)sizeof(tmp) - 1));
        }

        const char nl = '\n';
        lg->cb(lg->cb_user, &nl, 1);
    }
}

/* ------------------------------------------------------------
 * the big one: now short and readable
 * ------------------------------------------------------------ */

void osh_logger_logv_ex(struct osh_logger *lg, int level, unsigned flags_override, const char *file, int line,
                        const char *function, const char *fmt, va_list ap) {
    char prefix[512]; /* Prefix is small; safe on stack */
    size_t prefix_len;

    if (!lg || !fmt)
        return;

    /* Fast reject (best-effort; races with setters are acceptable) */
    if (lg->closed)
        return;
    if (level >= OSH_LOG_OFF)
        return;
    if (level < lg->level)
        return;

    unsigned flags = _effective_flags(lg, flags_override);

    prefix_len = _format_prefix(prefix, sizeof(prefix), level, flags, file, line, function);

    _mutex_lock(&lg->lock);
    if (!lg->closed) {
        _write_record_unlocked(lg, level, prefix, prefix_len, fmt, ap);
    }
    _mutex_unlock(&lg->lock);
}

/* -----------------------------
 * exported helpers
 * ----------------------------- */

const char *osh_log_level_name(int level) {
    /* Keep consistent with internal _level_name() */
    return _level_name(level);
}

/* -----------------------------
 * sinks / callbacks
 * ----------------------------- */

static FILE *_open_log_file(const char *path, int append) {
    if (!path || !*path)
        return NULL;

    /* On Windows, text mode is fine; if you ever need strict \n behavior, use "ab"/"wb". */
    return fopen(path, append ? "a" : "w");
}

int osh_logger_add_file(struct osh_logger *lg, const char *path, int append) {
    if (!lg)
        return -1;

    FILE *fp = _open_log_file(path, append);
    if (!fp)
        return -1;

    _mutex_lock(&lg->lock);
    if (lg->fp_file) {
        fclose(lg->fp_file);
    }
    lg->fp_file = fp;
    _mutex_unlock(&lg->lock);

    return 0;
}

int osh_log_add_file(const char *path, int append) {
    return osh_logger_add_file(osh_log_default(), path, append);
}

int osh_logger_enable_stdout(struct osh_logger *lg, int enable) {
    if (!lg)
        return -1;

    _mutex_lock(&lg->lock);
    lg->use_stdout = enable ? 1 : 0;
    _mutex_unlock(&lg->lock);

    return 0;
}

int osh_log_enable_stdout(int enable) {
    return osh_logger_enable_stdout(osh_log_default(), enable);
}

int osh_logger_set_callback(struct osh_logger *lg, osh_log_write_cb cb, void *user) {
    if (!lg)
        return -1;

    _mutex_lock(&lg->lock);
    lg->cb = cb;
    lg->cb_user = user;
    _mutex_unlock(&lg->lock);

    return 0;
}

/* -----------------------------
 * wrappers around logv_ex
 * ----------------------------- */

void osh_logger_log_ex(struct osh_logger *lg, int level, unsigned flags_override, const char *file, int line,
                       const char *function, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    osh_logger_logv_ex(lg, level, flags_override, file, line, function, fmt, ap);
    va_end(ap);
}

void osh_logger_logv(struct osh_logger *lg, int level, const char *fmt, va_list ap) {
    osh_logger_logv_ex(lg, level, 0u, NULL, 0, NULL, fmt, ap);
}

void osh_logger_log(struct osh_logger *lg, int level, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    osh_logger_logv_ex(lg, level, 0u, NULL, 0, NULL, fmt, ap);
    va_end(ap);
}

/* -----------------------------
 * default logger convenience functions
 * ----------------------------- */

void osh_trace(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    osh_logger_logv_ex(osh_log_default(), OSH_LOG_TRACE, 0u, __FILE__, __LINE__, __func__, fmt, ap);
    va_end(ap);
}

void osh_debug(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    osh_logger_logv_ex(osh_log_default(), OSH_LOG_DEBUG, 0u, __FILE__, __LINE__, __func__, fmt, ap);
    va_end(ap);
}

void osh_info(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    osh_logger_logv_ex(osh_log_default(), OSH_LOG_INFO, 0u, __FILE__, __LINE__, __func__, fmt, ap);
    va_end(ap);
}

void osh_warn(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    osh_logger_logv_ex(osh_log_default(), OSH_LOG_WARN, 0u, __FILE__, __LINE__, __func__, fmt, ap);
    va_end(ap);
}
void osh_error(int exit_code, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    osh_logger_logv_ex(osh_log_default(), OSH_LOG_FATAL, 0u, __FILE__, __LINE__, __func__, fmt, ap);
    va_end(ap);

    osh_log_flush();
    exit(exit_code);
}

void osh_alloc_failed(size_t size) {
    int err = errno;

    osh_error(EXIT_FAILURE, "memory allocation failed (size=%zu): %s", size, err ? strerror(err) : "out of memory");
}