#include "osh_logger.h"

#include <stdio.h>
#include <stdlib.h>

char const *osh_diag_level_name(int level) {
    switch (level) {
    case OSH_DIAG_LEVEL_TRACE:
        return "TRACE";
    case OSH_DIAG_LEVEL_DEBUG:
        return "DEBUG";
    case OSH_DIAG_LEVEL_INFO:
        return "INFO";
    case OSH_DIAG_LEVEL_WARN:
        return "WARN";
    case OSH_DIAG_LEVEL_ERROR:
        return "ERROR";
    case OSH_DIAG_LEVEL_FATAL:
        return "FATAL";
    case OSH_DIAG_LEVEL_OFF:
        return "OFF";
    default:
        return "UNKNOWN";
    }
}

void osh_diag_emitfv(struct osh_diag_sink const *diag,
                     int level,
                     char const *file,
                     int line,
                     char const *function,
                     char const *fmt,
                     va_list ap) {
    char stackbuf[512];
    char *heapbuf = NULL;
    char const *msg = NULL;
    va_list ap_copy;
    int needed;

    if (!diag || !diag->emit) {
        return;
    }
    if (diag->min_level >= OSH_DIAG_LEVEL_OFF || level < diag->min_level) {
        return;
    }

    va_copy(ap_copy, ap);
    needed = vsnprintf(stackbuf, sizeof(stackbuf), fmt, ap_copy);
    va_end(ap_copy);

    if (needed < 0) {
        msg = "diagnostics formatting failed";
    } else if ((size_t) needed < sizeof(stackbuf)) {
        msg = stackbuf;
    } else {
        heapbuf = (char *) malloc((size_t) needed + 1u);
        if (!heapbuf) {
            msg = "diagnostics allocation failed";
        } else {
            va_copy(ap_copy, ap);
            (void) vsnprintf(heapbuf, (size_t) needed + 1u, fmt, ap_copy);
            va_end(ap_copy);
            msg = heapbuf;
        }
    }

    diag->emit(diag->user, level, file, line, function, msg ? msg : "");
    free(heapbuf);
}

void osh_diag_emitf(struct osh_diag_sink const *diag,
                    int level,
                    char const *file,
                    int line,
                    char const *function,
                    char const *fmt,
                    ...) {
    va_list ap;

    va_start(ap, fmt);
    osh_diag_emitfv(diag, level, file, line, function, fmt, ap);
    va_end(ap);
}

void osh_alloc_failed(char const *fmt, ...) {
    va_list ap;

    fputs("[FATAL] ", stderr);
    va_start(ap, fmt);
    if (fmt && *fmt) {
        (void) vfprintf(stderr, fmt, ap);
    } else {
        fputs("memory allocation failed", stderr);
    }
    va_end(ap);
    fputc('\n', stderr);
    fflush(stderr);
    abort();
}
