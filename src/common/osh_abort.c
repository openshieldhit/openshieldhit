#include "osh_abort.h"

#include <stdio.h>
#include <stdlib.h>

void osh_abort_oomf(char const *fmt, ...) {
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
