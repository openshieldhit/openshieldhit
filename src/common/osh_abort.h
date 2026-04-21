#ifndef OSH_ABORT_H
#define OSH_ABORT_H

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Write an OOM message to stderr and abort the process.
 *
 * @details
 * For use only by low-level library utilities that cannot propagate an error
 * return because they have no owning context in scope. External callers should
 * use abort() directly or handle allocation failures through their own code.
 */
void osh_abort_oomf(char const *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2))) __attribute__((noreturn))
#endif
    ;

#ifdef __cplusplus
}
#endif

#endif /* OSH_ABORT_H */
