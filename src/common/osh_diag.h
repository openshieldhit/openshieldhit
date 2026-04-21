#ifndef OSH_DIAG_H
#define OSH_DIAG_H

#include "openshieldhit/diag.h"

#define OSH_DIAGF(diag, lvl, fmt, ...) osh_diag_emitf((diag), (lvl), __FILE__, __LINE__, __func__, (fmt), ##__VA_ARGS__)

#define OSH_DIAG_TRACEF(diag, fmt, ...) OSH_DIAGF((diag), OSH_DIAG_LEVEL_TRACE, (fmt), ##__VA_ARGS__)
#define OSH_DIAG_DEBUGF(diag, fmt, ...) OSH_DIAGF((diag), OSH_DIAG_LEVEL_DEBUG, (fmt), ##__VA_ARGS__)
#define OSH_DIAG_INFOF(diag, fmt, ...) OSH_DIAGF((diag), OSH_DIAG_LEVEL_INFO, (fmt), ##__VA_ARGS__)
#define OSH_DIAG_WARNF(diag, fmt, ...) OSH_DIAGF((diag), OSH_DIAG_LEVEL_WARN, (fmt), ##__VA_ARGS__)
#define OSH_DIAG_ERRORF(diag, fmt, ...) OSH_DIAGF((diag), OSH_DIAG_LEVEL_ERROR, (fmt), ##__VA_ARGS__)
#define OSH_DIAG_FATALF(diag, fmt, ...) OSH_DIAGF((diag), OSH_DIAG_LEVEL_FATAL, (fmt), ##__VA_ARGS__)

#endif /* OSH_DIAG_H */
