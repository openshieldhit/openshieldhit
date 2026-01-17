#include "osh_rc.h"

const char *osh_strerr(int c) {
    switch (c) {
    case OSH_OK:
        return "OK";
    case OSH_EINVAL:
        return "Invalid argument";
    case OSH_ENOMEM:
        return "Out of memory";
    case OSH_EIO:
        return "I/O error";
    case OSH_EPARSE:
        return "Parse error";
    case OSH_EINCOMPLETE:
        return "Incomplete configuration";
    case OSH_ENOTSUP:
        return "Not supported";
    case OSH_ESTATE:
        return "Invalid internal state";
    default:
        return "Unknown error";
    }
}
