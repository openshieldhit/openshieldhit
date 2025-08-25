#ifndef _OSH_RC_H
#define _OSH_RC_H

/* Library-level status codes */
#define OSH_OK             0
#define OSH_EINVAL         1   /* bad argument/config */
#define OSH_ENOMEM         2    /* out of memory */
#define OSH_EIO            3    /* I/O error */
#define OSH_EPARSE         4    /* parse error */
#define OSH_EINCOMPLETE    5    /* incomplete data */
#define OSH_ENOTSUP        6    /* operation not supported */
#define OSH_ESTATE         7   /* inconsistent state */

const char *osh_strerr(int code);

#endif /* !_OSH_RC_H */
