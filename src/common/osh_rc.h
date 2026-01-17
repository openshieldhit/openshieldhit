#ifndef _OSH_RC_H
#define _OSH_RC_H

/* Library-level status return codes */
enum osh_status {
    OSH_OK = 0,
    OSH_EINVAL,      /* bad argument/config */
    OSH_ENOMEM,      /* out of memory */
    OSH_EIO,         /* I/O error */
    OSH_EPARSE,      /* parse error */
    OSH_EINCOMPLETE, /* incomplete data */
    OSH_ENOTSUP,     /* operation not supported */
    OSH_ESTATE       /* inconsistent state */
};

const char *osh_strerr(int code);

#endif /* !_OSH_RC_H */
