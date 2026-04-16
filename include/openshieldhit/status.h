#ifndef OPENSHIELDHIT_STATUS_H
#define OPENSHIELDHIT_STATUS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Component-level status codes shared by core modules and apps.
 *
 * @details
 * This is the status enum used throughout osh_core library components —
 * beam, geometry, material, transport, scoring — and application/parser layers.
 * Application code maps these codes to process exit values or user-facing
 * error messages as appropriate.
 */
enum osh_status { OSH_OK = 0, OSH_EINVAL, OSH_ENOMEM, OSH_EIO, OSH_EPARSE, OSH_EINCOMPLETE, OSH_ENOTSUP, OSH_ESTATE };

char const *osh_strerr(int code);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_STATUS_H */
