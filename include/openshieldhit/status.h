#ifndef OPENSHIELDHIT_STATUS_H
#define OPENSHIELDHIT_STATUS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Component-level status codes shared by core modules and apps.
 *
 * @details
 * This is the low-level status enum used by modular OpenShieldHIT components
 * such as beam, material, geometry, transport, and application/parser layers.
 * It is distinct from the higher-level
 * @ref openshieldhit_status enum exposed by `openshieldhit/openshieldhit.h`,
 * which belongs to the opaque application-style context API.
 */
enum osh_status { OSH_OK = 0, OSH_EINVAL, OSH_ENOMEM, OSH_EIO, OSH_EPARSE, OSH_EINCOMPLETE, OSH_ENOTSUP, OSH_ESTATE };

char const *osh_strerr(int code);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_STATUS_H */
