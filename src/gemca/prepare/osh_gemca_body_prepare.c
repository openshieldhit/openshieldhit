#include "gemca/prepare/osh_gemca_body_prepare.h"

#include <stdlib.h>

/**
 * @brief Allocate and zero-initialise one internal GEMCA body object.
 *
 * @details
 * The public cold geometry API stores bodies as flat value objects in
 * @ref osh_geometry_workspace. During prepare those are copied into the
 * internal compatibility workspace @ref osh_gemca_prepared, which still uses
 * pointer-owned @ref body objects. This helper allocates one such internal
 * object.
 *
 * @param[out] body Receives the newly allocated internal body.
 *
 * @returns OSH_OK on success, OSH_EINVAL when @p body is NULL,
 *          OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_gemca_body_init(struct body **body) {
    if (body == NULL) {
        return OSH_EINVAL;
    }

    *body = NULL;
    *body = calloc(1, sizeof(struct body));
    if (*body == NULL) {
        return OSH_ENOMEM;
    }
    return OSH_OK;
}
