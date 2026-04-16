#include "gemca/prepare/osh_gemca_body_prepare.h"

#include <stdlib.h>

#include "common/osh_logger.h"

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
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_gemca_body_init(struct body **body) {

    *body = calloc(1, sizeof(struct body));
    if (*body == NULL) {
        osh_error("osh_gemca_body_init() cannot allocate memory");
        return OSH_ENOMEM;
    }
    return OSH_OK;
}
