#ifndef OSH_GEMCA_BODY_PREPARE_H
#define OSH_GEMCA_BODY_PREPARE_H

#include "gemca/osh_gemca2.h"

/**
 * @file osh_gemca_body_prepare.h
 * @brief Internal GEMCA body-allocation helpers used by the prepare layer.
 *
 * @details
 * These helpers operate on the internal compatibility workspace
 * @ref gemca_workspace built from the public cold
 * @ref osh_geometry_workspace. They are not file parsers.
 */

enum osh_status osh_gemca_body_init(struct body **body);

#endif /* OSH_GEMCA_BODY_PREPARE_H */
