#ifndef OSH_GEMCA2_INTERNAL_H
#define OSH_GEMCA2_INTERNAL_H

#include "gemca/osh_gemca2.h"

/**
 * @file osh_gemca2_internal.h
 * @brief Internal gemca workspace functions not part of the public API.
 *
 * Include this header within src/gemca/ only.
 */

enum osh_status osh_gemca_prepared_free(struct osh_gemca_prepared *wg);

#endif /* OSH_GEMCA2_INTERNAL_H */
