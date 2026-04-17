#ifndef OSH_FRONTEND_OPENSHIELDHIT_GEOMETRY_PARSE_H
#define OSH_FRONTEND_OPENSHIELDHIT_GEOMETRY_PARSE_H

#include "openshieldhit/file.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/status.h"

/**
 * @brief Parse a geo.dat file into a cold geometry workspace.
 *
 * @details
 * Reads from @p oshf and performs a three-phase parse:
 *
 *   1. **Body section** (before the first @c END card):
 *      Reads body-type keys (SPH, RCC, RPP, …), body names, and raw
 *      argument arrays into @ref osh_geometry_body entries.
 *
 *   2. **Zone section** (between the first and second @c END cards):
 *      Reads zone names and accumulates boolean body-expression strings
 *      into @ref osh_geometry_zone.expr.  No tokenisation or CSG tree
 *      construction happens here — that is the job of
 *      @ref osh_geometry_workspace_prepare().
 *
 *   3. **Material section** (after the second @c END card):
 *      Reads @c ASSIGNMAT (and the legacy positional) cards and stores
 *      material names as strings in @ref osh_geometry_zone.material_name.
 *      No cross-reference against mat.dat is done here.
 *
 * @p ws must be pre-allocated by the caller.  On failure the workspace may
 * be partially populated; the caller is responsible for freeing it via
 * @ref osh_geometry_workspace_free().
 *
 * @param[in]     oshf  Open geometry file positioned anywhere (rewound internally).
 * @param[in,out] ws    Pre-allocated workspace to fill.
 *
 * @returns OSH_OK on success, OSH_EPARSE on format errors, OSH_ENOMEM on
 *          allocation failure.
 */
enum osh_status osh_geometry_parse(struct oshfile *oshf, struct osh_geometry_workspace *ws);

#endif /* OSH_FRONTEND_OPENSHIELDHIT_GEOMETRY_PARSE_H */
