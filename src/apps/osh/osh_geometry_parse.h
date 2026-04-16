#ifndef OSH_FRONTEND_OPENSHIELDHIT_GEOMETRY_PARSE_H
#define OSH_FRONTEND_OPENSHIELDHIT_GEOMETRY_PARSE_H

#include "openshieldhit/geometry.h"
#include "openshieldhit/status.h"

/**
 * @brief Parse a geo.dat file into a cold geometry workspace.
 *
 * @details
 * Opens @p path and performs a three-phase parse:
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
 * On success, @p *ws_out is populated with a newly allocated workspace
 * that the caller owns and must eventually free with
 * @ref osh_geometry_workspace_free().  The workspace is not yet prepared;
 * call @ref osh_geometry_workspace_prepare() before passing it to the
 * geometry runtime.
 *
 * On failure, @p *ws_out is set to NULL and the partially built workspace
 * is freed internally.
 *
 * @param[in]  path    Path to the geo.dat file.
 * @param[out] ws_out  Receives the newly allocated cold workspace.
 *
 * @returns OSH_OK on success, OSH_EIO if the file cannot be opened,
 *          OSH_EPARSE on format errors, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_geometry_parse_file(char const *path, struct osh_geometry_workspace **ws_out);

#endif /* OSH_FRONTEND_OPENSHIELDHIT_GEOMETRY_PARSE_H */
