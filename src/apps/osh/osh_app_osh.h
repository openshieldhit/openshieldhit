#ifndef OSH_APP_OSH_H
#define OSH_APP_OSH_H

#include "openshieldhit/beam.h"
#include "openshieldhit/diag.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"

/**
 * @brief Load, parse, and prepare a beam workspace from @p path.
 *
 * @details
 * Opens @p path, parses it into a cold beam workspace, loads any referenced
 * spot list (USECBEAM card), and calls @ref osh_beam_workspace_prepare(). On
 * success the caller owns the workspace and must release it with
 * @ref osh_beam_workspace_free().
 *
 * @param[in]  path    Path to the beam input file (beam.dat).
 * @param[in]  diag    Borrowed diagnostics sink for parser/setup messages, or NULL.
 * @param[out] wb_out  Receives the prepared workspace.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status
osh_beam_setup_from_path(char const *path, struct osh_diag_sink const *diag, struct osh_beam_workspace **wb_out);

/**
 * @brief Load, parse, and prepare a geometry workspace from @p path.
 *
 * @details
 * Opens @p path, parses body/zone/material-assignment sections into a cold
 * geometry workspace, and calls @ref osh_geometry_workspace_prepare(). On
 * success the caller owns the workspace and must release it with
 * @ref osh_geometry_workspace_free().
 *
 * @param[in]  path    Path to the geometry input file (geo.dat).
 * @param[in]  diag    Borrowed diagnostics sink for parser/setup messages, or NULL.
 * @param[out] ws_out  Receives the prepared workspace.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_geometry_setup_from_path(char const *path,
                                             struct osh_diag_sink const *diag,
                                             struct osh_geometry_workspace **ws_out);

/**
 * @brief Load, parse, and prepare a material workspace from @p path.
 *
 * @details
 * Opens @p path, parses material definitions into a cold material workspace,
 * and calls @ref osh_material_workspace_prepare(). On success the caller owns
 * the workspace and must release it with @ref osh_material_workspace_free().
 *
 * @param[in]  path    Path to the material input file (mat.dat).
 * @param[in]  diag    Borrowed diagnostics sink for parser/setup messages, or NULL.
 * @param[out] wm_out  Receives the prepared workspace.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_material_setup_from_path(char const *path,
                                             struct osh_diag_sink const *diag,
                                             struct osh_material_workspace **wm_out);

/**
 * @brief Load and parse a scoring workspace from @p path.
 *
 * @details
 * Opens @p path and parses scoring detector definitions into a cold scoring
 * workspace. There is no prepare step for scoring; the workspace is ready for
 * @ref osh_scoring_compile() after this call. On success the caller owns the
 * workspace and must release it with @ref osh_scoring_workspace_free().
 *
 * @param[in]  path    Path to the detect input file (detect.dat).
 * @param[in]  diag    Borrowed diagnostics sink for parser/setup messages, or NULL.
 * @param[out] ws_out  Receives the workspace.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status
osh_scoring_setup_from_path(char const *path, struct osh_diag_sink const *diag, struct osh_scoring_workspace **ws_out);

/**
 * @brief Resolve `Geometry Zone` selectors (zone names) to transport zone indices.
 *
 * @details
 * The scoring parser stores `Zone <name>` selectors verbatim because it has no
 * access to geo.dat. This app-level step maps each name to its 0-based transport
 * zone index using @p geom's zone table (populating @c zone_indices for every
 * Zone scoring geometry), so the library scoring/compile path only ever sees
 * resolved indices — never names. Must be called after both geo.dat and
 * detect.dat are parsed and before @ref osh_scoring_compile().
 *
 * An unknown zone name is a hard error. A zone with no `Volume` card is a warning
 * (volume-normalized quantities on that zone default to 1.0 cm3).
 *
 * @param[in,out] scoring  Scoring workspace whose Zone geometries are resolved.
 * @param[in]     geom     Prepared geometry workspace providing the zone-name table.
 * @param[in]     diag     Borrowed diagnostics sink, or NULL.
 *
 * @returns OSH_OK on success; OSH_EINVAL on an unknown zone name; OSH_ENOMEM on
 *          allocation failure.
 */
enum osh_status osh_scoring_resolve_zone_names(struct osh_scoring_workspace *scoring,
                                               struct osh_geometry_workspace const *geom,
                                               struct osh_diag_sink const *diag);

#endif /* OSH_APP_OSH_H */
