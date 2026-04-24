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
 * @brief Parse a geometry workspace from @p path without preparing it.
 *
 * @details
 * Opens @p path and parses body/zone/material-assignment sections into a cold
 * geometry workspace.  @ref osh_geometry_workspace_prepare() is NOT called; the
 * caller must do so explicitly after setting any fields that prepare reads (e.g.
 * @c hu_table_type). On success the caller owns the workspace and must release
 * it with @ref osh_geometry_workspace_free().
 *
 * @param[in]  path    Path to the geometry input file (geo.dat).
 * @param[in]  diag    Borrowed diagnostics sink for parser messages, or NULL.
 * @param[out] ws_out  Receives the unprepared workspace.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_geometry_parse_from_path(char const *path,
                                             struct osh_diag_sink const *diag,
                                             struct osh_geometry_workspace **ws_out);

/**
 * @brief Load, parse, and prepare a geometry workspace from @p path.
 *
 * @details
 * Convenience wrapper around @ref osh_geometry_parse_from_path() followed by
 * @ref osh_geometry_workspace_prepare().  Suitable for non-CT geometry where
 * @c hu_table_type does not need to be propagated before prepare. On success
 * the caller owns the workspace and must release it with
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
 * @brief Parse a material workspace from @p path without preparing it.
 *
 * @details
 * Opens @p path, parses material definitions into a cold material workspace,
 * and returns without calling @ref osh_material_workspace_prepare(). The caller
 * must call prepare explicitly when ready. On success the caller owns the
 * workspace and must release it with @ref osh_material_workspace_free().
 *
 * @param[in]  path    Path to the material input file (mat.dat).
 * @param[in]  diag    Borrowed diagnostics sink for parser messages, or NULL.
 * @param[out] wm_out  Receives the unprepared workspace.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_material_parse_from_path(char const *path,
                                             struct osh_diag_sink const *diag,
                                             struct osh_material_workspace **wm_out);

/**
 * @brief Load, parse, and prepare a material workspace from @p path.
 *
 * @details
 * Convenience wrapper around @ref osh_material_parse_from_path() followed by
 * @ref osh_material_workspace_prepare(). On success the caller owns the
 * workspace and must release it with @ref osh_material_workspace_free().
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

#endif /* OSH_APP_OSH_H */
