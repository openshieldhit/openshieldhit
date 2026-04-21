#ifndef OSH_FRONTEND_OPENSHIELDHIT_BEAM_PARSE_H
#define OSH_FRONTEND_OPENSHIELDHIT_BEAM_PARSE_H

#include "common/osh_file.h"
#include "openshieldhit/beam.h"
#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"

/**
 * @brief Parse a beam.dat file into a beam workspace.
 *
 * Reads key-value lines from @p oshf and fills raw workspace fields plus one
 * caller-owned template spot. No derived quantities are computed here — that
 * is the job of the post-parse step in osh_beam_setup_from_path().
 *
 * Recognised keys include PRIMARY, TMAX0, BEAMPOS, BEAMSIGMA, NSTAT,
 * and others listed in osh_beam_parse_keys.h. Unknown keys and malformed
 * values both log a diagnostic and return an OSH_E* code to the
 * caller.
 *
 * @param[in]     oshf  Open file handle positioned at the start of the beam
 *                      data.  filename and lineno are used in diagnostics.
 * @param[in]     diag              Borrowed diagnostics sink for parse messages, or NULL.
 * @param[in,out] beam              Workspace to fill.  Must be pre-allocated
 *                                  and zero-initialised by the caller.
 * @param[out]    spot_out          Receives the parsed single-spot template.
 *                                  Spot-list imports may later use it as the
 *                                  fallback source for omitted columns.
 * @param[out]    spotlist_path_out Receives an owned, resolved USECBEAM path
 *                                  when that card is present. May be NULL if
 *                                  the caller does not need spot-list loading.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_beam_parse(struct oshfile *oshf,
                               struct osh_diag_sink const *diag,
                               struct osh_beam_workspace *beam,
                               struct osh_beam_spot *spot_out,
                               char **spotlist_path_out);

#endif /* OSH_FRONTEND_OPENSHIELDHIT_BEAM_PARSE_H */
