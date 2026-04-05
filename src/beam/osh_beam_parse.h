#ifndef OSH_BEAM_PARSE_H
#define OSH_BEAM_PARSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "beam/osh_beam.h"
#include "common/osh_file.h"

/**
 * @brief Parse a beam.dat file into a beam workspace.
 *
 * Reads key-value lines from @p oshf and fills raw fields in @p beam.
 * No derived quantities are computed here — that is the job of the
 * post-parse step in osh_beam_setup_from_path().
 *
 * Recognised keys include PRIMARY, TMAX0, BEAMPOS, BEAMSIGMA, NSTAT,
 * and others listed in osh_beam_parse_keys.h.  Unknown keys produce a
 * warning and are skipped; malformed values call osh_error() and do not
 * return.
 *
 * @param[in]     oshf  Open file handle positioned at the start of the beam
 *                      data.  filename and lineno are used in diagnostics.
 * @param[in,out] beam  Workspace to fill.  Must be pre-allocated and
 *                      zero-initialised by the caller.
 *
 * @returns OSH_OK on success.
 */
int osh_beam_parse(struct oshfile *oshf, struct beam_workspace *beam);

#endif /* OSH_BEAM_PARSE_H */