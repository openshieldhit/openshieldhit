#ifndef OSH_FRONTEND_OPENSHIELDHIT_SCORING_PARSE_H
#define OSH_FRONTEND_OPENSHIELDHIT_SCORING_PARSE_H

#include "common/osh_file.h"
#include "openshieldhit/scoring.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Parse grouped scoring definitions from an already-open oshfile into
 *        an existing (caller-owned) workspace.
 *
 * @details
 * Pure parser: no file open/close, no workspace allocation or deallocation.
 * On failure the workspace is left in a partially-filled state; the caller is
 * responsible for freeing it.
 *
 * This parser targets the newer grouped `detect.dat` format:
 *
 * - `Filter`
 * - `Settings`
 * - `Geometry <type>`
 * - `Output`
 *
 * The legacy fixed-column SHIELD-HIT Fortran/USRBIN-style format is out of
 * scope here.
 */
enum osh_status osh_scoring_parse(struct oshfile *oshf, struct osh_scoring_workspace *ws);

#ifdef __cplusplus
}
#endif /* OSH_FRONTEND_OPENSHIELDHIT_SCORING_PARSE_H */

#endif
