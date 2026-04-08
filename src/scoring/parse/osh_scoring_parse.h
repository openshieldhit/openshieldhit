#ifndef OSH_SCORING_PARSE_H
#define OSH_SCORING_PARSE_H

#include "scoring/osh_scoring.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Parse grouped scoring definitions from an already-open file path.
 *
 * @details
 * This parser currently targets the newer grouped `detect.dat` format used by
 * `tests/fixtures/test01/detect.dat`:
 *
 * - `Filter`
 * - `Settings`
 * - `Geometry <type>`
 * - `Output`
 *
 * The legacy fixed-column detector format from SHIELD-HIT Fortran/USRBIN-style
 * `detect.dat` files is intentionally out of scope here.
 *
 * It intentionally produces only a raw parsed workspace. Runtime scorer
 * compilation and BDO output metadata setup should happen in a later phase.
 */
enum osh_status osh_scoring_parse_file(char const *path, struct osh_scoring_workspace **ws_out);

#ifdef __cplusplus
}
#endif

#endif
