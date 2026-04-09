#ifndef OSH_SCORING_SAVE_BDO2019_H
#define OSH_SCORING_SAVE_BDO2019_H

#include "common/osh_rc.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save one compiled output in BDO 2019 format.
 *
 * @details
 * BDO 2019 is the current target binary format and should remain compatible
 * with existing external readers. Page data must be serialized in the
 * canonical flat runtime order without reshaping:
 *
 *   idx = ix + nx * (iy + ny * iz)
 *
 * That matches the original SHIELD-HIT mesh layout and keeps the on-disk
 * representation stable for existing readers.
 *
 * TODO: design this writer so a fully postprocessed save request can later be
 * queued to a dedicated CPU thread for asynchronous disk output.
 */
enum osh_status osh_scoring_save_bdo2019_output(struct osh_scoring_save_request const *req, size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_BDO2019_H */
