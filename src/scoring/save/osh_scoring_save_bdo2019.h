#ifndef OSH_SCORING_SAVE_BDO2019_H
#define OSH_SCORING_SAVE_BDO2019_H

#include "openshieldhit/status.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save one compiled output in BDO 2019 format.
 *
 * @details
 * BDO 2019 is the current target binary format and should remain compatible
 * with existing external readers. Page data are written exactly as stored in
 * the runtime page buffer (no per-primary normalisation); the actual simulated
 * particle count is embedded in the file as OSHBDO_RT_NSTAT so readers can
 * normalise on load.
 *
 * Data order follows the canonical flat mesh layout:
 *   idx = ix + nx * (iy + ny * iz)
 *
 * TODO: design this writer so a fully postprocessed save call can later be
 * queued to a dedicated CPU thread for asynchronous disk output.
 *
 * @param[in] ws          Scoring workspace with output metadata and file paths.
 * @param[in] rt          Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat       Actual number of primary particles simulated; embedded
 *                        in the BDO file as OSHBDO_RT_NSTAT.
 * @param[in] output_idx  Zero-based index into ws->outputs / rt->outputs.
 */
enum osh_status osh_scoring_save_bdo2019_output(struct osh_scoring_workspace const *ws,
                                                struct osh_scoring_runtime const *rt,
                                                unsigned long long nstat,
                                                size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_BDO2019_H */
