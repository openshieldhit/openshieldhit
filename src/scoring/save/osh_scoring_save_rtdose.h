#ifndef OSH_SCORING_SAVE_RTDOSE_H
#define OSH_SCORING_SAVE_RTDOSE_H

#include "openshieldhit/status.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save one compiled output in DICOM RTDOSE format.
 *
 * @details
 * Reads the RTDOSE template stored in the runtime geometry's
 * rtdose_template_path, overlays the scored data (normalised by @p nstat),
 * and writes a modified copy to the output filename.  All DICOM metadata is
 * preserved unchanged; only the pixel data is replaced.
 *
 * The output stores raw scored values divided by nstat and by the template's
 * dose_grid_scaling factor, yielding uint32_t pixel values.  For ENERGY
 * scoring the result is energy-per-primary in dose_grid_scaling units, not
 * absorbed dose; proper DOSE scoring (with density weighting) is deferred.
 *
 * Only single-page outputs are supported.  Multi-page RTDOSE outputs return
 * OSH_ENOTSUP.
 *
 * @param[in] ws          Scoring workspace with output metadata and file paths.
 * @param[in] rt          Compiled scoring runtime with accumulated data.
 * @param[in] nstat       Actual number of primary particles simulated; must be > 0.
 * @param[in] output_idx  Zero-based index into ws->outputs / rt->outputs.
 */
enum osh_status osh_scoring_save_rtdose_output(struct osh_scoring_workspace const *ws,
                                               struct osh_scoring_runtime const *rt,
                                               unsigned long long nstat,
                                               size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_RTDOSE_H */
