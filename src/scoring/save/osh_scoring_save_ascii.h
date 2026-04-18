#ifndef OSH_SCORING_SAVE_ASCII_H
#define OSH_SCORING_SAVE_ASCII_H

#include "openshieldhit/status.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save one compiled output in plain-text (ASCII) format.
 *
 * @details
 * Writes a human-readable columnar file with a version header line followed by
 * one data line per bin. Each line contains the bin-centre coordinates and the
 * scored value normalised per primary (value / @p nstat).
 *
 * @param[in] ws          Scoring workspace with output metadata and file paths.
 * @param[in] rt          Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat       Actual number of primary particles simulated; used for
 *                        per-primary normalisation.
 * @param[in] output_idx  Zero-based index into ws->outputs / rt->outputs.
 */
enum osh_status osh_scoring_save_ascii_output(struct osh_scoring_workspace const *ws,
                                              struct osh_scoring_runtime const *rt,
                                              unsigned long long nstat,
                                              size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_ASCII_H */
