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
 * one data line per bin.  Each line contains the bin-centre coordinates and the
 * scored value, normalised by @p req->nstat when @p req->has_nstat is set.
 */
enum osh_status osh_scoring_save_ascii_output(struct osh_scoring_save_request const *req, size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_ASCII_H */
