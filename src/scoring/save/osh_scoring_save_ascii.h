#ifndef OSH_SCORING_SAVE_ASCII_H
#define OSH_SCORING_SAVE_ASCII_H

#include "common/osh_rc.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

enum osh_status osh_scoring_save_ascii_output(struct osh_scoring_save_request const *req, size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_ASCII_H */
