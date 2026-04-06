#ifndef OSH_MATERIAL_PARSE_H
#define OSH_MATERIAL_PARSE_H

#include "common/osh_file.h"
#include "common/osh_rc.h"
#include "material/osh_material.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Parse a material input file into a material workspace.
 *
 * @details
 * Reads key-value lines from @p oshf and appends raw material definitions to
 * @p wm. No derived material properties or transport tables are computed here.
 *
 * @param[in,out] oshf  Open material file positioned at the start.
 * @param[in,out] wm    Workspace to populate.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_material_parse(struct oshfile *oshf, struct material_workspace *wm);

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_PARSE_H */
