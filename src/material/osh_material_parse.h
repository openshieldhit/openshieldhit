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
 * Element composition cards use a natural element when only Z is provided:
 * `NUCLID <Z> <amount>` stores A=0. An explicit isotope is selected by adding
 * a mass number: `NUCLID <Z> <A> <amount>`. The same Z/[A]/amount rule applies
 * to ELEMENT, ELEMENTBYNUMBER, and ELEMENTBYMASS.
 *
 * MEE cards are parsed raw here and resolved later during material assembly.
 * For compounds, users may specify either one material-level MEE or element-
 * level MEE values, but not both in the same material definition.
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
