#ifndef OSH_MATERIAL_ATOMIC_DATA_H
#define OSH_MATERIAL_ATOMIC_DATA_H

#include "common/osh_rc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Look up the natural atomic weight for an element.
 *
 * @details
 * The returned value is a standard atomic weight for the natural element, in
 * Dalton (Da). This is intentionally separate from isotope masses in
 * particle/, which are used when a material element has an explicit mass number
 * A > 0.
 *
 * @param[in]  z         Atomic number Z.
 * @param[out] mass_out  Receives natural atomic weight [Da].
 *
 * @returns OSH_OK on success, or OSH_EINVAL if Z is outside the table.
 */
enum osh_status osh_material_natural_atomic_mass_da(unsigned int z, double *mass_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_ATOMIC_DATA_H */
