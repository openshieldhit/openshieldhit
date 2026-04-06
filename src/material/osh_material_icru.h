#ifndef OSH_MATERIAL_ICRU_H
#define OSH_MATERIAL_ICRU_H

#include <stddef.h>

#include "common/osh_rc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Physical state of a material.
 *
 * The ICRU database is the authoritative source for gas/condensed classification,
 * so the enum lives here rather than in osh_material.h.
 */
enum osh_material_state { OSH_MATERIAL_STATE_UNSET = 0, OSH_MATERIAL_STATE_CONDENSED = 1, OSH_MATERIAL_STATE_GAS = 2 };

#define OSH_MATERIAL_ICRU_MAX_ELEMENTS 16u

struct osh_material_icru_element {
    double mass_fraction;
    unsigned int z;
    unsigned int a;
};

struct osh_material_icru_entry {
    struct osh_material_icru_element elements[OSH_MATERIAL_ICRU_MAX_ELEMENTS];
    double rho;
    double mean_excitation_energy;
    size_t nelements;
    int icru_id;
    int state;
};

/**
 * @brief Look up an embedded ICRU material definition.
 *
 * @param[in]  icru_id    ICRU material id.
 * @param[out] entry_out  Receives the immutable DB values copied into caller-owned storage.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_material_icru_lookup(int icru_id, struct osh_material_icru_entry *entry_out);

/**
 * @brief Return the ICRU compound element mean excitation energy when tabulated.
 *
 * @details
 * ICRU elemental entries store pure-element mean excitation energies. ICRU 49
 * Table 2.11 also lists state-dependent values for selected elements when they
 * occur in compounds. Material assembly uses those tabulated values only for
 * unset element-level MEE fields; explicit user IVALUE/IAV/ELEMENTI cards still
 * take precedence and are never overwritten.
 *
 * @param[in] z                 Atomic number.
 * @param[in] state             enum osh_material_state value.
 * @param[in] pure_element_mee  Pure-element mean excitation energy [eV].
 *
 * @returns Compound element MEE [eV] if tabulated, otherwise @p pure_element_mee.
 */
double osh_material_icru_compound_element_mean_excitation_energy(unsigned int z, int state, double pure_element_mee);

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_ICRU_H */
