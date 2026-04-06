#ifndef OSH_MATERIAL_ICRU_H
#define OSH_MATERIAL_ICRU_H

#include <stddef.h>

#include "common/osh_rc.h"

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_ICRU_H */
