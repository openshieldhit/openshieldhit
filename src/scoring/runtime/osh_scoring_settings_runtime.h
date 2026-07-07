#ifndef OSH_SCORING_SETTINGS_RUNTIME_H
#define OSH_SCORING_SETTINGS_RUNTIME_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compiled settings block, mirroring the parsed form.
 *
 * @details
 * Optional fields are guarded by `has_*` flags; a field should only be applied
 * when its flag is non-zero.  This mirrors the parsed @ref osh_scoring_settings_def
 * but lives in the runtime layer so the hot path does not need to touch the
 * parse workspace.
 */
struct osh_scoring_settings_runtime {
    char *name;              /* Settings name (owned). */
    double rescale;          /* Multiplicative output rescale factor. */
    double offset;           /* Additive output offset. */
    double site_diameter_um; /* Site diameter for microdosimetry [um]. */
    double density_g_cm3;    /* Local density override [g/cm^3]. */
    size_t npart;            /* Maximum particle count cap. */
    int medium;              /* Transport medium override index. */
    int nkmedium;            /* Neutron-kerma medium override index. */
    int variance;            /* Monte-Carlo standard-error tracking: 1 = on, 0 = off. */
    char has_rescale;
    char has_offset;
    char has_site_diameter_um;
    char has_density_g_cm3;
    char has_npart;
    char has_medium;
    char has_nkmedium;
    char has_variance;
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SETTINGS_RUNTIME_H */
