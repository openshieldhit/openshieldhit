#ifndef OSH_SCORING_SETTINGS_RUNTIME_H
#define OSH_SCORING_SETTINGS_RUNTIME_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct osh_scoring_settings_runtime {
    char *name;
    double rescale;
    double offset;
    double site_diameter_um;
    double density_g_cm3;
    size_t npart;
    int medium;
    int nkmedium;
    char has_rescale;
    char has_offset;
    char has_site_diameter_um;
    char has_density_g_cm3;
    char has_npart;
    char has_medium;
    char has_nkmedium;
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SETTINGS_RUNTIME_H */
