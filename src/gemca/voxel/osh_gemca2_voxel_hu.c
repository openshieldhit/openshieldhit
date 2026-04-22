#include "gemca/voxel/osh_gemca2_voxel_hu.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_interpolate.h"
#include "gemca/voxel/osh_gemca2_voxel_defines.h"
#include "gemca/voxel/osh_voxel_mat_schneider2000.h"
#include "openshieldhit/material.h"

static inline float _wepl_minohara1993(int16_t hu);
static inline float _wepl_jacob1996(int16_t hu);
static inline float _wepl_geiss1999(int16_t hu);
static int _hu2idx_slow(int16_t hu);

/*
 * HU → density [g/cm³], Schneider 2000, Eqs. 20–23.
 *
 * This continuous piecewise-linear function gives the *actual* voxel density
 * used for stopping-power density scaling during transport:
 *
 *   dE/dx_actual = (rho_actual / rho_nominal) * dE/dx_table
 *
 * where rho_nominal is the representative bin density from _ct_hu_rho[] and
 * dE/dx_table is the stopping power tabulated for that material bin.
 *
 * The function is NOT the same as _ct_hu_rho[]: those are discrete bin
 * averages used to characterise material composition; this function varies
 * continuously within a bin.
 */
float osh_gemca_voxel_hu2rho(int16_t hu, char alg) {
    float ret = 0.0;

    (void) (alg); /* unused; reserved for future algorithm selection */

    if (hu < -1000) {
        hu = -1000;
    }
    if (hu > 1600) {
        hu = 1600;
    }

    /* Schneider 2000 breakpoints: -98, 14, 23, 100 */
    if (hu <= -98) {
        ret = 1.03091f + 0.0010297f * hu; /* Eq. 20: air / lung */
    } else if (hu <= 14) {
        ret = 1.018f + 0.000893f * hu; /* Eq. 21: adipose–water */
    } else if (hu <= 23) {
        ret = 1.03f; /* Eq. 22: water transition */
    } else if (hu <= 100) {
        ret = 1.003f + 0.001169f * hu; /* Eq. 23: soft tissue */
    } else {
        ret = 1.017f + 0.000592f * hu; /* Eq. 24: bone */
    }
    return ret;
}

int osh_gemca_voxel_hu2idx(int16_t hu) {
    if (hu < -1000) {
        return 0;
    }
    if (hu > 1600) {
        return 0;
    }
    return _hu2idx_slow(hu);
}

/*
 * WEPL (water-equivalent path length) conversion functions.
 *
 * These are NOT used in the main transport path. Range calculation in
 * transport uses hu2rho() with stopping-power density scaling instead.
 * The WEPL functions are provided as optional utilities for external use
 * (e.g. range probing, treatment planning cross-checks).
 */
float osh_gemca_voxel_hu2wepl(int16_t hu, char alg) {

    float wepl = 0.0;

    if (hu < -1000) {
        return 0.0;
    }

    if (hu > 4000) {
        return 0.0;
    }

    switch (alg) {
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG1:
        wepl = _wepl_minohara1993(hu);
        break;
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG2:
        wepl = _wepl_jacob1996(hu);
        break;
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG3:
        wepl = _wepl_geiss1999(hu);
        break;
    default:
        return 0.0;
        break;
    }
    return wepl * 1000.0;
}

enum osh_status osh_gemca_voxel_register_schneider_materials(struct osh_material_workspace *wm) {
    struct osh_material *materials;
    struct osh_material *mat;
    struct osh_material_element *elements;
    size_t first_idx;
    size_t i;
    size_t j;
    size_t k;
    size_t nelm;
    double frac_total;
    char name[16];

    if (!wm) {
        return OSH_EINVAL;
    }

    first_idx = wm->nmaterials;
    materials = (struct osh_material *) realloc(wm->materials, (first_idx + (size_t) _nmat) * sizeof(*materials));
    if (!materials) {
        return OSH_ENOMEM;
    }
    wm->materials = materials;
    wm->nmaterials = first_idx + (size_t) _nmat;

    i = 0;
    while (i < (size_t) _nmat) {
        mat = &wm->materials[first_idx + i];

        nelm = 0;
        j = 0;
        while (j < (size_t) _nelm) {
            if (_ct_relm[i][j] > 0.0f && (unsigned int) _ct_elmz[j] > 0u) {
                nelm++;
            }
            j++;
        }

        mat->elements = NULL;
        mat->dedx_overrides = NULL;
        mat->name = NULL;
        mat->rho = (double) _ct_hu_rho[i];
        mat->mean_excitation_energy = -1.0;
        mat->rgba[0] = 0.8f;
        mat->rgba[1] = 0.8f;
        mat->rgba[2] = 0.8f;
        mat->rgba[3] = 1.0f;
        mat->nelements = 0u;
        mat->ndedx_overrides = 0u;
        mat->lineno = 0u;
        mat->index = first_idx + i;
        mat->icru_id = 0;
        mat->state = (i == 0u) ? OSH_MATERIAL_STATE_GAS : OSH_MATERIAL_STATE_CONDENSED;

        snprintf(name, sizeof(name), "schneider_%02zu", i);
        mat->name = (char *) malloc(strlen(name) + 1u);
        if (!mat->name) {
            wm->nmaterials = first_idx + i;
            return OSH_ENOMEM;
        }
        memcpy(mat->name, name, strlen(name) + 1u);

        elements = (struct osh_material_element *) calloc(nelm, sizeof(*elements));
        if (!elements) {
            free(mat->name);
            mat->name = NULL;
            wm->nmaterials = first_idx + i;
            return OSH_ENOMEM;
        }

        /* sum active fractions for normalization (table values may not sum to 100.0 exactly) */
        frac_total = 0.0;
        j = 0;
        while (j < (size_t) _nelm) {
            if (_ct_relm[i][j] > 0.0f && (unsigned int) _ct_elmz[j] > 0u) {
                frac_total += (double) _ct_relm[i][j];
            }
            j++;
        }

        k = 0;
        j = 0;
        while (j < (size_t) _nelm) {
            if (_ct_relm[i][j] > 0.0f && (unsigned int) _ct_elmz[j] > 0u) {
                elements[k].z = (unsigned int) _ct_elmz[j];
                elements[k].a = 0u;
                elements[k].mass_fraction = (double) _ct_relm[i][j] / frac_total;
                elements[k].atom_count = -1.0;
                elements[k].mean_excitation_energy = -1.0;
                elements[k].lineno = 0u;
                k++;
            }
            j++;
        }

        mat->elements = elements;
        mat->nelements = nelm;
        i++;
    }

    return OSH_OK;
}

void osh_gemca_voxel_build_hu_lut(uint8_t lut[2601]) {
    int hu;

    hu = -1000;
    while (hu <= 1600) {
        lut[hu + 1000] = (uint8_t) _hu2idx_slow((int16_t) hu);
        hu++;
    }
}

static int _hu2idx_slow(int16_t hu) {
    return (int) osh_binary_search_i2(hu, _ct_hu, (unsigned long int) (_nmat + 1));
}

/* Optional WEPL conversion algorithms — not used in the transport path */

static inline float _wepl_minohara1993(int16_t hu) {
    if (hu < -49) {
        return 1.075e-3 * hu + 1.050;
    } else {
        return 4.597e-4 * hu + 1.019;
    }
}

static inline float _wepl_jacob1996(int16_t hu) {
    if (hu < -60.81) { /* float to int comparison */
        return 1.011e-3 * hu + 1.052;
    } else {
        return 4.190e-4 * hu + 1.016;
    }
}

static inline float _wepl_geiss1999(int16_t hu) {
    /* TODO */
    (void) (hu);
    return 0.0;
}
