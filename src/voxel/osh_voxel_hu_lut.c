#include "voxel/osh_voxel_hu_lut.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_interpolate.h"
#include "gemca/voxel/osh_gemca2_voxel_defines.h"
#include "openshieldhit/material.h"
#include "voxel/osh_voxel_mat_schneider2000.h"

static int _hu2idx_slow(int16_t hu);
static inline float _wepl_minohara1993(int16_t hu);
static inline float _wepl_jacob1996(int16_t hu);

/* ---- Schneider 2000 ------------------------------------------------------ */

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
float osh_voxel_hu2rho_schneider2000(int16_t hu) {
    float ret;

    if (hu < -1000) {
        hu = -1000;
    }
    if (hu > 1600) {
        hu = 1600;
    }

    /* Schneider 2000 breakpoints: -98, 14, 23, 100 */
    if (hu <= -98) {
        ret = 1.03091f + 0.0010297f * (float) hu; /* Eq. 20: air / lung */
    } else if (hu <= 14) {
        ret = 1.018f + 0.000893f * (float) hu; /* Eq. 21: adipose–water */
    } else if (hu <= 23) {
        ret = 1.03f; /* Eq. 22: water transition */
    } else if (hu <= 100) {
        ret = 1.003f + 0.001169f * (float) hu; /* Eq. 23: soft tissue */
    } else {
        ret = 1.017f + 0.000592f * (float) hu; /* Eq. 24: bone */
    }
    return ret;
}

int osh_voxel_hu2idx_schneider2000(int16_t hu) {
    if (hu < -1000) {
        return 0;
    }
    if (hu > 1600) {
        return _nmat - 1;
    }
    return _hu2idx_slow(hu);
}

/**
 * @brief Register Schneider 2000 tissue bins as material workspace entries.
 *
 * @details
 * Appends 24 material definitions to @p wm, named "schneider_00" through
 * "schneider_23".  The representative density, elemental composition, and
 * physical state of each bin are taken directly from the Schneider 2000 table.
 *
 * @param[in,out] wm  Material workspace to extend.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_voxel_register_schneider_materials(struct osh_material_workspace *wm) {
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
            if (_ct_relm[i][j] > 0.0f && _ct_elmz[j] > 0u) {
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

        /* Sum active fractions for normalisation (table values may not sum to 100.0 exactly). */
        frac_total = 0.0;
        j = 0;
        while (j < (size_t) _nelm) {
            if (_ct_relm[i][j] > 0.0f && _ct_elmz[j] > 0u) {
                frac_total += (double) _ct_relm[i][j];
            }
            j++;
        }

        k = 0;
        j = 0;
        while (j < (size_t) _nelm) {
            if (_ct_relm[i][j] > 0.0f && _ct_elmz[j] > 0u) {
                elements[k].z = _ct_elmz[j];
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

/**
 * @brief Fill the HU→bin-index LUT for the Schneider 2000 scheme.
 *
 * @details
 * Fills @p lut[hu+1000] for hu in [-1000, 1600] with the material bin index
 * according to the Schneider 2000 breakpoints.  Used by the geometry runtime.
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_bin_lut_schneider2000(uint8_t lut[OSH_VOXEL_HU_LUT_SIZE]) {
    int hu;

    hu = -1000;
    while (hu <= 1600) {
        lut[hu + 1000] = (uint8_t) _hu2idx_slow((int16_t) hu);
        hu++;
    }
}

/**
 * @brief Fill the HU→density LUT for the Schneider 2000 scheme.
 *
 * @details
 * Fills @p lut[hu+1000] for hu in [-1000, 1600] with the voxel density
 * [g/cm³] from the Schneider 2000 piecewise-linear fit.  Used by the
 * material runtime.
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_rho_lut_schneider2000(float lut[OSH_VOXEL_HU_LUT_SIZE]) {
    int hu;

    hu = -1000;
    while (hu <= 1600) {
        lut[hu + 1000] = osh_voxel_hu2rho_schneider2000((int16_t) hu);
        hu++;
    }
}

/* ---- WEPL conversions ---------------------------------------------------- */

/**
 * @brief Convert a Hounsfield unit to a water-equivalent path-length ratio.
 *
 * @details
 * Not used in the main Monte Carlo transport path.  Provided as an optional
 * utility for external callers.  @p alg selects among three empirical
 * conversion formulae (1=Minohara 1993, 2=Jacob 1996, 3=Geiss 1999 stub).
 *
 * @param[in] hu   Hounsfield unit; returns 0 outside [-1000, 4000].
 * @param[in] alg  Algorithm selector.
 *
 * @returns WEPL multiplied by 1000 (per the legacy convention), or 0.
 */
float osh_voxel_hu2wepl(int16_t hu, char alg) {
    float wepl;

    if (hu < -1000 || hu > 4000) {
        return 0.0;
    }

    switch (alg) {
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG1:
        wepl = _wepl_minohara1993(hu);
        break;
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG2:
        wepl = _wepl_jacob1996(hu);
        break;
    default:
        return 0.0; /* ALG3 (Geiss 1999) is not yet implemented */
    }

    return wepl * 1000.0f;
}

/* ---- Private helpers ----------------------------------------------------- */

static int _hu2idx_slow(int16_t hu) {
    return (int) osh_binary_search_i2(hu, _ct_hu, (unsigned long int) _nmat + 1ul);
}

static inline float _wepl_minohara1993(int16_t hu) {
    if (hu < -49) {
        return 1.075e-3f * (float) hu + 1.050f;
    }
    return 4.597e-4f * (float) hu + 1.019f;
}

static inline float _wepl_jacob1996(int16_t hu) {
    if ((float) hu < -60.81f) {
        return 1.011e-3f * (float) hu + 1.052f;
    }
    return 4.190e-4f * (float) hu + 1.016f;
}
