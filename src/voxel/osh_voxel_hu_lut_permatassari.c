#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/material.h"
#include "voxel/osh_voxel_hu_lut.h"
#include "voxel/osh_voxel_mat_permatassari2020.h"

static uint8_t hu_to_bin_from_breakpoints(int16_t hu, int16_t const *breakpoints, size_t nbreakpoints);

/* ---- Permatassari 2020 --------------------------------------------------- */

/**
 * @brief Register Permatassari 2020 tissue bins as material workspace entries.
 *
 * @details
 * Appends 40 material definitions (named "permatassari2020_00" … "_39") to
 * @p wm.  The representative density is derived from the per-bin linear
 * density factor at the bin midpoint HU; mean excitation energy is taken from
 * the published table.
 *
 * @param[in,out] wm  Material workspace to extend.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_voxel_register_permatassari_materials(struct osh_material_workspace *wm) {
    struct osh_material *materials;
    struct osh_material *mat;
    struct osh_material_element *elements;
    size_t first_idx;
    size_t i;
    size_t j;
    size_t k;
    size_t nelm;
    double frac_total;
    double hu_mid;
    double rho_nominal;
    char name[24];

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
        /* Nominal density at bin midpoint HU, from the Permatassari linear model. */
        hu_mid = 0.5 * ((double) _ct_hu[i] + (double) _ct_hu[i + 1]);
        rho_nominal = _ct_density_factor[i] * (1000.0 + hu_mid);
        if (rho_nominal < 0.0001) {
            rho_nominal = 0.0001;
        }
        mat->rho = rho_nominal;
        mat->mean_excitation_energy = _ct_mean_excitation_energy[i];
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

        snprintf(name, sizeof(name), "permatassari2020_%02zu", i);
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
 * @brief Map a raw HU value to a Permatassari 2020 density [g/cm³].
 *
 * @details
 * Applies the linear density model @c rho = density_factor[bin] * (1000 + hu)
 * with a lower clamp of 0.0001 g/cm³.  Both @p hu and @p bin are clamped to
 * their valid ranges before use.
 *
 * @param[in] hu   Hounsfield unit; clamped to [-1000, 1600].
 * @param[in] bin  Calibration bin index; clamped to [0, _nmat-1].
 *
 * @returns Density [g/cm³].
 */
float osh_voxel_hu2rho_permatassari2020(int16_t hu, int bin) {
    double rho;

    if (hu < -1000) {
        hu = -1000;
    }
    if (hu > 1600) {
        hu = 1600;
    }
    if (bin < 0) {
        bin = 0;
    }
    if (bin >= _nmat) {
        bin = _nmat - 1;
    }

    rho = _ct_density_factor[bin] * (1000.0 + (double) hu);
    if (rho < 0.0001) {
        rho = 0.0001;
    }

    return (float) rho;
}

/**
 * @brief Fill the HU→bin-index LUT for the Permatassari 2020 scheme.
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_bin_lut_permatassari2020(uint8_t lut[OSH_VOXEL_HU_LUT_SIZE]) {
    int hu;

    hu = -1000;
    while (hu <= 1600) {
        lut[hu + 1000] = hu_to_bin_from_breakpoints((int16_t) hu, _ct_hu, (size_t) _nmat + 1u);
        hu++;
    }
}

/**
 * @brief Fill the HU→density LUT for the Permatassari 2020 scheme.
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_rho_lut_permatassari2020(float lut[OSH_VOXEL_HU_LUT_SIZE]) {
    int hu;
    int bin;

    hu = -1000;
    while (hu <= 1600) {
        bin = (int) hu_to_bin_from_breakpoints((int16_t) hu, _ct_hu, (size_t) _nmat + 1u);
        lut[hu + 1000] = osh_voxel_hu2rho_permatassari2020((int16_t) hu, bin);
        hu++;
    }
}

/* ---- Private helpers ----------------------------------------------------- */

/*
 * Return the bin index for @p hu given a sorted breakpoint array.
 * Scans left-to-right; the result is the last bin whose lower edge is <= hu.
 */
static uint8_t hu_to_bin_from_breakpoints(int16_t hu, int16_t const *breakpoints, size_t nbreakpoints) {
    size_t bin;

    bin = 0u;
    while ((bin + 1u) < nbreakpoints && hu >= breakpoints[bin + 1u]) {
        bin++;
    }
    if (bin >= nbreakpoints - 1u) {
        bin = nbreakpoints - 2u;
    }

    return (uint8_t) bin;
}
