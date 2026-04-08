#include "material/osh_material.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "material/osh_material_atomic_data.h"
#include "material/osh_material_icru.h"
#include "material/osh_material_parse.h"
#include "particle/osh_isotope_db.h"

static void material_defaults(struct material *mat);
static void material_set_rgba(struct material *mat, float r, float g, float b, float a);
static enum osh_status material_workspace_init_reserved(struct material_workspace *wm);
static enum osh_status material_set_name(struct material *mat, char const *name);
static enum osh_status material_workspace_complete(struct material_workspace *wm);
static enum osh_status material_workspace_validate_completed(struct material_workspace const *wm);
static enum osh_status material_workspace_validate_raw(struct material_workspace const *wm);
static enum osh_status material_complete(struct material *mat);
static enum osh_status material_complete_composition(struct material *mat);
static enum osh_status material_complete_icru(struct material *mat);
static enum osh_status material_derive_atom_counts_from_mass_fractions(struct material *mat);
static enum osh_status material_derive_mass_fractions_from_atom_counts(struct material *mat);
static void material_default_state_if_unset(struct material *mat);
static enum osh_status material_element_mass(unsigned int z, unsigned int a, double *mass_out);
static enum osh_status material_set_elements_from_icru(struct material *mat,
                                                       struct osh_material_icru_entry const *entry);
static enum osh_status material_complete_mean_excitation_energy(struct material *mat);
static enum osh_status material_fill_element_mee_defaults(struct material *mat);
static enum osh_status material_derive_compound_mee(struct material *mat);
static void material_free_fields(struct material *mat);
static char const *material_state_name(int state);
static int material_has_any_element_mee(struct material const *mat);

enum osh_status
osh_material_setup_from_path(char const *path, struct osh_logger *lg, struct material_workspace **wm_out) {
    enum osh_status rc;
    struct oshfile *sf;
    struct material_workspace *wm;
    char *wdir;

    (void) lg;

    if (!path || !wm_out) {
        return OSH_EINVAL;
    }
    *wm_out = NULL;

    wdir = osh_path_dirname(path);
    sf = osh_fopen(path);
    if (!sf) {
        free(wdir);
        return OSH_EIO;
    }

    wm = (struct material_workspace *) calloc(1, sizeof(*wm));
    if (!wm) {
        free(wdir);
        osh_fclose(sf);
        return OSH_ENOMEM;
    }

    wm->wdir = wdir;
    wdir = NULL;
    wm->fname = strdup(path);
    if (!wm->fname) {
        osh_material_workspace_free(wm);
        osh_fclose(sf);
        return OSH_ENOMEM;
    }

    rc = material_workspace_init_reserved(wm);
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        osh_fclose(sf);
        return rc;
    }

    rc = osh_material_parse(sf, wm);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    rc = material_workspace_validate_raw(wm);
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    rc = material_workspace_complete(wm);
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    rc = material_workspace_validate_completed(wm);
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    if (osh_log_get_level() <= OSH_LOG_INFO) {
        osh_material_print(wm);
    }

    *wm_out = wm;
    return OSH_OK;
}

enum osh_status osh_material_workspace_free(struct material_workspace *wm) {
    size_t i;

    if (!wm) {
        return OSH_OK;
    }

    i = 0;
    while (i < wm->nmaterials) {
        material_free_fields(&wm->materials[i]);
        i++;
    }

    free(wm->materials);
    free(wm->wdir);
    free(wm->fname);
    free(wm);

    return OSH_OK;
}

struct material const *osh_material_by_index(struct material_workspace const *wm, size_t index) {
    if (!wm || index >= wm->nmaterials) {
        return NULL;
    }

    return &wm->materials[index];
}

struct material const *osh_material_by_name(struct material_workspace const *wm, char const *name) {
    size_t i;

    if (!wm || !name) {
        return NULL;
    }

    i = 0;
    while (i < wm->nmaterials) {
        if (wm->materials[i].name && strcmp(wm->materials[i].name, name) == 0) {
            return &wm->materials[i];
        }
        i++;
    }

    return NULL;
}

void osh_material_print(struct material_workspace const *wm) {
    size_t i;
    size_t j;
    struct material const *mat;
    struct material_element const *elem;

    if (!wm) {
        return;
    }

    osh_info("%s", "");
    osh_info("Material configuration:");
    osh_info(OSH_LOG_HLINE);
    osh_info("%-24s : %zu", "Number of materials", wm->nmaterials);

    i = 0;
    while (i < wm->nmaterials) {
        mat = &wm->materials[i];
        osh_info("%s", "");
        osh_info("Material[%zu]: index=%zu name=%s", i, mat->index, mat->name ? mat->name : "(unset)");
        if (mat->icru_id == 0) {
            osh_info("%-24s : %s", "ICRU id", "(unset)");
        } else {
            osh_info("%-24s : %i", "ICRU id", mat->icru_id);
        }
        if (mat->rho < 0.0) {
            osh_info("%-24s : %s", "Density RHO", "(N/A)");
        } else {
            osh_info("%-24s : %.6f g/cm^3", "Density RHO", mat->rho);
        }
        if (mat->mean_excitation_energy < 0.0) {
            osh_info("%-24s : %s", "Mean excitation energy", "(N/A)");
        } else {
            osh_info("%-24s : %.4f eV", "Mean excitation energy", mat->mean_excitation_energy);
        }
        osh_info("%-24s : %s", "dE/dx table", mat->dedx_table_path ? mat->dedx_table_path : "(unset)");
        osh_info("%-24s : %.3g %.3g %.3g %.3g", "RGBA", mat->rgba[0], mat->rgba[1], mat->rgba[2], mat->rgba[3]);
        osh_info("%-24s : %s", "State", material_state_name(mat->state));
        osh_info("%-24s : %zu", "Elements", mat->nelements);

        j = 0;
        while (j < mat->nelements) {
            elem = &mat->elements[j];
            if (elem->mean_excitation_energy < 0.0) {
                osh_info("    element[%zu]: Z=%u A=%u atoms=%.8g mass_fraction=%.8g mean_excitation_energy=(N/A)",
                         j,
                         elem->z,
                         elem->a,
                         elem->atom_count,
                         elem->mass_fraction);
            } else {
                osh_info("    element[%zu]: Z=%u A=%u atoms=%.8g mass_fraction=%.8g mean_excitation_energy=%.4f eV",
                         j,
                         elem->z,
                         elem->a,
                         elem->atom_count,
                         elem->mass_fraction,
                         elem->mean_excitation_energy);
            }
            j++;
        }

        i++;
    }
}

static void material_defaults(struct material *mat) {
    mat->elements = NULL;
    mat->name = NULL;
    mat->dedx_table_path = NULL;
    mat->rho = -1.0;
    mat->mean_excitation_energy = -1.0;
    material_set_rgba(mat, 0.8f, 0.8f, 0.8f, 1.0f);
    mat->nelements = 0u;
    mat->lineno = 0u;
    mat->index = 0u;
    mat->icru_id = 0;
    mat->state = OSH_MATERIAL_STATE_UNSET;
}

static void material_set_rgba(struct material *mat, float r, float g, float b, float a) {
    mat->rgba[0] = r;
    mat->rgba[1] = g;
    mat->rgba[2] = b;
    mat->rgba[3] = a;
}

/**
 * @brief Allocate and populate the two reserved material slots (blackhole, vacuum).
 *
 * @details
 * Reserved slots occupy indices 0 and 1 so that user materials always start at
 * OSH_MATERIAL_INDEX_FIRST_USER. Blackhole is opaque black (alpha=1) with no
 * physical density; vacuum is fully transparent (alpha=0) with rho=0.
 *
 * @param[in,out] wm  Workspace to initialize.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status material_workspace_init_reserved(struct material_workspace *wm) {
    struct material *materials;
    enum osh_status rc;

    materials = (struct material *) calloc(OSH_MATERIAL_INDEX_FIRST_USER, sizeof(*materials));
    if (!materials) {
        return OSH_ENOMEM;
    }

    wm->materials = materials;
    wm->nmaterials = OSH_MATERIAL_INDEX_FIRST_USER;

    material_defaults(&wm->materials[OSH_MATERIAL_INDEX_BLACKHOLE]);
    wm->materials[OSH_MATERIAL_INDEX_BLACKHOLE].index = OSH_MATERIAL_INDEX_BLACKHOLE;
    material_set_rgba(&wm->materials[OSH_MATERIAL_INDEX_BLACKHOLE], 0.0f, 0.0f, 0.0f, 1.0f);
    rc = material_set_name(&wm->materials[OSH_MATERIAL_INDEX_BLACKHOLE], OSH_MATERIAL_NAME_BLACKHOLE);
    if (rc != OSH_OK) {
        return rc;
    }

    material_defaults(&wm->materials[OSH_MATERIAL_INDEX_VACUUM]);
    wm->materials[OSH_MATERIAL_INDEX_VACUUM].index = OSH_MATERIAL_INDEX_VACUUM;
    wm->materials[OSH_MATERIAL_INDEX_VACUUM].rho = 0.0;
    material_set_rgba(&wm->materials[OSH_MATERIAL_INDEX_VACUUM], 0.0f, 0.0f, 0.0f, 0.0f);
    rc = material_set_name(&wm->materials[OSH_MATERIAL_INDEX_VACUUM], OSH_MATERIAL_NAME_VACUUM);
    if (rc != OSH_OK) {
        return rc;
    }

    return OSH_OK;
}

/**
 * @brief Copy @p name into owned heap storage and assign it to @p mat.
 *
 * @param[in,out] mat   Material to update.
 * @param[in]     name  NUL-terminated name string to copy.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status material_set_name(struct material *mat, char const *name) {
    char *copy;
    size_t len;

    len = strlen(name);
    copy = (char *) malloc(len + 1u);
    if (!copy) {
        return OSH_ENOMEM;
    }
    memcpy(copy, name, len + 1u);

    free(mat->name);
    mat->name = copy;

    return OSH_OK;
}

/**
 * @brief Complete parsed material definitions before transport setup.
 *
 * @details
 * This assembly pass expands ICRU convenience definitions and derives the
 * complementary atom-count or mass-fraction fields. Explicit scalar cards such
 * as RHO, STATE, and MATERIALI remain user overrides and are not overwritten.
 */
static enum osh_status material_workspace_complete(struct material_workspace *wm) {
    enum osh_status rc;
    size_t i;

    if (!wm) {
        return OSH_EINVAL;
    }

    i = OSH_MATERIAL_INDEX_FIRST_USER;
    while (i < wm->nmaterials) {
        rc = material_complete(&wm->materials[i]);
        if (rc != OSH_OK) {
            return rc;
        }
        i++;
    }

    return OSH_OK;
}

/**
 * @brief Validate raw parser output before ICRU expansion.
 *
 * @details
 * Raw validation catches user inputs that would become ambiguous after
 * assembly, especially mixing ICRU with explicit element cards in one material.
 */
static enum osh_status material_workspace_validate_raw(struct material_workspace const *wm) {
    size_t i;
    size_t j;
    struct material const *mat;
    struct material_element const *elem;

    if (!wm) {
        return OSH_EINVAL;
    }
    if (wm->nmaterials == 0u) {
        osh_error("material: no materials found");
        return OSH_EPARSE;
    }

    i = 0;
    while (i < wm->nmaterials) {
        mat = &wm->materials[i];
        if (mat->index != i) {
            osh_error("material: material at array slot %zu has inconsistent index %zu", i, mat->index);
            return OSH_ESTATE;
        }
        if (!mat->name) {
            osh_error("material: material index %zu has no name", mat->index);
            return OSH_EPARSE;
        }
        if (mat->index < OSH_MATERIAL_INDEX_FIRST_USER) {
            i++;
            continue;
        }
        if (mat->icru_id == 0 && mat->nelements == 0u) {
            osh_error("material: material '%s' defines neither ICRU nor elemental composition", mat->name);
            return OSH_EPARSE;
        }
        if (mat->icru_id != 0 && mat->nelements > 0u) {
            osh_error("material: material '%s' mixes ICRU and explicit elemental composition", mat->name);
            return OSH_EPARSE;
        }
        if (mat->icru_id == 0 && mat->nelements > 0u && mat->rho < 0.0) {
            osh_error("material: material '%s' requires RHO for elemental composition", mat->name);
            return OSH_EPARSE;
        }
        if (mat->nelements > 1u && mat->mean_excitation_energy >= 0.0 && material_has_any_element_mee(mat)) {
            osh_error("material: compound material '%s' cannot define both material and element mean excitation energy",
                      mat->name);
            return OSH_EPARSE;
        }
        if (mat->rho == 0.0) {
            osh_error("material: material '%s' has invalid RHO %.8g", mat->name, mat->rho);
            return OSH_EPARSE;
        }
        if (mat->icru_id < 0) {
            osh_error("material: material '%s' has invalid ICRU id %i", mat->name, mat->icru_id);
            return OSH_EPARSE;
        }
        if (mat->state != OSH_MATERIAL_STATE_UNSET && mat->state != OSH_MATERIAL_STATE_CONDENSED
            && mat->state != OSH_MATERIAL_STATE_GAS) {
            osh_error("material: material '%s' has invalid STATE %i", mat->name, mat->state);
            return OSH_EPARSE;
        }

        j = 0;
        while (j < mat->nelements) {
            elem = &mat->elements[j];
            if (elem->z == 0u) {
                osh_error("material: material '%s' element %zu has invalid Z=0", mat->name, j);
                return OSH_EPARSE;
            }
            if (elem->atom_count <= 0.0 && elem->mass_fraction <= 0.0) {
                osh_error("material: material '%s' element %zu has neither atom count nor mass fraction", mat->name, j);
                return OSH_EPARSE;
            }
            j++;
        }

        i++;
    }

    return OSH_OK;
}

/**
 * @brief Validate assembled material definitions.
 */
static enum osh_status material_workspace_validate_completed(struct material_workspace const *wm) {
    size_t i;
    size_t j;
    struct material const *mat;
    struct material_element const *elem;
    double mass_fraction_sum; /* dimensionless, must sum to 1 */

    if (!wm) {
        return OSH_EINVAL;
    }

    i = 0;
    while (i < wm->nmaterials) {
        mat = &wm->materials[i];
        if (mat->index != i) {
            osh_error("material: material at array slot %zu has inconsistent index %zu", i, mat->index);
            return OSH_ESTATE;
        }
        if (!mat->name) {
            osh_error("material: material index %zu has no name", mat->index);
            return OSH_EPARSE;
        }
        if (mat->index < OSH_MATERIAL_INDEX_FIRST_USER) {
            i++;
            continue;
        }
        if (mat->rho <= 0.0) {
            osh_error("material: material '%s' has invalid completed RHO %.8g", mat->name, mat->rho);
            return OSH_EPARSE;
        }
        if (mat->nelements == 0u) {
            osh_error("material: material '%s' has no completed elemental composition", mat->name);
            return OSH_EPARSE;
        }
        if (mat->icru_id < 0) {
            osh_error("material: material '%s' has invalid ICRU id %i", mat->name, mat->icru_id);
            return OSH_EPARSE;
        }
        if (mat->state != OSH_MATERIAL_STATE_UNSET && mat->state != OSH_MATERIAL_STATE_CONDENSED
            && mat->state != OSH_MATERIAL_STATE_GAS) {
            osh_error("material: material '%s' has invalid STATE %i", mat->name, mat->state);
            return OSH_EPARSE;
        }

        j = 0;
        mass_fraction_sum = 0.0;
        while (j < mat->nelements) {
            elem = &mat->elements[j];
            if (elem->z == 0u) {
                osh_error("material: material '%s' element %zu has invalid Z=0", mat->name, j);
                return OSH_EPARSE;
            }
            if (elem->atom_count <= 0.0 || elem->mass_fraction <= 0.0) {
                osh_error("material: material '%s' element %zu has incomplete assembled composition", mat->name, j);
                return OSH_EPARSE;
            }
            mass_fraction_sum += elem->mass_fraction;
            j++;
        }
        if (mass_fraction_sum < 0.999999 || mass_fraction_sum > 1.000001) {
            osh_error("material: material '%s' mass fractions sum to %.8g", mat->name, mass_fraction_sum);
            return OSH_EPARSE;
        }

        i++;
    }

    return OSH_OK;
}

/**
 * @brief Complete one user material in place.
 *
 * @details
 * Runs four sub-steps in order:
 *  1. ICRU expansion fills rho, state, MEE and element list from the ICRU
 *     database when an ICRU card is present.
 *  2. Composition derivation computes the missing complementary field
 *     (atom counts from mass fractions or vice versa).
 *  3. State default makes materials condensed unless set explicitly or
 *     by ICRU.
 *  4. Mean-excitation completion preserves fixed material MEE values or, when
 *     material MEE is unset, fills element defaults and derives the material
 *     value with Bragg additivity.
 *
 * @param[in,out] mat  Material to complete.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
static enum osh_status material_complete(struct material *mat) {
    enum osh_status rc;

    if (!mat) {
        return OSH_EINVAL;
    }

    rc = material_complete_icru(mat);
    if (rc != OSH_OK) {
        return rc;
    }

    rc = material_complete_composition(mat);
    if (rc != OSH_OK) {
        return rc;
    }

    material_default_state_if_unset(mat);

    return material_complete_mean_excitation_energy(mat);
}

/**
 * @brief Default user-defined materials to condensed matter.
 */
static void material_default_state_if_unset(struct material *mat) {
    if (mat->state == OSH_MATERIAL_STATE_UNSET) {
        mat->state = OSH_MATERIAL_STATE_CONDENSED;
    }
}

/**
 * @brief Expand an ICRU material while preserving explicit scalar overrides.
 */
static enum osh_status material_complete_icru(struct material *mat) {
    enum osh_status rc;
    struct osh_material_icru_entry entry;

    if (mat->icru_id == 0) {
        return OSH_OK;
    }

    rc = osh_material_icru_lookup(mat->icru_id, &entry);
    if (rc != OSH_OK) {
        osh_error("material: material '%s' references unavailable ICRU id %i", mat->name, mat->icru_id);
        return rc;
    }

    if (mat->rho < 0.0) {
        mat->rho = entry.rho;
    }
    if (mat->state == OSH_MATERIAL_STATE_UNSET) {
        mat->state = entry.state;
    }
    if (mat->mean_excitation_energy < 0.0) {
        mat->mean_excitation_energy = entry.mean_excitation_energy;
    }

    return material_set_elements_from_icru(mat, &entry);
}

/**
 * @brief Fill missing complementary composition fields.
 */
static enum osh_status material_complete_composition(struct material *mat) {
    if (mat->nelements == 0u) {
        return OSH_OK;
    }

    if (mat->elements[0].mass_fraction > 0.0 && mat->elements[0].atom_count < 0.0) {
        return material_derive_atom_counts_from_mass_fractions(mat);
    }
    if (mat->elements[0].atom_count > 0.0 && mat->elements[0].mass_fraction < 0.0) {
        return material_derive_mass_fractions_from_atom_counts(mat);
    }

    return OSH_OK;
}

/**
 * @brief Replace material element composition from an ICRU entry.
 */
static enum osh_status material_set_elements_from_icru(struct material *mat,
                                                       struct osh_material_icru_entry const *entry) {
    struct material_element *elements;
    size_t i;

    if (mat->nelements != 0u) {
        return OSH_ESTATE;
    }

    elements = (struct material_element *) calloc(entry->nelements, sizeof(*elements));
    if (!elements) {
        return OSH_ENOMEM;
    }

    i = 0;
    while (i < entry->nelements) {
        elements[i].atom_count = -1.0;
        elements[i].mass_fraction = entry->elements[i].mass_fraction;
        elements[i].mean_excitation_energy = -1.0;
        elements[i].lineno = mat->lineno;
        elements[i].z = entry->elements[i].z;
        elements[i].a = entry->elements[i].a;
        i++;
    }

    free(mat->elements);
    mat->elements = elements;
    mat->nelements = entry->nelements;

    return OSH_OK;
}

/**
 * @brief Complete material and element mean excitation energies.
 *
 * @details
 * Material-level MEE is authoritative for compounds. If it is already known
 * from ICRU or a user MATERIALI/MIVALUE/MIAV card, compound element MEE values
 * are left unset because the inverse Bragg problem is underdetermined. Raw
 * validation rejects compound input that tries to define both material- and
 * element-level MEE explicitly. For a single-element material, the element can
 * safely inherit the material value.
 *
 * If material-level MEE is unset, missing element MEE values are filled from
 * ICRU defaults with ICRU-49 in-compound corrections where available, then the
 * material value is derived using Bragg additivity.
 *
 * @param[in,out] mat  Material to complete.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
static enum osh_status material_complete_mean_excitation_energy(struct material *mat) {
    enum osh_status rc;
    size_t i;

    if (mat->nelements == 0u) {
        return OSH_OK;
    }

    if (mat->mean_excitation_energy >= 0.0) {
        if (mat->nelements == 1u) {
            mat->elements[0].mean_excitation_energy = mat->mean_excitation_energy;
            return OSH_OK;
        }

        i = 0;
        while (i < mat->nelements) {
            mat->elements[i].mean_excitation_energy = -1.0;
            i++;
        }
        return OSH_OK;
    }

    rc = material_fill_element_mee_defaults(mat);
    if (rc != OSH_OK) {
        return rc;
    }

    return material_derive_compound_mee(mat);
}

/**
 * @brief Fill per-element mean excitation energies from ICRU data.
 *
 * @details
 * For each element whose mean_excitation_energy is unset (< 0), looks up the
 * ICRU elemental entry (ICRU id == Z for pure elements). For compounds, ICRU 49
 * Table 2.11 gives adjusted element mean excitation energies for selected
 * elements in gaseous or condensed compounds; those are used when available.
 * This matters for common compounds such as water, where pure oxygen (95 eV)
 * would otherwise under-estimate the compound-derived material MEE.
 *
 * Elements with an explicit user IVALUE override are left unchanged. If the
 * lookup fails, the element MEE remains -1; the subsequent Bragg-additivity
 * step will then skip compound MEE derivation, and the user or a dE/dx table
 * must supply it.
 *
 * @param[in,out] mat  Material whose element MEE fields are to be filled.
 *
 * @returns OSH_OK on success.
 */
static enum osh_status material_fill_element_mee_defaults(struct material *mat) {
    struct osh_material_icru_entry entry;
    size_t i;

    i = 0;
    while (i < mat->nelements) {
        if (mat->elements[i].mean_excitation_energy < 0.0) {
            if (osh_material_icru_lookup((int) mat->elements[i].z, &entry) == OSH_OK) {
                if (mat->nelements > 1u) {
                    mat->elements[i].mean_excitation_energy = osh_material_icru_compound_element_mean_excitation_energy(
                        mat->elements[i].z, mat->state, entry.mean_excitation_energy);
                } else {
                    mat->elements[i].mean_excitation_energy = entry.mean_excitation_energy;
                }
            }
        }
        i++;
    }

    return OSH_OK;
}

/**
 * @brief Derive the material mean excitation energy via Bragg additivity.
 *
 * @details
 * When the material-level mean excitation energy is unset (< 0) and all
 * elements have known MEE values, this function applies the Bragg-additivity
 * rule from ICRU Report 37:
 *
 *   ln(I) = sum_i(w_i * (Z_i/A_i) * ln(I_i)) / sum_i(w_i * Z_i/A_i)
 *
 * where w_i is the mass fraction, Z_i the atomic number, A_i the atomic mass
 * [Da] from either the natural element or explicit isotope, and I_i the mean
 * excitation energy of element i [eV].
 *
 * If any element still has an unset MEE, the derivation is skipped and the
 * material MEE remains -1 (to be supplied by the user or a dE/dx table).
 * A user-set or ICRU-set material-level MEE is never overwritten.
 *
 * @param[in,out] mat  Material whose mean_excitation_energy is to be derived.
 *
 * @returns OSH_OK on success, or an error code from the atomic-mass lookup.
 */
static enum osh_status material_derive_compound_mee(struct material *mat) {
    enum osh_status rc;
    double numerator;
    double denominator;
    double mass;       /* atomic mass of natural element or explicit isotope [Da] */
    double zi_over_ai; /* Z_i / A_i  [1/Da] */
    size_t i;

    if (mat->mean_excitation_energy >= 0.0 || mat->nelements == 0u) {
        return OSH_OK;
    }

    numerator = 0.0;
    denominator = 0.0;
    i = 0;
    while (i < mat->nelements) {
        if (mat->elements[i].mean_excitation_energy < 0.0 || mat->elements[i].mass_fraction < 0.0) {
            return OSH_OK; /* incomplete data; skip derivation */
        }
        rc = material_element_mass(mat->elements[i].z, mat->elements[i].a, &mass);
        if (rc != OSH_OK) {
            return rc;
        }
        zi_over_ai = (double) mat->elements[i].z / mass;
        numerator += mat->elements[i].mass_fraction * zi_over_ai * log(mat->elements[i].mean_excitation_energy);
        denominator += mat->elements[i].mass_fraction * zi_over_ai;
        i++;
    }

    if (denominator > 0.0) {
        mat->mean_excitation_energy = exp(numerator / denominator);
    }

    return OSH_OK;
}

/**
 * @brief Derive relative atom counts from mass fractions.
 *
 * @details
 * For each element i with mass fraction w_i and atomic mass M_i [Da], the
 * normalized atom count is:
 *
 *   n_i = (w_i / M_i) / sum_j(w_j / M_j)
 *
 * The resulting n_i sum to 1. They are relative mole fractions; exact
 * stoichiometric integer ratios are not required for transport.
 *
 * @param[in,out] mat  Material whose atom_count fields are to be filled.
 *
 * @returns OSH_OK on success, OSH_EINVAL/OSH_ESTATE on invalid element data.
 */
static enum osh_status material_derive_atom_counts_from_mass_fractions(struct material *mat) {
    enum osh_status rc;
    double mass;      /* atomic mass of current element [Da] */
    double sum_moles; /* sum of mass_fraction/mass across all elements [1/Da] */
    size_t i;

    sum_moles = 0.0;
    i = 0;
    while (i < mat->nelements) {
        rc = material_element_mass(mat->elements[i].z, mat->elements[i].a, &mass);
        if (rc != OSH_OK) {
            osh_error("material: material '%s' element %zu uses unsupported Z=%u A=%u",
                      mat->name,
                      i,
                      mat->elements[i].z,
                      mat->elements[i].a);
            return rc;
        }
        sum_moles += mat->elements[i].mass_fraction / mass;
        i++;
    }

    if (sum_moles <= 0.0) {
        return OSH_ESTATE;
    }

    i = 0;
    while (i < mat->nelements) {
        rc = material_element_mass(mat->elements[i].z, mat->elements[i].a, &mass);
        if (rc != OSH_OK) {
            return rc;
        }
        mat->elements[i].atom_count = (mat->elements[i].mass_fraction / mass) / sum_moles;
        i++;
    }

    return OSH_OK;
}

/**
 * @brief Derive normalized mass fractions from relative atom counts.
 *
 * @details
 * For each element i with relative atom count n_i and atomic mass M_i [Da],
 * the mass fraction is:
 *
 *   w_i = (n_i * M_i) / sum_j(n_j * M_j)
 *
 * The resulting w_i sum to 1.
 *
 * @param[in,out] mat  Material whose mass_fraction fields are to be filled.
 *
 * @returns OSH_OK on success, OSH_EINVAL/OSH_ESTATE on invalid element data.
 */
static enum osh_status material_derive_mass_fractions_from_atom_counts(struct material *mat) {
    enum osh_status rc;
    double mass;       /* atomic mass of current element [Da] */
    double total_mass; /* sum of atom_count * mass across all elements [Da] */
    size_t i;

    total_mass = 0.0;
    i = 0;
    while (i < mat->nelements) {
        rc = material_element_mass(mat->elements[i].z, mat->elements[i].a, &mass);
        if (rc != OSH_OK) {
            osh_error("material: material '%s' element %zu uses unsupported Z=%u A=%u",
                      mat->name,
                      i,
                      mat->elements[i].z,
                      mat->elements[i].a);
            return rc;
        }
        total_mass += mat->elements[i].atom_count * mass;
        i++;
    }

    if (total_mass <= 0.0) {
        return OSH_ESTATE;
    }

    i = 0;
    while (i < mat->nelements) {
        rc = material_element_mass(mat->elements[i].z, mat->elements[i].a, &mass);
        if (rc != OSH_OK) {
            return rc;
        }
        mat->elements[i].mass_fraction = (mat->elements[i].atom_count * mass) / total_mass;
        i++;
    }

    return OSH_OK;
}

/**
 * @brief Resolve an element mass for natural elements or explicit isotopes.
 */
static enum osh_status material_element_mass(unsigned int z, unsigned int a, double *mass_out) {
    struct isotope iso;

    if (!mass_out || z == 0u) {
        return OSH_EINVAL;
    }

    if (a == 0u) {
        return osh_material_natural_atomic_mass_da(z, mass_out);
    }

    if (!osh_isotope_from_za(&iso, z, a)) {
        return OSH_EINVAL;
    }
    *mass_out = iso.amass;
    return OSH_OK;
}

/**
 * @brief Free all heap-owned fields of @p mat and reset it to defaults.
 *
 * @details Does not free @p mat itself; the caller owns the struct storage.
 *
 * @param[in,out] mat  Material whose fields are to be released.
 */
static void material_free_fields(struct material *mat) {
    if (!mat) {
        return;
    }

    free(mat->elements);
    free(mat->name);
    free(mat->dedx_table_path);
    material_defaults(mat);
}

static char const *material_state_name(int state) {
    switch (state) {
    case OSH_MATERIAL_STATE_UNSET:
        return "unset";
    case OSH_MATERIAL_STATE_CONDENSED:
        return "condensed";
    case OSH_MATERIAL_STATE_GAS:
        return "gas";
    default:
        return "invalid";
    }
}

/**
 * @brief Return 1 when any element has an explicit/raw mean excitation value.
 */
static int material_has_any_element_mee(struct material const *mat) {
    size_t i;

    if (!mat) {
        return 0;
    }

    for (i = 0; i < mat->nelements; ++i) {
        if (mat->elements[i].mean_excitation_energy >= 0.0) {
            return 1;
        }
    }

    return 0;
}
