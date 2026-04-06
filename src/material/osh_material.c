#include "material/osh_material.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "material/osh_material_parse.h"

static void material_defaults(struct material *mat);
static void material_set_rgba(struct material *mat, float r, float g, float b, float a);
static enum osh_status material_workspace_init_reserved(struct material_workspace *wm);
static enum osh_status material_set_name(struct material *mat, char const *name);
static enum osh_status material_workspace_validate(struct material_workspace const *wm);
static void material_free_fields(struct material *mat);
static char const *material_state_name(int state);

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

    rc = material_workspace_validate(wm);
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    if (osh_log_get_level() <= OSH_LOG_DEBUG) {
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
        osh_info("%-24s : %i", "ICRU id", mat->icru_id);
        osh_info("%-24s : %.8g g/cm^3", "Density RHO", mat->rho);
        osh_info("%-24s : %.8g eV", "Mean excitation energy", mat->mean_excitation_energy);
        osh_info("%-24s : %.3g %.3g %.3g %.3g", "RGBA", mat->rgba[0], mat->rgba[1], mat->rgba[2], mat->rgba[3]);
        osh_info("%-24s : %s", "State", material_state_name(mat->state));
        osh_info("%-24s : %zu", "Elements", mat->nelements);

        j = 0;
        while (j < mat->nelements) {
            elem = &mat->elements[j];
            osh_info("    element[%zu]: Z=%u A=%u atoms=%.8g mass_fraction=%.8g mean_excitation_energy=%.8g eV",
                     j,
                     elem->z,
                     elem->a,
                     elem->atom_count,
                     elem->mass_fraction,
                     elem->mean_excitation_energy);
            j++;
        }

        i++;
    }
}

static void material_defaults(struct material *mat) {
    mat->elements = NULL;
    mat->name = NULL;
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

static enum osh_status material_workspace_validate(struct material_workspace const *wm) {
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
        if (mat->rho == 0.0) {
            osh_error("material: material '%s' has invalid RHO %.8g", mat->name, mat->rho);
            return OSH_EPARSE;
        }
        if (mat->icru_id == 0 && mat->nelements == 0u) {
            osh_error("material: material '%s' defines neither ICRU nor elemental composition", mat->name);
            return OSH_EPARSE;
        }
        if (mat->icru_id == 0 && mat->nelements > 0u && mat->rho < 0.0) {
            osh_error("material: material '%s' requires RHO for elemental composition", mat->name);
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

static void material_free_fields(struct material *mat) {
    if (!mat) {
        return;
    }

    free(mat->elements);
    free(mat->name);
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
