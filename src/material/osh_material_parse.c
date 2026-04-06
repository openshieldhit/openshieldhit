#include "material/osh_material_parse.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "material/osh_material_parse_keys.h"

static enum osh_status parse_density(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_color(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_element_by_mass(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_element_by_number(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status
parse_element_mean_excitation_energy(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_end(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_loaddedx(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status
parse_material_mean_excitation_energy(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_icru(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_material_start(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_state(struct material_workspace *wm, struct oshfile *oshf, char const *args);

static struct material *material_current(struct material_workspace *wm);
static void material_defaults(struct material *mat);
static enum osh_status parse_mean_excitation_energy_value(double *mean_excitation_energy_out,
                                                          struct oshfile *oshf,
                                                          char const *args,
                                                          char const *key_name);
static enum osh_status material_push(struct material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status material_push_element(struct material *mat,
                                             struct oshfile *oshf,
                                             unsigned int z,
                                             unsigned int a,
                                             double atom_count,
                                             double mass_fraction);
static enum osh_status copy_string(char **dst, char const *src);
static enum osh_status parse_state_value(int *state_out, char const *token);

struct material_dispatch_entry {
    char const *key;
    enum osh_status (*handler)(struct material_workspace *wm, struct oshfile *oshf, char const *args);
};

/*
 * Preferred mean-excitation keys are MATERIALI for material-level fallback
 * data and ELEMENTI for the last parsed element. MIVALUE/MIAV and IVALUE/IAV
 * are kept as short aliases for those two explicit scopes.
 *
 * ICRU and explicit element cards are alternative composition sources. Scalar
 * cards such as RHO, STATE, MATERIALI, and LOADDEDX remain user overrides that
 * later material assembly should preserve when it fills unset defaults.
 */
static struct material_dispatch_entry dispatch_table[] = {
    {OSH_MATERIAL_KEY_COLOR, parse_color},
    {OSH_MATERIAL_KEY_COLOUR, parse_color},
    {OSH_MATERIAL_KEY_DENSITY, parse_density},
    {OSH_MATERIAL_KEY_ELEMENT, parse_element_by_number},
    {OSH_MATERIAL_KEY_ELEMENTBYMASS, parse_element_by_mass},
    {OSH_MATERIAL_KEY_ELEMENTBYNUMBER, parse_element_by_number},
    {OSH_MATERIAL_KEY_ELEMENTI, parse_element_mean_excitation_energy},
    {OSH_MATERIAL_KEY_END, parse_end},
    {OSH_MATERIAL_KEY_IAV, parse_element_mean_excitation_energy},
    {OSH_MATERIAL_KEY_ICRU, parse_icru},
    {OSH_MATERIAL_KEY_IVALUE, parse_element_mean_excitation_energy},
    {OSH_MATERIAL_KEY_LOADDEDX, parse_loaddedx},
    {OSH_MATERIAL_KEY_MATERIAL, parse_material_start},
    {OSH_MATERIAL_KEY_MATERIALI, parse_material_mean_excitation_energy},
    {OSH_MATERIAL_KEY_MIAV, parse_material_mean_excitation_energy},
    {OSH_MATERIAL_KEY_MIVALUE, parse_material_mean_excitation_energy},
    {OSH_MATERIAL_KEY_NUCLID, parse_element_by_number},
    {OSH_MATERIAL_KEY_RHO, parse_density},
    {OSH_MATERIAL_KEY_STATE, parse_state},
    {NULL, NULL}};

enum osh_status osh_material_parse(struct oshfile *oshf, struct material_workspace *wm) {
    enum osh_status rc;
    char *line;
    char *key;
    char *args;
    int lineno;
    int i;
    int found;
    int material_active;

    if (!oshf || !wm) {
        return OSH_EINVAL;
    }

    line = NULL;
    material_active = 0;
    while (osh_readline_key(oshf, &line, &key, &args, &lineno) != -1) {
        i = 0;
        while (key[i] != '\0') {
            key[i] = (char) tolower((unsigned char) key[i]);
            i++;
        }

        if (strcmp(OSH_MATERIAL_KEY_MATERIAL, key) == 0) {
            rc = parse_material_start(wm, oshf, args);
            free(line);
            line = NULL;
            if (rc != OSH_OK) {
                return rc;
            }
            material_active = 1;
            continue;
        }

        if (strcmp(OSH_MATERIAL_KEY_END, key) == 0) {
            if (!material_active) {
                osh_error("in %s line %d: END outside MATERIAL block", oshf->filename, lineno);
                free(line);
                return OSH_EPARSE;
            }
            rc = parse_end(wm, oshf, args);
            free(line);
            line = NULL;
            if (rc != OSH_OK) {
                return rc;
            }
            material_active = 0;
            continue;
        }

        if (!material_active) {
            osh_error("in %s line %d: material key '%s' outside MATERIAL block", oshf->filename, lineno, key);
            free(line);
            return OSH_EPARSE;
        }

        found = 0;
        i = 0;
        while (dispatch_table[i].key != NULL) {
            if (strcmp(dispatch_table[i].key, key) == 0) {
                rc = dispatch_table[i].handler(wm, oshf, args);
                free(line);
                line = NULL;
                if (rc != OSH_OK) {
                    return rc;
                }
                found = 1;
                break;
            }
            i++;
        }

        if (!found) {
            osh_error("in %s line %d: unknown material key '%s'", oshf->filename, lineno, key);
            free(line);
            return OSH_EPARSE;
        }
    }

    if (material_active) {
        osh_error("in %s: MATERIAL block missing END", oshf->filename);
        free(line);
        return OSH_EPARSE;
    }

    return OSH_OK;
}

static enum osh_status parse_density(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    double rho;
    char extra;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: RHO outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%lf %c", &rho, &extra) != 1) {
        osh_error("in %s line %i: RHO expects one floating point value", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (rho <= 0.0) {
        osh_error("in %s line %i: RHO must be > 0", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    mat->rho = rho;
    return OSH_OK;
}

static enum osh_status parse_color(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    float rgba[4];
    char extra;
    int i;

    /* TODO: accept named colors following a common convention, e.g. Python/matplotlib color names. */

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: COLOR outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%f %f %f %f %c", &rgba[0], &rgba[1], &rgba[2], &rgba[3], &extra) != 4) {
        osh_error("in %s line %i: COLOR expects four floating point values", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    i = 0;
    while (i < 4) {
        if (rgba[i] < 0.0f || rgba[i] > 1.0f) {
            osh_error("in %s line %i: COLOR components must be within [0,1]", oshf->filename, oshf->lineno);
            return OSH_EPARSE;
        }
        mat->rgba[i] = rgba[i];
        i++;
    }

    return OSH_OK;
}

static enum osh_status parse_element_by_mass(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    unsigned int z;
    double mass_fraction;
    char extra;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: ELEMENTBYMASS outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%u %lf %c", &z, &mass_fraction, &extra) != 2) {
        osh_error("in %s line %i: ELEMENTBYMASS expects '<Z> <fraction>'", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    return material_push_element(mat, oshf, z, 0u, -1.0, mass_fraction);
}

static enum osh_status parse_element_by_number(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    unsigned int z;
    double atom_count;
    char extra;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: NUCLID outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%u %lf %c", &z, &atom_count, &extra) != 2) {
        osh_error("in %s line %i: NUCLID expects '<Z> <amount>'", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    return material_push_element(mat, oshf, z, 0u, atom_count, -1.0);
}

static enum osh_status parse_end(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    (void) wm;

    if (args && args[0] != '\0') {
        osh_error("in %s line %i: END expects no arguments", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    return OSH_OK;
}

static enum osh_status
parse_element_mean_excitation_energy(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    double mean_excitation_energy;
    enum osh_status rc;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: ELEMENTI/IVALUE/IAV outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (mat->nelements == 0u) {
        osh_error("in %s line %i: ELEMENTI/IVALUE/IAV requires a preceding element card", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    rc = parse_mean_excitation_energy_value(&mean_excitation_energy, oshf, args, "ELEMENTI/IVALUE/IAV");
    if (rc != OSH_OK) {
        return rc;
    }

    mat->elements[mat->nelements - 1u].mean_excitation_energy = mean_excitation_energy;

    return OSH_OK;
}

static enum osh_status parse_loaddedx(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    char path[512];
    char extra;
    char *resolved_path;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: LOADDEDX outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%511s %c", path, &extra) != 1) {
        osh_error("in %s line %i: LOADDEDX expects one path", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    /* Store the table path only. Table loading and projectile coverage are resolved during setup. */
    resolved_path = NULL;
    if (osh_relative_path_to_file(&resolved_path, wm->wdir, path) != 0) {
        return OSH_ENOMEM;
    }

    free(mat->dedx_table_path);
    mat->dedx_table_path = resolved_path;

    return OSH_OK;
}

static enum osh_status
parse_material_mean_excitation_energy(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    double mean_excitation_energy;
    enum osh_status rc;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: MATERIALI outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    rc = parse_mean_excitation_energy_value(&mean_excitation_energy, oshf, args, "MATERIALI/MIVALUE/MIAV");
    if (rc != OSH_OK) {
        return rc;
    }

    /* This value remains the fallback for projectiles not covered by external dE/dx tables. */
    mat->mean_excitation_energy = mean_excitation_energy;

    return OSH_OK;
}

static enum osh_status parse_icru(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *mat;
    int icru_id;
    char extra;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: ICRU outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%i %c", &icru_id, &extra) != 1) {
        osh_error("in %s line %i: ICRU expects one integer value", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (icru_id < 0) {
        osh_error("in %s line %i: ICRU must be >= 0", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    /* ICRU selects predefined material composition/density defaults; mean excitation is configured separately. */
    mat->icru_id = icru_id;
    return OSH_OK;
}

static enum osh_status parse_material_start(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    return material_push(wm, oshf, args);
}

static enum osh_status parse_state(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    enum osh_status rc;
    struct material *mat;
    char state_token[32];
    char extra;
    int state;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: STATE outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%31s %c", state_token, &extra) != 1) {
        osh_error("in %s line %i: STATE expects an integer or keyword", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    rc = parse_state_value(&state, state_token);
    if (rc != OSH_OK) {
        osh_error("in %s line %i: unknown STATE '%s'", oshf->filename, oshf->lineno, state_token);
        return rc;
    }

    mat->state = state;
    return OSH_OK;
}

static struct material *material_current(struct material_workspace *wm) {
    if (!wm || wm->nmaterials == 0u) {
        return NULL;
    }

    return &wm->materials[wm->nmaterials - 1u];
}

static void material_defaults(struct material *mat) {
    mat->elements = NULL;
    mat->name = NULL;
    mat->dedx_table_path = NULL;
    mat->rho = -1.0;
    mat->mean_excitation_energy = -1.0;
    mat->rgba[0] = 0.8f;
    mat->rgba[1] = 0.8f;
    mat->rgba[2] = 0.8f;
    mat->rgba[3] = 1.0f;
    mat->nelements = 0u;
    mat->lineno = 0u;
    mat->index = 0u;
    mat->icru_id = 0;
    mat->state = OSH_MATERIAL_STATE_UNSET;
}

static enum osh_status parse_mean_excitation_energy_value(double *mean_excitation_energy_out,
                                                          struct oshfile *oshf,
                                                          char const *args,
                                                          char const *key_name) {
    double mean_excitation_energy;
    char extra;

    if (!args || sscanf(args, "%lf %c", &mean_excitation_energy, &extra) != 1) {
        osh_error("in %s line %i: %s expects one floating point value", oshf->filename, oshf->lineno, key_name);
        return OSH_EPARSE;
    }
    if (mean_excitation_energy < 0.0) {
        osh_error("in %s line %i: %s must be >= 0", oshf->filename, oshf->lineno, key_name);
        return OSH_EPARSE;
    }

    *mean_excitation_energy_out = mean_excitation_energy;
    return OSH_OK;
}

static enum osh_status material_push(struct material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct material *materials;
    struct material *mat;
    size_t nmaterials;
    char name[128];
    char extra;

    if (!args || sscanf(args, "%127s %c", name, &extra) != 1) {
        osh_error("in %s line %i: MATERIAL expects one name", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    if (osh_material_by_name(wm, name) != NULL) {
        osh_error("in %s line %i: duplicate MATERIAL name '%s'", oshf->filename, oshf->lineno, name);
        return OSH_EPARSE;
    }

    nmaterials = wm->nmaterials + 1u;
    materials = (struct material *) realloc(wm->materials, nmaterials * sizeof(*materials));
    if (!materials) {
        return OSH_ENOMEM;
    }

    wm->materials = materials;
    wm->nmaterials = nmaterials;

    mat = &wm->materials[nmaterials - 1u];
    material_defaults(mat);
    mat->index = nmaterials - 1u;
    mat->lineno = (size_t) oshf->lineno;

    return copy_string(&mat->name, name);
}

static enum osh_status material_push_element(struct material *mat,
                                             struct oshfile *oshf,
                                             unsigned int z,
                                             unsigned int a,
                                             double atom_count,
                                             double mass_fraction) {
    struct material_element *elements;
    struct material_element *elem;
    size_t nelements;
    size_t i;

    if (z == 0u) {
        osh_error("in %s line %i: element Z must be > 0", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (atom_count <= 0.0 && mass_fraction <= 0.0) {
        osh_error("in %s line %i: element amount must be > 0", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    i = 0;
    while (i < mat->nelements) {
        if ((atom_count > 0.0 && mat->elements[i].mass_fraction > 0.0)
            || (mass_fraction > 0.0 && mat->elements[i].atom_count > 0.0)) {
            osh_error(
                "in %s line %i: mixed material element input modes are not supported", oshf->filename, oshf->lineno);
            return OSH_EPARSE;
        }
        i++;
    }

    nelements = mat->nelements + 1u;
    elements = (struct material_element *) realloc(mat->elements, nelements * sizeof(*elements));
    if (!elements) {
        return OSH_ENOMEM;
    }
    mat->elements = elements;

    elem = &mat->elements[mat->nelements];
    elem->atom_count = atom_count;
    elem->mass_fraction = mass_fraction;
    elem->mean_excitation_energy = -1.0;
    elem->lineno = (size_t) oshf->lineno;
    elem->z = z;
    elem->a = a;

    mat->nelements = nelements;

    return OSH_OK;
}

static enum osh_status copy_string(char **dst, char const *src) {
    char *copy;
    size_t len;

    len = strlen(src);
    copy = (char *) malloc(len + 1u);
    if (!copy) {
        return OSH_ENOMEM;
    }
    memcpy(copy, src, len + 1u);

    free(*dst);
    *dst = copy;

    return OSH_OK;
}

static enum osh_status parse_state_value(int *state_out, char const *token) {
    int state;
    char extra;

    if (sscanf(token, "%i%c", &state, &extra) == 1) {
        if (state == OSH_MATERIAL_STATE_CONDENSED || state == OSH_MATERIAL_STATE_GAS) {
            *state_out = state;
            return OSH_OK;
        }
        return OSH_EPARSE;
    }

    if (strcasecmp(token, "condensed") == 0 || strcasecmp(token, "solid") == 0 || strcasecmp(token, "liquid") == 0) {
        *state_out = OSH_MATERIAL_STATE_CONDENSED;
        return OSH_OK;
    }
    if (strcasecmp(token, "gas") == 0) {
        *state_out = OSH_MATERIAL_STATE_GAS;
        return OSH_OK;
    }

    return OSH_EPARSE;
}
