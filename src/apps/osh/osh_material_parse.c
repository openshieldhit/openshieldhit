#include "apps/osh/osh_material_parse.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_material_loaddedx.h"
#include "apps/osh/osh_material_parse_keys.h"
#include "common/osh_file.h"
#include "common/osh_readline.h"
#include "openshieldhit/logger.h"

static enum osh_status parse_density(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_color(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_element_by_mass(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status
parse_element_by_number(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status
parse_element_mean_excitation_energy(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_end(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_loaddedx(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status
parse_material_mean_excitation_energy(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_icru(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_material_start(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status parse_state(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);

static struct osh_material *material_current(struct osh_material_workspace *wm);
static void material_defaults(struct osh_material *mat);
static enum osh_status parse_mean_excitation_energy_value(double *mean_excitation_energy_out,
                                                          struct oshfile *oshf,
                                                          char const *args,
                                                          char const *key_name);
static enum osh_status parse_element_card_args(unsigned int *z_out,
                                               unsigned int *a_out,
                                               double *amount_out,
                                               struct oshfile *oshf,
                                               char const *args,
                                               char const *key_name);
static enum osh_status parse_double_token(double *value_out, char const *token);
static enum osh_status parse_uint_token(unsigned int *value_out, char const *token);
static enum osh_status material_push(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
static enum osh_status material_push_element(struct osh_material *mat,
                                             struct oshfile *oshf,
                                             unsigned int z,
                                             unsigned int a,
                                             double atom_count,
                                             double mass_fraction);
static enum osh_status copy_string(char **dst, char const *src);
static enum osh_status parse_state_value(int *state_out, char const *token);

struct material_dispatch_entry {
    char const *key;
    enum osh_status (*handler)(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args);
};

/*
 * Preferred mean-excitation keys are MATERIALI for material-level fallback
 * data and ELEMENTI for the last parsed element. MIVALUE/MIAV and IVALUE/IAV
 * are kept as short aliases for those two explicit scopes.
 *
 * ICRU and explicit element cards are alternative composition sources. Scalar
 * cards such as RHO, STATE, MATERIALI, and LOADDEDX remain user overrides that
 * later material assembly should preserve when it fills unset defaults.
 *
 * Element cards use `Z amount` for natural elements (A=0 internally) and
 * `Z A amount` for explicit isotopes.
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

enum osh_status osh_material_parse(struct oshfile *oshf, struct osh_material_workspace *wm) {
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

    return OSH_OK;
}

/**
 * @brief Handle the RHO / DENSITY card.
 *
 * @param[in,out] wm    Material workspace; the current material receives the density.
 * @param[in]     oshf  Open input file (used for error location).
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static enum osh_status parse_density(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
    double rho; /* [g/cm³] */
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

/**
 * @brief Handle the COLOR / COLOUR card.
 *
 * @details Expects four floating-point values in [0,1]: red, green, blue, alpha.
 *
 * @param[in,out] wm    Material workspace.
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static enum osh_status parse_color(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
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

/**
 * @brief Handle the ELEMENTBYMASS card.
 *
 * @details Appends an element specified by Z (and optionally A) with a mass
 * fraction. Atom count is derived during post-parse completion.
 *
 * @param[in,out] wm    Material workspace.
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static enum osh_status
parse_element_by_mass(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
    unsigned int z;
    unsigned int a;
    double mass_fraction;
    enum osh_status rc;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: ELEMENTBYMASS outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    rc = parse_element_card_args(&z, &a, &mass_fraction, oshf, args, "ELEMENTBYMASS");
    if (rc != OSH_OK) {
        return rc;
    }

    return material_push_element(mat, oshf, z, a, -1.0, mass_fraction);
}

/**
 * @brief Handle the NUCLID / ELEMENT / ELEMENTBYNUMBER card.
 *
 * @details Appends an element specified by Z (and optionally A) with a relative
 * atom count. Mass fraction is derived during post-parse completion.
 *
 * @param[in,out] wm    Material workspace.
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static enum osh_status
parse_element_by_number(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
    unsigned int z;
    unsigned int a;
    double atom_count;
    enum osh_status rc;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: NUCLID outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    rc = parse_element_card_args(&z, &a, &atom_count, oshf, args, "NUCLID/ELEMENT");
    if (rc != OSH_OK) {
        return rc;
    }

    return material_push_element(mat, oshf, z, a, atom_count, -1.0);
}

/**
 * @brief Handle the END card closing a MATERIAL block.
 *
 * @details
 * END is optional. A new MATERIAL card or EOF also terminates the previous
 * material block. END remains accepted as an explicit delimiter for backward
 * compatibility and readability. Outside a material block it acts as a no-op.
 *
 * @param[in,out] wm    Material workspace (unused; present for dispatch signature).
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line; must be empty.
 *
 * @returns OSH_OK on success, OSH_EPARSE if trailing arguments are present.
 */
static enum osh_status parse_end(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    (void) wm;

    if (args && args[0] != '\0') {
        osh_error("in %s line %i: END expects no arguments", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    return OSH_OK;
}

/**
 * @brief Handle the ELEMENTI / IVALUE / IAV card.
 *
 * @details Sets the mean excitation energy [eV] on the most recently parsed
 * element. Must follow a NUCLID/ELEMENT/ELEMENTBYMASS card.
 *
 * @param[in,out] wm    Material workspace.
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input or ordering error.
 */
static enum osh_status
parse_element_mean_excitation_energy(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
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

/**
 * @brief Handle the LOADDEDX card.
 *
 * @details
 * Relative paths are resolved against the directory of the material input
 * file. The legacy LOADDEDX text table is read immediately by the app parser
 * and translated into material-owned dE/dx overrides in the public cold model.
 * Repeated LOADDEDX cards replace the earlier override set for the current
 * material, matching the old "one path per material" behavior.
 *
 * @param[in,out] wm    Material workspace (wdir used for path resolution).
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE or OSH_ENOMEM on failure.
 */
static enum osh_status parse_loaddedx(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
    char path[512];
    char extra;
    char *resolved_path;
    struct osh_material_loaddedx_table table;
    double *dedx_values;
    unsigned int z;
    size_t proj_idx;
    size_t e_idx;
    enum osh_status rc;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: LOADDEDX outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%511s %c", path, &extra) != 1) {
        osh_error("in %s line %i: LOADDEDX expects one path", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    resolved_path = NULL;
    if (osh_relative_path_to_file(&resolved_path, wm->wdir, path) != 0) {
        return OSH_ENOMEM;
    }

    memset(&table, 0, sizeof(table));
    rc = osh_material_loaddedx_table_load(resolved_path, &table);
    free(resolved_path);
    if (rc != OSH_OK) {
        return rc;
    }

    osh_material_dedx_clear(mat);

    proj_idx = 0u;
    while (proj_idx < table.nprojectiles) {
        rc = osh_material_loaddedx_projectile_z(&table, proj_idx, &z);
        if (rc != OSH_OK) {
            osh_material_loaddedx_table_free(&table);
            return rc;
        }
        dedx_values = (double *) malloc(table.nenergy * sizeof(*dedx_values));
        if (!dedx_values) {
            osh_material_loaddedx_table_free(&table);
            return OSH_ENOMEM;
        }
        e_idx = 0u;
        while (e_idx < table.nenergy) {
            dedx_values[e_idx] = (double) osh_material_loaddedx_mass_stopping_power(&table, proj_idx, e_idx);
            e_idx++;
        }
        rc = osh_material_dedx_set(mat, z, table.energy_grid, dedx_values, table.nenergy);
        free(dedx_values);
        if (rc != OSH_OK) {
            osh_material_loaddedx_table_free(&table);
            return rc;
        }
        proj_idx++;
    }

    osh_material_loaddedx_table_free(&table);

    return OSH_OK;
}

/**
 * @brief Handle the MATERIALI / MIVALUE / MIAV card.
 *
 * @details Sets the material-level mean excitation energy [eV]. This acts as a
 * fallback for transport projectiles not covered by an external dE/dx table.
 *
 * @param[in,out] wm    Material workspace.
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static enum osh_status
parse_material_mean_excitation_energy(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
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

/**
 * @brief Handle the ICRU card.
 *
 * @details Records an ICRU material id for lookup during post-parse completion.
 * Density, state, and elemental composition are filled from the embedded ICRU
 * database unless overridden by explicit cards.
 *
 * @param[in,out] wm    Material workspace.
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static enum osh_status parse_icru(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *mat;
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

/**
 * @brief Handle the MATERIAL card; push a new material slot onto the workspace.
 *
 * @param[in,out] wm    Material workspace to extend.
 * @param[in]     oshf  Open input file (for error location).
 * @param[in]     args  MATERIAL card argument: one name token.
 *
 * @returns OSH_OK on success, OSH_EPARSE on missing/duplicate name, OSH_ENOMEM on failure.
 */
static enum osh_status parse_material_start(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    return material_push(wm, oshf, args);
}

/**
 * @brief Handle the STATE card.
 *
 * @details Accepts an integer (1=condensed, 2=gas) or a keyword
 * ("condensed", "solid", "liquid", "gas"). Overrides the ICRU default state
 * when both are present.
 *
 * @param[in,out] wm    Material workspace.
 * @param[in]     oshf  Open input file.
 * @param[in]     args  Remainder of the input line after the key.
 *
 * @returns OSH_OK on success, OSH_EPARSE on unknown value.
 */
static enum osh_status parse_state(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    enum osh_status rc;
    struct osh_material *mat;
    char state_token[32];
    char extra;
    int state;
    int i;

    mat = material_current(wm);
    if (!mat) {
        osh_error("in %s line %i: STATE outside MATERIAL block", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (!args || sscanf(args, "%31s %c", state_token, &extra) != 1) {
        osh_error("in %s line %i: STATE expects an integer or keyword", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    i = 0;
    while (state_token[i] != '\0') {
        state_token[i] = (char) tolower((unsigned char) state_token[i]);
        i++;
    }

    rc = parse_state_value(&state, state_token);
    if (rc != OSH_OK) {
        osh_error("in %s line %i: unknown STATE '%s'", oshf->filename, oshf->lineno, state_token);
        return rc;
    }

    mat->state = state;
    return OSH_OK;
}

/**
 * @brief Return a pointer to the material currently being built (the last one).
 *
 * @param[in] wm  Material workspace.
 *
 * @returns Pointer to the last material in the array, or NULL if none exist.
 */
static struct osh_material *material_current(struct osh_material_workspace *wm) {
    if (!wm || wm->nmaterials == 0u) {
        return NULL;
    }

    return &wm->materials[wm->nmaterials - 1u];
}

/**
 * @brief Reset all fields of @p mat to their unset sentinel values.
 *
 * @details Negative doubles signal "unset" for rho and mean_excitation_energy
 * so that explicit user cards take precedence over ICRU defaults.
 *
 * @param[out] mat  Material to reset. Must point to allocated storage.
 */
static void material_defaults(struct osh_material *mat) {
    mat->elements = NULL;
    mat->dedx_overrides = NULL;
    mat->name = NULL;
    mat->rho = -1.0;
    mat->mean_excitation_energy = -1.0;
    mat->rgba[0] = 0.8f;
    mat->rgba[1] = 0.8f;
    mat->rgba[2] = 0.8f;
    mat->rgba[3] = 1.0f;
    mat->nelements = 0u;
    mat->ndedx_overrides = 0u;
    mat->lineno = 0u;
    mat->index = 0u;
    mat->icru_id = 0;
    mat->state = OSH_MATERIAL_STATE_UNSET;
}

/**
 * @brief Parse a single non-negative floating-point mean excitation energy [eV].
 *
 * @param[out] mean_excitation_energy_out  Receives the parsed value [eV].
 * @param[in]  oshf                        Open input file (for error location).
 * @param[in]  args                        Argument string to parse.
 * @param[in]  key_name                    Card name used in error messages.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static enum osh_status parse_mean_excitation_energy_value(double *mean_excitation_energy_out,
                                                          struct oshfile *oshf,
                                                          char const *args,
                                                          char const *key_name) {
    double mean_excitation_energy; /* [eV] */
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

/**
 * @brief Parse element composition arguments with optional isotope mass number.
 *
 * @details
 * `Z amount` means natural element and stores A=0. `Z A amount` means an
 * explicit isotope and requires A>0.
 */
static enum osh_status parse_element_card_args(unsigned int *z_out,
                                               unsigned int *a_out,
                                               double *amount_out,
                                               struct oshfile *oshf,
                                               char const *args,
                                               char const *key_name) {
    enum osh_status rc;
    char tok0[32];
    char tok1[32];
    char tok2[32];
    char extra[32];
    int ntok;
    unsigned int z;
    unsigned int a;
    double amount;

    if (!args) {
        osh_error(
            "in %s line %i: %s expects '<Z> <amount>' or '<Z> <A> <amount>'", oshf->filename, oshf->lineno, key_name);
        return OSH_EPARSE;
    }

    ntok = sscanf(args, "%31s %31s %31s %31s", tok0, tok1, tok2, extra);
    if (ntok != 2 && ntok != 3) {
        osh_error(
            "in %s line %i: %s expects '<Z> <amount>' or '<Z> <A> <amount>'", oshf->filename, oshf->lineno, key_name);
        return OSH_EPARSE;
    }

    rc = parse_uint_token(&z, tok0);
    if (rc != OSH_OK || z == 0u) {
        osh_error("in %s line %i: %s element Z must be a positive integer", oshf->filename, oshf->lineno, key_name);
        return OSH_EPARSE;
    }

    if (ntok == 2) {
        a = 0u;
        rc = parse_double_token(&amount, tok1);
    } else {
        rc = parse_uint_token(&a, tok1);
        if (rc != OSH_OK || a == 0u) {
            osh_error("in %s line %i: %s isotope A must be a positive integer when provided",
                      oshf->filename,
                      oshf->lineno,
                      key_name);
            return OSH_EPARSE;
        }
        rc = parse_double_token(&amount, tok2);
    }
    if (rc != OSH_OK || amount <= 0.0) {
        osh_error("in %s line %i: %s amount must be a positive number", oshf->filename, oshf->lineno, key_name);
        return OSH_EPARSE;
    }

    *z_out = z;
    *a_out = a;
    *amount_out = amount;
    return OSH_OK;
}

/**
 * @brief Parse a single double from a token string, rejecting trailing garbage.
 *
 * @param[out] value_out  Receives the parsed value.
 * @param[in]  token      NUL-terminated token to parse.
 *
 * @returns OSH_OK on success, OSH_EPARSE if the token is not a bare number.
 */
static enum osh_status parse_double_token(double *value_out, char const *token) {
    double value;
    char extra;

    if (sscanf(token, "%lf%c", &value, &extra) != 1) {
        return OSH_EPARSE;
    }
    *value_out = value;
    return OSH_OK;
}

/**
 * @brief Parse a single unsigned int from a token string, rejecting trailing garbage.
 *
 * @param[out] value_out  Receives the parsed value.
 * @param[in]  token      NUL-terminated token to parse.
 *
 * @returns OSH_OK on success, OSH_EPARSE if the token is not a bare integer.
 */
static enum osh_status parse_uint_token(unsigned int *value_out, char const *token) {
    unsigned int value;
    char extra;

    if (sscanf(token, "%u%c", &value, &extra) != 1) {
        return OSH_EPARSE;
    }
    *value_out = value;
    return OSH_OK;
}

/**
 * @brief Append a new material slot to the workspace and set its name.
 *
 * @details All fields are set to unset sentinels via material_defaults(). The
 * dense index is assigned as the current array length before growing it.
 *
 * @param[in,out] wm    Material workspace to grow.
 * @param[in]     oshf  Open input file (for error location).
 * @param[in]     args  MATERIAL card argument: one name token.
 *
 * @returns OSH_OK on success, OSH_EPARSE on missing/duplicate name, OSH_ENOMEM on failure.
 */
static enum osh_status material_push(struct osh_material_workspace *wm, struct oshfile *oshf, char const *args) {
    struct osh_material *materials;
    struct osh_material *mat;
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
    materials = (struct osh_material *) realloc(wm->materials, nmaterials * sizeof(*materials));
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

/**
 * @brief Append one element to a material's composition array.
 *
 * @details Exactly one of @p atom_count or @p mass_fraction must be positive;
 * the other must be negative (sentinel for "unset"). Mixed input modes within
 * a single material are rejected.
 *
 * @param[in,out] mat           Material to extend.
 * @param[in]     oshf          Open input file (for error location).
 * @param[in]     z             Atomic number (must be > 0).
 * @param[in]     a             Mass number; 0 = natural element.
 * @param[in]     atom_count    Relative atom count (dimensionless); < 0 if unset.
 * @param[in]     mass_fraction Mass fraction (dimensionless, 0–1); < 0 if unset.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input, OSH_ENOMEM on failure.
 */
static enum osh_status material_push_element(struct osh_material *mat,
                                             struct oshfile *oshf,
                                             unsigned int z,
                                             unsigned int a,
                                             double atom_count,
                                             double mass_fraction) {
    struct osh_material_element *elements;
    struct osh_material_element *elem;
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
    elements = (struct osh_material_element *) realloc(mat->elements, nelements * sizeof(*elements));
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

/**
 * @brief Duplicate @p src into a new heap allocation and store it at @p dst.
 *
 * @details Any previous allocation at @p *dst is freed before overwriting.
 *
 * @param[in,out] dst  Receives the newly allocated copy.
 * @param[in]     src  NUL-terminated source string.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
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

/**
 * @brief Interpret a lowercase state token as an osh_material_state value.
 *
 * @details Accepts numeric codes (1=condensed, 2=gas) and keywords
 * "condensed", "solid", "liquid", and "gas".
 *
 * @param[out] state_out  Receives the matched enum osh_material_state value.
 * @param[in]  token      Lowercase token to interpret.
 *
 * @returns OSH_OK on success, OSH_EPARSE if the token is not recognized.
 */
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

    if (strcmp(token, "condensed") == 0 || strcmp(token, "solid") == 0 || strcmp(token, "liquid") == 0) {
        *state_out = OSH_MATERIAL_STATE_CONDENSED;
        return OSH_OK;
    }
    if (strcmp(token, "gas") == 0) {
        *state_out = OSH_MATERIAL_STATE_GAS;
        return OSH_OK;
    }

    return OSH_EPARSE;
}
