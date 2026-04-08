#include "material/osh_material_loaddedx.h"

#include <errno.h>
#include <stdlib.h>

#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "particle/osh_particle.h"

enum { OSH_MATERIAL_LOADDEDX_MINPROJECTILES = 18 };

static void table_reset(struct osh_material_loaddedx_table *table) {
    if (!table) {
        return;
    }
    table->energy_grid = NULL;
    table->mass_stopping_power = NULL;
    table->projectile_z = NULL;
    table->projectile_a = NULL;
    table->nprojectiles = 0u;
    table->nenergy = 0u;
}

static enum osh_status
parse_numeric_row(char const *line, double *energy_out, double **values_out, size_t *nvalues_out) {
    char *cursor;
    char *endptr;
    double *values;
    double value;
    size_t nvalues;

    if (!line || !energy_out || !values_out || !nvalues_out) {
        return OSH_EINVAL;
    }

    errno = 0;
    cursor = (char *) line;
    value = strtod(cursor, &endptr);
    if (endptr == cursor || errno != 0) {
        return OSH_EPARSE;
    }
    *energy_out = value;
    cursor = endptr;

    values = NULL;
    nvalues = 0u;

    while (1) {
        errno = 0;
        value = strtod(cursor, &endptr);
        if (errno != 0 || endptr == cursor) {
            break;
        }
        {
            double *new_values;

            new_values = realloc(values, (nvalues + 1u) * sizeof(double));
            if (!new_values) {
                free(values);
                return OSH_ENOMEM;
            }
            values = new_values;
        }
        values[nvalues] = value;
        nvalues += 1u;
        cursor = endptr;
    }

    while (*cursor == ' ' || *cursor == '\t') {
        cursor++;
    }
    if (*cursor != '\0') {
        free(values);
        return OSH_EPARSE;
    }
    if (nvalues < OSH_MATERIAL_LOADDEDX_MINPROJECTILES) {
        free(values);
        return OSH_EPARSE;
    }

    *values_out = values;
    *nvalues_out = nvalues;
    return OSH_OK;
}

static enum osh_status append_row(
    double **energy_rows, float **value_rows, size_t nprojectiles, size_t *nrows, double energy, double const *values) {
    double *new_energy_rows;
    float *new_value_rows;
    size_t i;

    new_energy_rows = realloc(*energy_rows, (*nrows + 1u) * sizeof(double));
    if (!new_energy_rows) {
        return OSH_ENOMEM;
    }
    *energy_rows = new_energy_rows;

    new_value_rows = realloc(*value_rows, (*nrows + 1u) * nprojectiles * sizeof(float));
    if (!new_value_rows) {
        return OSH_ENOMEM;
    }
    *value_rows = new_value_rows;

    (*energy_rows)[*nrows] = energy;
    for (i = 0; i < nprojectiles; i++) {
        (*value_rows)[(*nrows) * nprojectiles + i] = (float) values[i];
    }
    *nrows += 1u;

    return OSH_OK;
}

static enum osh_status init_projectile_map(struct osh_material_loaddedx_table *table, size_t nprojectiles) {
    size_t i;

    table->projectile_z = malloc(nprojectiles * sizeof(unsigned int));
    table->projectile_a = malloc(nprojectiles * sizeof(unsigned int));
    if (!table->projectile_z || !table->projectile_a) {
        return OSH_ENOMEM;
    }

    for (i = 0; i < nprojectiles; i++) {
        unsigned int z;
        unsigned int a;

        z = (unsigned int) (i + 1u);
        if (!osh_particle_default_isotope_a(z, &a)) {
            return OSH_EPARSE;
        }
        table->projectile_z[i] = z;
        table->projectile_a[i] = a;
    }

    return OSH_OK;
}

enum osh_status osh_material_loaddedx_table_load(char const *path, struct osh_material_loaddedx_table *table) {
    struct oshfile *oshf;
    char *line;
    double *energy_rows;
    float *value_rows;
    double *row_values;
    double row_energy;
    float *projectile_major_values;
    size_t nprojectiles;
    size_t i;
    size_t j;
    size_t nrows;
    int lineno;
    int rc_read;
    enum osh_status rc;

    if (!path || !table) {
        return OSH_EINVAL;
    }

    osh_material_loaddedx_table_free(table);
    table_reset(table);

    oshf = osh_fopen(path);
    if (!oshf) {
        return OSH_EIO;
    }

    line = NULL;
    energy_rows = NULL;
    value_rows = NULL;
    row_values = NULL;
    nprojectiles = 0u;
    nrows = 0u;
    rc = OSH_OK;

    while ((rc_read = osh_readline(oshf, &line, &lineno)) >= 0) {
        if (rc_read == 0) {
            free(line);
            line = NULL;
            continue;
        }

        rc = parse_numeric_row(line, &row_energy, &row_values, &nprojectiles);
        if (rc == OSH_OK && table->nprojectiles == 0u) {
            table->nprojectiles = nprojectiles;
            rc = init_projectile_map(table, nprojectiles);
        } else if (rc == OSH_OK && nprojectiles != table->nprojectiles) {
            rc = OSH_EPARSE;
        }
        free(line);
        line = NULL;
        if (rc != OSH_OK) {
            osh_error("in %s line %d: expected 1 energy column and at least %d contiguous stopping-power columns",
                      path,
                      lineno,
                      OSH_MATERIAL_LOADDEDX_MINPROJECTILES);
            free(row_values);
            row_values = NULL;
            break;
        }
        if (nrows > 0u && row_energy <= energy_rows[nrows - 1u]) {
            osh_error("in %s line %d: energy grid must be strictly increasing", path, lineno);
            rc = OSH_EPARSE;
            free(row_values);
            row_values = NULL;
            break;
        }

        rc = append_row(&energy_rows, &value_rows, table->nprojectiles, &nrows, row_energy, row_values);
        free(row_values);
        row_values = NULL;
        if (rc != OSH_OK) {
            break;
        }
    }

    if (rc == OSH_OK && nrows == 0u) {
        osh_error("in %s: no numeric stopping-power rows found", path);
        rc = OSH_EPARSE;
    }

    if (rc == OSH_OK) {
        projectile_major_values = malloc(table->nprojectiles * nrows * sizeof(float));
        if (!projectile_major_values) {
            rc = OSH_ENOMEM;
        } else {
            for (j = 0; j < table->nprojectiles; j++) {
                for (i = 0; i < nrows; i++) {
                    projectile_major_values[j * nrows + i] = value_rows[i * table->nprojectiles + j];
                }
            }
            table->energy_grid = energy_rows;
            table->mass_stopping_power = projectile_major_values;
            table->nenergy = nrows;
            energy_rows = NULL;
        }
    }

    free(line);
    free(energy_rows);
    free(value_rows);
    free(row_values);
    if (rc != OSH_OK) {
        osh_material_loaddedx_table_free(table);
    }
    osh_fclose(oshf);

    return rc;
}

void osh_material_loaddedx_table_free(struct osh_material_loaddedx_table *table) {
    if (!table) {
        return;
    }
    free(table->energy_grid);
    free(table->mass_stopping_power);
    free(table->projectile_z);
    free(table->projectile_a);
    table_reset(table);
}

enum osh_status osh_material_loaddedx_projectile_za(struct osh_material_loaddedx_table const *table,
                                                    size_t projectile_idx,
                                                    unsigned int *z_out,
                                                    unsigned int *a_out) {
    if (!table || !z_out || !a_out) {
        return OSH_EINVAL;
    }
    if (projectile_idx >= table->nprojectiles) {
        return OSH_EINVAL;
    }
    *z_out = table->projectile_z[projectile_idx];
    *a_out = table->projectile_a[projectile_idx];
    return OSH_OK;
}
