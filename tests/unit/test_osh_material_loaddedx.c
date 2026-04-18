/*
 * test_osh_material_loaddedx.c
 *
 * Unit tests for legacy/external LOADDEDX table loading.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_material_loaddedx.h"
#include "openshieldhit/status.h"

#define FTOL 1e-5

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((double) (a) - (double) (b)) < (tol))

static void fixture_path(char *buf, size_t buflen, char const *fname) {
    snprintf(buf, buflen, "%s/tests/cases/04_simple_loaddedx/%s", OSH_PROJECT_SOURCE_DIR, fname);
}

static void project_fixture_path(char *buf, size_t buflen, char const *fname) {
    snprintf(buf, buflen, "%s/tests/fixtures/%s", OSH_PROJECT_SOURCE_DIR, fname);
}

static void test_load_lucite_numeric_only(void) {
    struct osh_material_loaddedx_table table;
    char path[512];
    unsigned int z;

    memset(&table, 0, sizeof(table));
    fixture_path(path, sizeof(path), "Lucite.txt");

    ASSERT_TRUE(osh_material_loaddedx_table_load(path, &table) == OSH_OK);
    ASSERT_TRUE(table.nprojectiles == 18u);
    ASSERT_TRUE(table.nenergy == 53u);
    ASSERT_NEAR(table.energy_grid[0], 0.025, 1e-12);
    ASSERT_NEAR(table.energy_grid[table.nenergy - 1u], 1000.0, 1e-9);
    ASSERT_NEAR(osh_material_loaddedx_mass_stopping_power(&table, 0u, 0u), 743.21, 1e-3);
    ASSERT_NEAR(osh_material_loaddedx_mass_stopping_power(&table, 5u, 0u), 3958.0, 1e-3);
    ASSERT_NEAR(osh_material_loaddedx_mass_stopping_power(&table, 17u, table.nenergy - 1u), 714.4, 1e-3);
    ASSERT_TRUE(osh_material_loaddedx_projectile_z(&table, 5u, &z) == OSH_OK);
    ASSERT_TRUE(z == 6u);

    osh_material_loaddedx_table_free(&table);
}

static void test_load_water_with_comment_header(void) {
    struct osh_material_loaddedx_table table;
    char path[512];

    memset(&table, 0, sizeof(table));
    fixture_path(path, sizeof(path), "Water.txt");

    ASSERT_TRUE(osh_material_loaddedx_table_load(path, &table) == OSH_OK);
    ASSERT_TRUE(table.nprojectiles == 18u);
    ASSERT_TRUE(table.nenergy == 53u);
    ASSERT_NEAR(osh_material_loaddedx_mass_stopping_power(&table, 0u, 0u), 619.82, 1e-3);
    ASSERT_NEAR(osh_material_loaddedx_mass_stopping_power(&table, 1u, 0u), 1131.0, 1e-3);
    ASSERT_NEAR(osh_material_loaddedx_mass_stopping_power(&table, 17u, table.nenergy - 1u), 725.8, 1e-3);

    osh_material_loaddedx_table_free(&table);
}

static void test_all_fixture_grids_match(void) {
    struct osh_material_loaddedx_table water;
    struct osh_material_loaddedx_table lucite;
    struct osh_material_loaddedx_table air;
    struct osh_material_loaddedx_table alanine;
    char path[512];
    size_t i;

    memset(&water, 0, sizeof(water));
    memset(&lucite, 0, sizeof(lucite));
    memset(&air, 0, sizeof(air));
    memset(&alanine, 0, sizeof(alanine));

    fixture_path(path, sizeof(path), "Water.txt");
    ASSERT_TRUE(osh_material_loaddedx_table_load(path, &water) == OSH_OK);
    fixture_path(path, sizeof(path), "Lucite.txt");
    ASSERT_TRUE(osh_material_loaddedx_table_load(path, &lucite) == OSH_OK);
    fixture_path(path, sizeof(path), "Air.txt");
    ASSERT_TRUE(osh_material_loaddedx_table_load(path, &air) == OSH_OK);
    fixture_path(path, sizeof(path), "Alanine.txt");
    ASSERT_TRUE(osh_material_loaddedx_table_load(path, &alanine) == OSH_OK);

    ASSERT_TRUE(water.nenergy == lucite.nenergy);
    ASSERT_TRUE(water.nenergy == air.nenergy);
    ASSERT_TRUE(water.nenergy == alanine.nenergy);

    for (i = 0; i < water.nenergy; i++) {
        ASSERT_NEAR(water.energy_grid[i], lucite.energy_grid[i], 1e-12);
        ASSERT_NEAR(water.energy_grid[i], air.energy_grid[i], 1e-12);
        ASSERT_NEAR(water.energy_grid[i], alanine.energy_grid[i], 1e-12);
    }

    osh_material_loaddedx_table_free(&water);
    osh_material_loaddedx_table_free(&lucite);
    osh_material_loaddedx_table_free(&air);
    osh_material_loaddedx_table_free(&alanine);
}

static void test_load_extended_contiguous_columns(void) {
    struct osh_material_loaddedx_table table;
    char path[512];
    unsigned int z;

    memset(&table, 0, sizeof(table));
    project_fixture_path(path, sizeof(path), "loaddedx_extended.txt");

    ASSERT_TRUE(osh_material_loaddedx_table_load(path, &table) == OSH_OK);
    ASSERT_TRUE(table.nprojectiles == 20u);
    ASSERT_TRUE(table.nenergy == 2u);
    ASSERT_TRUE(osh_material_loaddedx_projectile_z(&table, 18u, &z) == OSH_OK);
    ASSERT_TRUE(z == 19u);
    ASSERT_TRUE(osh_material_loaddedx_projectile_z(&table, 19u, &z) == OSH_OK);
    ASSERT_TRUE(z == 20u);
    ASSERT_NEAR(osh_material_loaddedx_mass_stopping_power(&table, 19u, 1u), 2.0, 1e-6);

    osh_material_loaddedx_table_free(&table);
}

int main(void) {
    test_load_lucite_numeric_only();
    test_load_water_with_comment_header();
    test_all_fixture_grids_match();
    test_load_extended_contiguous_columns();
    return 0;
}
