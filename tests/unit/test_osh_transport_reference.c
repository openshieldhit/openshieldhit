#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_rc.h"
#include "material/osh_material.h"
#include "transport/prepare/osh_transport_prepare.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void print_projectile_slice(struct osh_transport_runtime_tables const *tables,
                                   size_t mat_idx,
                                   size_t proj_idx,
                                   char const *label,
                                   double const *energies,
                                   size_t nenergies) {
    size_t i;

    printf("%s\n", label);
    for (i = 0; i < nenergies; ++i) {
        double e;
        double sp;
        double range;

        e = energies[i];
        sp = osh_transport_sp_lookup(tables, mat_idx, proj_idx, e);
        range = osh_transport_range_lookup(tables, mat_idx, proj_idx, e);
        printf("  E=%9.4f MeV/u  dEdx=%12.6f MeV cm^2/g  range=%12.6f g/cm^2\n", e, sp, range);
    }
}

static void print_projectile_slice_with_diff(struct osh_transport_runtime_tables const *tables,
                                             size_t mat_idx,
                                             size_t ref_mat_idx,
                                             size_t proj_idx,
                                             char const *label,
                                             double const *energies,
                                             size_t nenergies) {
    size_t i;

    printf("%s\n", label);
    for (i = 0; i < nenergies; ++i) {
        double e;
        double sp;
        double range;
        double ref_sp;
        double ref_range;
        double d_sp;
        double d_range_pct;

        e = energies[i];
        sp = osh_transport_sp_lookup(tables, mat_idx, proj_idx, e);
        range = osh_transport_range_lookup(tables, mat_idx, proj_idx, e);
        ref_sp = osh_transport_sp_lookup(tables, ref_mat_idx, proj_idx, e);
        ref_range = osh_transport_range_lookup(tables, ref_mat_idx, proj_idx, e);
        d_sp = sp - ref_sp;
        d_range_pct = (ref_range != 0.0) ? 100.0 * (range - ref_range) / ref_range : 0.0;

        printf("  E=%9.4f MeV/u  dEdx=%12.6f MeV cm^2/g  range=%12.6f g/cm^2  dSP=%+10.6f  dR=%+8.3f %%\n",
               e,
               sp,
               range,
               d_sp,
               d_range_pct);
    }
}

int main(void) {
    char mat_path[512];
    char water_path[512];
    char mat_text[4096];
    struct material_workspace *wm;
    struct osh_transport_runtime_tables tables;
    struct material const *water_loaded;
    struct material const *water_bethe;
    enum osh_status rc;
    double const energies[] = {0.03, 0.10, 1.0, 10.0, 100.0, 400.0};

    wm = NULL;
    memset(&tables, 0, sizeof(tables));

    snprintf(water_path, sizeof(water_path), "%s/tests/cases/04_simple_loaddedx/Water.txt", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/fixtures/test_transport_reference_water.tmp", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_text,
             sizeof(mat_text),
             "MATERIAL WaterLoaded\n"
             "ICRU 276\n"
             "LOADDEDX %s\n",
             water_path);
    snprintf(mat_text + strlen(mat_text),
             sizeof(mat_text) - strlen(mat_text),
             "MATERIAL WaterBethe\n"
             "ICRU 276\n");

    {
        FILE *fp;
        fp = fopen(mat_path, "w");
        ASSERT_TRUE(fp != NULL);
        ASSERT_TRUE(fputs(mat_text, fp) >= 0);
        ASSERT_TRUE(fclose(fp) == 0);
    }

    rc = osh_material_setup_from_path(mat_path, NULL, &wm);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    rc = osh_transport_prepare(wm, 6u, &tables);
    ASSERT_TRUE(rc == OSH_OK);

    water_loaded = osh_material_by_name(wm, "WaterLoaded");
    water_bethe = osh_material_by_name(wm, "WaterBethe");
    ASSERT_TRUE(water_loaded != NULL);
    ASSERT_TRUE(water_bethe != NULL);

    printf("Water reference values from runtime transport tables\n");
    printf("Table grid: %zu points from %.4f to %.1f MeV/u\n", tables.nenergy, tables.emin, tables.emax);
    printf("Source 1: LOADDEDX %s\n", water_path);
    printf("Source 2: Bethe fallback from assembled material data\n");
    printf("\n");

    printf("Material index: %zu (%s)\n", water_loaded->index, water_loaded->name);
    print_projectile_slice(&tables,
                           water_loaded->index,
                           0u,
                           "Proton on WaterLoaded (H-1, LOADDEDX)",
                           energies,
                           sizeof(energies) / sizeof(energies[0]));
    print_projectile_slice(&tables,
                           water_loaded->index,
                           1u,
                           "Alpha on WaterLoaded (He-4, LOADDEDX)",
                           energies,
                           sizeof(energies) / sizeof(energies[0]));
    print_projectile_slice(&tables,
                           water_loaded->index,
                           5u,
                           "Carbon on WaterLoaded (C-12, LOADDEDX)",
                           energies,
                           sizeof(energies) / sizeof(energies[0]));

    printf("\n");
    printf("Material index: %zu (%s)\n", water_bethe->index, water_bethe->name);
    print_projectile_slice_with_diff(&tables,
                                     water_bethe->index,
                                     water_loaded->index,
                                     0u,
                                     "Proton on WaterBethe (H-1, Bethe vs LOADDEDX)",
                                     energies,
                                     sizeof(energies) / sizeof(energies[0]));
    print_projectile_slice_with_diff(&tables,
                                     water_bethe->index,
                                     water_loaded->index,
                                     1u,
                                     "Alpha on WaterBethe (He-4, Bethe vs LOADDEDX)",
                                     energies,
                                     sizeof(energies) / sizeof(energies[0]));
    print_projectile_slice_with_diff(&tables,
                                     water_bethe->index,
                                     water_loaded->index,
                                     5u,
                                     "Carbon on WaterBethe (C-12, Bethe vs LOADDEDX)",
                                     energies,
                                     sizeof(energies) / sizeof(energies[0]));

    osh_transport_runtime_tables_free(&tables);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(mat_path) == 0);
    return 0;
}
