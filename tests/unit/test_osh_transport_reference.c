#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_const.h"
#include "common/osh_rc.h"
#include "material/osh_material.h"
#include "material/osh_material_atomic_data.h"
#include "material/osh_material_icru.h"
#include "particle/osh_isotope_db.h"
#include "particle/osh_particle.h"
#include "physics/osh_physics_bethe.h"
#include "transport/prepare/osh_transport_material_prepare.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_transport_reference_test.tmp");
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static double test_material_element_mass_da(struct material_element const *el) {
    struct isotope iso;
    double mass_da;

    if (!el) {
        return 0.0;
    }

    if (el->a > 0u && osh_isotope_from_za(&iso, el->z, el->a)) {
        return iso.amass;
    }

    if (osh_material_natural_atomic_mass_da(el->z, &mass_da) == OSH_OK) {
        return mass_da;
    }

    return 2.0 * (double) el->z;
}

static void build_element_target(struct material_element const *el, struct osh_physics_bethe_target *tgt) {
    struct osh_material_icru_entry entry;
    double mass_da;

    mass_da = test_material_element_mass_da(el);

    tgt->z_mean = (double) el->z;
    tgt->a_mean = (mass_da > 0.0) ? mass_da : 2.0 * (double) el->z;
    tgt->i_value = el->mean_excitation_energy;

    if (osh_material_icru_lookup((int) el->z, &entry) == OSH_OK && entry.rho > 0.0) {
        tgt->rho = entry.rho;
    } else {
        tgt->rho = 1.0;
    }
}

static double
eval_element_sum_bethe(struct material const *mat, unsigned int proj_z, unsigned int proj_a, double energy_mev_per_u) {
    struct osh_physics_bethe_projectile proj;
    size_t i;
    double sp_sum;

    proj.z = (double) proj_z;
    proj.a = (double) proj_a;
    if (!osh_particle_nuclear_mass_mev_from_za(proj_z, proj_a, &proj.mass_mev)) {
        proj.mass_mev = (double) proj_a * OSH_AMU; /* fallback */
    }
    sp_sum = 0.0;

    for (i = 0; i < mat->nelements; ++i) {
        struct osh_physics_bethe_target tgt;
        struct osh_physics_bethe_sewn sewn;
        double sp_i;

        build_element_target(&mat->elements[i], &tgt);
        osh_physics_bethe_sewn_compute(&proj, &tgt, &sewn);
        sp_i = osh_physics_bethe_eval(energy_mev_per_u, &proj, &tgt, &sewn);
        if (sp_i > 0.0) {
            sp_sum += mat->elements[i].mass_fraction * sp_i;
        }
    }

    return sp_sum;
}

static void print_projectile_slice(struct osh_transport_material_runtime const *tables,
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
        sp = osh_transport_material_sp_lookup(tables, mat_idx, proj_idx, e);
        range = osh_transport_material_range_lookup(tables, mat_idx, proj_idx, e);
        printf("  E=%9.4f MeV/u  dEdx=%12.6f MeV cm^2/g  range=%12.6f g/cm^2\n", e, sp, range);
    }
}

static void print_projectile_slice_with_diff(struct osh_transport_material_runtime const *tables,
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
        double d_sp_pct;
        double d_range_pct;

        e = energies[i];
        sp = osh_transport_material_sp_lookup(tables, mat_idx, proj_idx, e);
        range = osh_transport_material_range_lookup(tables, mat_idx, proj_idx, e);
        ref_sp = osh_transport_material_sp_lookup(tables, ref_mat_idx, proj_idx, e);
        ref_range = osh_transport_material_range_lookup(tables, ref_mat_idx, proj_idx, e);
        d_sp = sp - ref_sp;
        d_sp_pct = (ref_sp != 0.0) ? 100.0 * d_sp / ref_sp : 0.0;
        d_range_pct = (ref_range != 0.0) ? 100.0 * (range - ref_range) / ref_range : 0.0;

        printf("  E=%9.4f MeV/u  dEdx=%12.6f MeV cm^2/g  range=%12.6f g/cm^2  dSP=%+10.6f  dSP%%=%+8.3f %%  dR=%+8.3f "
               "%%\n",
               e,
               sp,
               range,
               d_sp,
               d_sp_pct,
               d_range_pct);
    }
}

static void print_projectile_direct_bethe_with_diff(struct osh_transport_material_runtime const *tables,
                                                    struct material const *mat,
                                                    size_t ref_mat_idx,
                                                    size_t proj_idx,
                                                    char const *label,
                                                    double const *energies,
                                                    size_t nenergies) {
    size_t i;
    unsigned int proj_z;
    unsigned int proj_a;

    proj_z = tables->projectile_z[proj_idx];
    proj_a = tables->projectile_a[proj_idx];

    printf("%s\n", label);
    for (i = 0; i < nenergies; ++i) {
        double e;
        double sp_direct;
        double ref_sp;
        double d_sp;
        double d_sp_pct;

        e = energies[i];
        sp_direct = eval_element_sum_bethe(mat, proj_z, proj_a, e);
        ref_sp = osh_transport_material_sp_lookup(tables, ref_mat_idx, proj_idx, e);
        d_sp = sp_direct - ref_sp;
        d_sp_pct = (ref_sp != 0.0) ? 100.0 * d_sp / ref_sp : 0.0;

        printf("  E=%9.4f MeV/u  dEdx=%12.6f MeV cm^2/g  dSP=%+10.6f  dSP%%=%+8.3f %%\n", e, sp_direct, d_sp, d_sp_pct);
    }
}

int main(void) {
    char mat_path[512];
    char water_path[512];
    char mat_text[4096];
    struct material_workspace *wm;
    struct osh_transport_material_runtime tables;
    struct material const *water_loaded;
    struct material const *water_bethe;
    enum osh_status rc;
    double const energies[] = {0.03, 0.10, 1.0, 10.0, 100.0, 400.0};

    wm = NULL;
    memset(&tables, 0, sizeof(tables));

    snprintf(water_path, sizeof(water_path), "%s/tests/cases/04_simple_loaddedx/Water.txt", OSH_PROJECT_SOURCE_DIR);
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
    write_temp_file(mat_path, sizeof(mat_path), mat_text);

    rc = osh_material_setup_from_path(mat_path, NULL, &wm);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    rc = osh_transport_material_prepare(wm, 6u, &tables);
    ASSERT_TRUE(rc == OSH_OK);

    water_loaded = osh_material_by_name(wm, "WaterLoaded");
    water_bethe = osh_material_by_name(wm, "WaterBethe");
    ASSERT_TRUE(water_loaded != NULL);
    ASSERT_TRUE(water_bethe != NULL);

    printf("Water reference values from runtime transport tables\n");
    printf("Table grid: %zu points from %.4f to %.1f MeV/u\n", tables.nenergy, tables.emin, tables.emax);
    printf("Source 1: LOADDEDX %s\n", water_path);
    printf("Source 2: Bethe fallback from assembled material data\n");
    printf("WaterBethe inferred element MEE values:\n");
    printf("  element[0]: Z=%u A=%u I=%.6f eV\n",
           water_bethe->elements[0].z,
           water_bethe->elements[0].a,
           water_bethe->elements[0].mean_excitation_energy);
    printf("  element[1]: Z=%u A=%u I=%.6f eV\n",
           water_bethe->elements[1].z,
           water_bethe->elements[1].a,
           water_bethe->elements[1].mean_excitation_energy);
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

    printf("\n");
    printf("Material index: %zu (%s direct)\n", water_bethe->index, water_bethe->name);
    print_projectile_direct_bethe_with_diff(&tables,
                                            water_bethe,
                                            water_bethe->index,
                                            0u,
                                            "Proton on WaterBethe (H-1, direct Bethe vs runtime Bethe)",
                                            energies,
                                            sizeof(energies) / sizeof(energies[0]));
    print_projectile_direct_bethe_with_diff(&tables,
                                            water_bethe,
                                            water_bethe->index,
                                            1u,
                                            "Alpha on WaterBethe (He-4, direct Bethe vs runtime Bethe)",
                                            energies,
                                            sizeof(energies) / sizeof(energies[0]));
    print_projectile_direct_bethe_with_diff(&tables,
                                            water_bethe,
                                            water_bethe->index,
                                            5u,
                                            "Carbon on WaterBethe (C-12, direct Bethe vs runtime Bethe)",
                                            energies,
                                            sizeof(energies) / sizeof(energies[0]));

    osh_transport_material_runtime_free(&tables);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(mat_path) == 0);
    return 0;
}
