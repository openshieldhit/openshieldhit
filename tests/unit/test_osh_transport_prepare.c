#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((double) (a) - (double) (b)) < (tol))

static int tmp_counter = 0;

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_transport_prepare_test_%d.tmp", tmp_counter++);
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

static void build_effective_target(struct material const *mat, struct osh_physics_bethe_target *tgt) {
    size_t i;
    double sum_wz_over_a;
    double sum_wz;
    double mee_ln_sum;

    sum_wz_over_a = 0.0;
    sum_wz = 0.0;
    mee_ln_sum = 0.0;

    for (i = 0; i < mat->nelements; ++i) {
        struct material_element const *el;
        double a_i;
        double z_i;
        double w_i;

        el = &mat->elements[i];
        a_i = test_material_element_mass_da(el);
        z_i = (double) el->z;
        w_i = el->mass_fraction;

        if (a_i <= 0.0) {
            continue;
        }

        sum_wz_over_a += w_i * z_i / a_i;
        sum_wz += w_i * z_i;

        if (el->mean_excitation_energy > 0.0) {
            mee_ln_sum += w_i * (z_i / a_i) * log(el->mean_excitation_energy);
        }
    }

    if (sum_wz_over_a <= 0.0) {
        sum_wz_over_a = 1.0;
    }
    if (sum_wz <= 0.0) {
        sum_wz = 1.0;
    }

    tgt->z_mean = sum_wz;
    tgt->a_mean = sum_wz / sum_wz_over_a;
    tgt->rho = (mat->rho > 0.0) ? mat->rho : 1.0;
    if (mat->mean_excitation_energy > 0.0) {
        tgt->i_value = mat->mean_excitation_energy;
    } else if (mee_ln_sum != 0.0) {
        tgt->i_value = exp(mee_ln_sum / sum_wz_over_a);
    } else {
        tgt->i_value = 78.0;
    }
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
eval_effective_bethe(struct material const *mat, unsigned int proj_z, unsigned int proj_a, double energy_mev_per_u) {
    struct osh_physics_bethe_target tgt;
    struct osh_physics_bethe_projectile proj;
    struct osh_physics_bethe_sewn sewn;

    build_effective_target(mat, &tgt);
    proj.z = (double) proj_z;
    proj.a = (double) proj_a;
    if (!osh_particle_nuclear_mass_mev_from_za(proj_z, proj_a, &proj.mass_mev)) {
        proj.mass_mev = 0.0;
    }
    ASSERT_TRUE(proj.mass_mev > 0.0);
    osh_physics_bethe_sewn_compute(&proj, &tgt, &sewn);
    return osh_physics_bethe_eval(energy_mev_per_u, &proj, &tgt, &sewn);
}

static double
eval_element_sum_bethe(struct material const *mat, unsigned int proj_z, unsigned int proj_a, double energy_mev_per_u) {
    struct osh_physics_bethe_projectile proj;
    size_t i;
    double sp_sum;

    proj.z = (double) proj_z;
    proj.a = (double) proj_a;
    if (!osh_particle_nuclear_mass_mev_from_za(proj_z, proj_a, &proj.mass_mev)) {
        proj.mass_mev = 0.0;
    }
    ASSERT_TRUE(proj.mass_mev > 0.0);
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

static void test_prepare_mixed_loaddedx_and_bethe_paths(void) {
    char mat_path[512];
    char water_path[512];
    char ext20_path[512];
    char mat_text[4096];
    struct material_workspace *wm;
    struct osh_transport_material_runtime tables;
    size_t idx_table18;
    size_t idx_table20;
    size_t idx_bethe;
    enum osh_status rc;

    wm = NULL;
    memset(&tables, 0, sizeof(tables));

    snprintf(water_path, sizeof(water_path), "%s/tests/cases/04_simple_loaddedx/Water.txt", OSH_PROJECT_SOURCE_DIR);
    snprintf(ext20_path, sizeof(ext20_path), "%s/tests/fixtures/loaddedx_extended.txt", OSH_PROJECT_SOURCE_DIR);

    snprintf(mat_text,
             sizeof(mat_text),
             "MATERIAL Table18\n"
             "RHO 1.0\n"
             "LOADDEDX %s\n"
             "ELEMENT 1 1\n"
             "MATERIAL Table20\n"
             "RHO 1.0\n"
             "LOADDEDX %s\n"
             "ELEMENT 1 1\n"
             "MATERIAL BetheWater\n"
             "ICRU 276\n",
             water_path,
             ext20_path);
    write_temp_file(mat_path, sizeof(mat_path), mat_text);

    rc = osh_material_setup_from_path(mat_path, NULL, &wm);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    rc = osh_transport_material_prepare(wm, 25u, &tables);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(tables.nmaterials == 5u);
    ASSERT_TRUE(tables.nprojectiles == 25u);
    ASSERT_TRUE(tables.nenergy == (size_t) OSH_TRANSPORT_MATERIAL_NENERGY);
    ASSERT_NEAR(tables.emin, OSH_TRANSPORT_MATERIAL_EMIN, 1e-12);
    ASSERT_NEAR(tables.emax, OSH_TRANSPORT_MATERIAL_EMAX, 1e-12);

    idx_table18 = osh_material_by_name(wm, "Table18")->index;
    idx_table20 = osh_material_by_name(wm, "Table20")->index;
    idx_bethe = osh_material_by_name(wm, "BetheWater")->index;

    ASSERT_NEAR(
        osh_transport_material_sp_lookup(&tables, OSH_MATERIAL_INDEX_BLACKHOLE, 0u, OSH_TRANSPORT_MATERIAL_EMIN),
        0.0,
        1e-12);
    ASSERT_NEAR(osh_transport_material_sp_lookup(&tables, OSH_MATERIAL_INDEX_VACUUM, 0u, OSH_TRANSPORT_MATERIAL_EMIN),
                0.0,
                1e-12);

    ASSERT_NEAR(osh_transport_material_sp_lookup(&tables, idx_table18, 0u, 0.025), 619.82, 1e-2);
    ASSERT_NEAR(osh_transport_material_sp_lookup(&tables, idx_table18, 0u, 0.0),
                osh_transport_material_sp_lookup(&tables, idx_table18, 0u, tables.emin),
                1e-12);
    ASSERT_NEAR(osh_transport_material_sp_lookup(&tables, idx_table18, 0u, -1.0),
                osh_transport_material_sp_lookup(&tables, idx_table18, 0u, tables.emin),
                1e-12);
    ASSERT_NEAR(osh_transport_material_sp_lookup(&tables, idx_table18, 0u, 2000.0),
                osh_transport_material_sp_lookup(&tables, idx_table18, 0u, tables.emax),
                1e-12);
    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_table18, 19u, 0.025) > 0.0);
    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_table18, 19u, 10.0) > 0.0);

    ASSERT_NEAR(osh_transport_material_sp_lookup(&tables, idx_table20, 19u, 0.025), 1.0, 1e-6);
    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_table20, 19u, 0.03) > 1.0);
    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_table20, 19u, 0.03) < 2.0);
    ASSERT_NEAR(osh_transport_material_sp_lookup(&tables, idx_table20, 19u, 1.0), 2.0, 1e-6);
    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_table18, 24u, 10.0) > 0.0);
    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_table20, 24u, 10.0) > 0.0);
    ASSERT_NEAR(osh_transport_material_range_lookup(&tables, idx_table20, 19u, 0.05)
                    / osh_transport_material_range_lookup(&tables, idx_table20, 0u, 0.05),
                (double) tables.projectile_a[19u] / (double) tables.projectile_a[0u],
                1e-4);

    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_bethe, 0u, 1.0) > 0.0);
    ASSERT_TRUE(osh_transport_material_sp_lookup(&tables, idx_bethe, 5u, 10.0) > 0.0);

    ASSERT_TRUE(tables.range_csda != NULL);
    ASSERT_TRUE(tables.range_csda[(idx_bethe * tables.nprojectiles + 0u) * tables.nenergy + 0u] == 0.0f);
    ASSERT_TRUE(tables.range_csda[(idx_bethe * tables.nprojectiles + 0u) * tables.nenergy + 1u] >= 0.0f);
    ASSERT_NEAR(osh_transport_material_range_lookup(&tables, idx_table18, 0u, 0.0),
                osh_transport_material_range_lookup(&tables, idx_table18, 0u, tables.emin),
                1e-12);
    ASSERT_NEAR(osh_transport_material_range_lookup(&tables, idx_table18, 0u, 2000.0),
                osh_transport_material_range_lookup(&tables, idx_table18, 0u, tables.emax),
                1e-12);

    osh_transport_material_runtime_free(&tables);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(mat_path) == 0);
}

static void test_prepare_bethe_infers_element_mee_from_material_mean_excitation(void) {
    char mat_path[512];
    char mat_text[4096];
    struct material_workspace *wm;
    struct osh_transport_material_runtime tables;
    struct material const *mat_elements;
    struct material const *mat_effective;
    double expected_elements;
    double expected_effective;
    double got_elements;
    double got_effective;
    enum osh_status rc;

    wm = NULL;
    memset(&tables, 0, sizeof(tables));

    snprintf(mat_text,
             sizeof(mat_text),
             "MATERIAL WaterElemI\n"
             "RHO 1.0\n"
             "ELEMENT 1 2\n"
             "IVALUE 22.9\n"
             "ELEMENT 8 1\n"
             "\n"
             "MATERIAL WaterMatI\n"
             "ICRU 276\n");
    write_temp_file(mat_path, sizeof(mat_path), mat_text);

    rc = osh_material_setup_from_path(mat_path, NULL, &wm);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    mat_elements = osh_material_by_name(wm, "WaterElemI");
    mat_effective = osh_material_by_name(wm, "WaterMatI");
    ASSERT_TRUE(mat_elements != NULL);
    ASSERT_TRUE(mat_effective != NULL);
    ASSERT_TRUE(mat_elements->nelements == 2u);
    ASSERT_TRUE(mat_effective->nelements == 2u);
    ASSERT_TRUE(mat_elements->elements[0].mean_excitation_energy > 0.0);
    ASSERT_TRUE(mat_elements->elements[1].mean_excitation_energy > 0.0);
    ASSERT_TRUE(mat_effective->elements[0].mean_excitation_energy > 0.0);
    ASSERT_TRUE(mat_effective->elements[1].mean_excitation_energy > 0.0);

    rc = osh_transport_material_prepare(wm, 1u, &tables);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(tables.nprojectiles == 1u);

    expected_elements = eval_element_sum_bethe(mat_elements, 1u, 1u, tables.emin);
    expected_effective = eval_element_sum_bethe(mat_effective, 1u, 1u, tables.emin);
    got_elements = osh_transport_material_sp_lookup(&tables, mat_elements->index, 0u, tables.emin);
    got_effective = osh_transport_material_sp_lookup(&tables, mat_effective->index, 0u, tables.emin);

    ASSERT_NEAR(got_elements, expected_elements, 5e-2);
    ASSERT_NEAR(got_effective, expected_effective, 5e-2);
    ASSERT_TRUE(fabs(got_effective - eval_effective_bethe(mat_effective, 1u, 1u, tables.emin)) > 1e-2);

    osh_transport_material_runtime_free(&tables);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(mat_path) == 0);
}

static void test_prepare_rejects_projectiles_beyond_isotope_db(void) {
    char mat_path[512];
    char mat_text[1024];
    struct material_workspace *wm;
    struct osh_transport_material_runtime tables;
    enum osh_status rc;

    wm = NULL;
    memset(&tables, 0, sizeof(tables));

    snprintf(mat_text,
             sizeof(mat_text),
             "MATERIAL Water\n"
             "ICRU 276\n");
    write_temp_file(mat_path, sizeof(mat_path), mat_text);

    rc = osh_material_setup_from_path(mat_path, NULL, &wm);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    rc = osh_transport_material_prepare(wm, (unsigned int) OSH_ISOTOPE_DB_NELEM, &tables);
    ASSERT_TRUE(rc == OSH_EINVAL);

    osh_transport_material_runtime_free(&tables);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(mat_path) == 0);
}

int main(void) {
    test_prepare_mixed_loaddedx_and_bethe_paths();
    test_prepare_bethe_infers_element_mee_from_material_mean_excitation();
    test_prepare_rejects_projectiles_beyond_isotope_db();
    return 0;
}
