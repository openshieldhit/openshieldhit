#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/material.h"
#include "openshieldhit/status.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void fixture_path(char *path, size_t cap, char const *rel) {
    snprintf(path, cap, "%s/material/%s", OSH_TEST_FIXTURES_DIR, rel);
}

static void
load_workspace_from_fixture(char const *rel, enum osh_status expected_rc, struct osh_material_workspace **wm_out) {
    char path[1024];
    struct osh_material_workspace *wm = NULL;
    enum osh_status rc;

    fixture_path(path, sizeof(path), rel);
    rc = osh_material_setup_from_path(path, NULL, &wm);
    ASSERT_TRUE(rc == expected_rc);
    if (expected_rc == OSH_OK) {
        ASSERT_TRUE(wm != NULL);
    } else {
        ASSERT_TRUE(wm == NULL);
    }
    *wm_out = wm;
}

static double test_material_element_mass_da(struct osh_material_element const *elem) {
    if (elem->a == 0u && elem->z == 1u) {
        return 1.00794;
    }
    if (elem->a == 0u && elem->z == 8u) {
        return 15.9994;
    }
    fprintf(stderr, "unsupported test element Z=%u A=%u\n", elem->z, elem->a);
    ASSERT_TRUE(0);
    return 0.0;
}

static double test_material_bragg_mean_excitation_energy(struct osh_material const *mat) {
    double numerator;
    double denominator;
    double mass;
    double weight;
    size_t i;

    numerator = 0.0;
    denominator = 0.0;

    i = 0;
    while (i < mat->nelements) {
        mass = test_material_element_mass_da(&mat->elements[i]);
        weight = mat->elements[i].mass_fraction * (double) mat->elements[i].z / mass;
        numerator += weight * log(mat->elements[i].mean_excitation_energy);
        denominator += weight;
        i++;
    }

    ASSERT_TRUE(denominator > 0.0);
    return exp(numerator / denominator);
}

static void test_valid_icru_override(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material const *mat;

    load_workspace_from_fixture("valid/water_icru_override.dat", OSH_OK, &wm);

    ASSERT_TRUE(wm->nmaterials == 3u);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_BLACKHOLE) != NULL);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_VACUUM) != NULL);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_BLACKHOLE)->name != NULL);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_VACUUM)->name != NULL);

    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->index == OSH_MATERIAL_INDEX_FIRST_USER);
    ASSERT_TRUE(mat->icru_id == 276);
    ASSERT_TRUE(fabs(mat->rho - 0.9987654321) < 1e-12);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 75.0) < 1e-12);
    ASSERT_TRUE(mat->nelements == 2u);
    ASSERT_TRUE(mat->elements[0].z == 1u);
    ASSERT_TRUE(mat->elements[0].a == 0u);
    ASSERT_TRUE(fabs(mat->elements[0].mass_fraction - 0.111894) < 1e-7);
    ASSERT_TRUE(mat->elements[0].atom_count > 0.0);
    ASSERT_TRUE(mat->elements[1].z == 8u);
    ASSERT_TRUE(mat->elements[1].a == 0u);
    ASSERT_TRUE(fabs(mat->elements[1].mass_fraction - 0.888106) < 1e-7);
    ASSERT_TRUE(mat->elements[1].atom_count > 0.0);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count + mat->elements[1].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(test_material_bragg_mean_excitation_energy(mat) - mat->mean_excitation_energy) < 1e-12);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_valid_explicit_composition(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material const *mat;

    load_workspace_from_fixture("valid/water_explicit.dat", OSH_OK, &wm);

    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_TRUE(fabs(mat->rho - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 78.02104963) < 1e-8);
    ASSERT_TRUE(mat->nelements == 2u);
    ASSERT_TRUE(mat->elements[0].z == 1u);
    ASSERT_TRUE(mat->elements[1].z == 8u);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count - 2.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[1].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[0].mean_excitation_energy - 22.9) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[1].mean_excitation_energy - 106.0) < 1e-12);
    ASSERT_TRUE(fabs(test_material_bragg_mean_excitation_energy(mat) - mat->mean_excitation_energy) < 1e-12);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_valid_case_and_loaddedx(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material const *mat;

    load_workspace_from_fixture("valid/case_and_loaddedx.dat", OSH_OK, &wm);

    mat = osh_material_by_name(wm, "WaterFancy");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_GAS);
    ASSERT_TRUE(fabs(mat->rgba[0] - 0.1f) < 1e-6);
    ASSERT_TRUE(fabs(mat->rgba[1] - 0.2f) < 1e-6);
    ASSERT_TRUE(fabs(mat->rgba[2] - 0.3f) < 1e-6);
    ASSERT_TRUE(fabs(mat->rgba[3] - 0.4f) < 1e-6);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 80.0) < 1e-12);
    ASSERT_TRUE(mat->ndedx_overrides == 18u);
    ASSERT_TRUE(mat->dedx_overrides != NULL);
    ASSERT_TRUE(mat->dedx_overrides[0].projectile_z == 1u);
    ASSERT_TRUE(mat->dedx_overrides[0].npoints == 2u);
    ASSERT_TRUE(fabs(mat->dedx_overrides[0].energy_mev_per_u[0] - 0.025) < 1e-12);
    ASSERT_TRUE(fabs(mat->dedx_overrides[0].dedx_mev_cm2_per_g[0] - 100.0) < 1e-12);
    ASSERT_TRUE(mat->dedx_overrides[17].projectile_z == 18u);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_material_dedx_set_replace_and_clear(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material *mat;
    double e1[] = {0.025, 1.0, 10.0};
    double sp1[] = {100.0, 10.0, 1.0};
    double e2[] = {0.025, 2.0, 20.0};
    double sp2[] = {200.0, 20.0, 2.0};

    ASSERT_TRUE(osh_material_workspace_create(&wm) == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    mat = &wm->materials[OSH_MATERIAL_INDEX_VACUUM];
    ASSERT_TRUE(osh_material_dedx_set(mat, 6u, e1, sp1, 3u) == OSH_OK);
    ASSERT_TRUE(mat->ndedx_overrides == 1u);
    ASSERT_TRUE(mat->dedx_overrides[0].projectile_z == 6u);
    ASSERT_TRUE(fabs(mat->dedx_overrides[0].dedx_mev_cm2_per_g[1] - 10.0) < 1e-12);

    ASSERT_TRUE(osh_material_dedx_set(mat, 6u, e2, sp2, 3u) == OSH_OK);
    ASSERT_TRUE(mat->ndedx_overrides == 1u);
    ASSERT_TRUE(fabs(mat->dedx_overrides[0].energy_mev_per_u[1] - 2.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->dedx_overrides[0].dedx_mev_cm2_per_g[0] - 200.0) < 1e-12);

    osh_material_dedx_clear(mat);
    ASSERT_TRUE(mat->ndedx_overrides == 0u);
    ASSERT_TRUE(mat->dedx_overrides == NULL);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_valid_mass_fraction_and_isotope(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material const *mat;

    load_workspace_from_fixture("valid/water_by_mass.dat", OSH_OK, &wm);
    mat = osh_material_by_name(wm, "WaterByMass");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->nelements == 2u);
    ASSERT_TRUE(fabs(mat->elements[0].mass_fraction - 0.1118983441) < 1e-10);
    ASSERT_TRUE(fabs(mat->elements[1].mass_fraction - 0.8881016559) < 1e-10);
    ASSERT_TRUE(mat->elements[0].atom_count > mat->elements[1].atom_count);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count + mat->elements[1].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);

    load_workspace_from_fixture("valid/li6_isotope.dat", OSH_OK, &wm);
    mat = osh_material_by_name(wm, "Li6");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->nelements == 1u);
    ASSERT_TRUE(mat->elements[0].z == 3u);
    ASSERT_TRUE(mat->elements[0].a == 6u);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[0].mass_fraction - 1.0) < 1e-12);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_valid_end_handling(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material const *mat_a;
    struct osh_material const *mat_b;

    load_workspace_from_fixture("valid/end_outside.dat", OSH_OK, &wm);
    ASSERT_TRUE(wm->nmaterials == 2u);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);

    load_workspace_from_fixture("valid/missing_end.dat", OSH_OK, &wm);
    ASSERT_TRUE(osh_material_by_name(wm, "Water") != NULL);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);

    load_workspace_from_fixture("valid/two_blocks_no_end.dat", OSH_OK, &wm);
    mat_a = osh_material_by_name(wm, "First");
    mat_b = osh_material_by_name(wm, "Second");
    ASSERT_TRUE(mat_a != NULL);
    ASSERT_TRUE(mat_b != NULL);
    ASSERT_TRUE(fabs(mat_a->rho - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat_b->rho - 2.0) < 1e-12);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_invalid_key_outside_material(void) {
    struct osh_material_workspace *wm = NULL;
    load_workspace_from_fixture("invalid/key_outside_material.dat", OSH_EPARSE, &wm);
}

static void test_invalid_conflicting_mee(void) {
    struct osh_material_workspace *wm = NULL;
    load_workspace_from_fixture("invalid/conflicting_mee.dat", OSH_EPARSE, &wm);
}

static void test_invalid_end_arguments(void) {
    struct osh_material_workspace *wm = NULL;
    load_workspace_from_fixture("invalid/end_with_arguments.dat", OSH_EPARSE, &wm);
}

static void test_invalid_elementi_before_element(void) {
    struct osh_material_workspace *wm = NULL;
    load_workspace_from_fixture("invalid/elementi_before_element.dat", OSH_EPARSE, &wm);
}

static void test_invalid_mixed_icru_and_elements(void) {
    struct osh_material_workspace *wm = NULL;
    load_workspace_from_fixture("invalid/mixed_icru_and_elements.dat", OSH_EPARSE, &wm);
}

int main(void) {
    test_valid_icru_override();
    test_valid_explicit_composition();
    test_valid_case_and_loaddedx();
    test_material_dedx_set_replace_and_clear();
    test_valid_mass_fraction_and_isotope();
    test_valid_end_handling();

    test_invalid_key_outside_material();
    test_invalid_conflicting_mee();
    test_invalid_end_arguments();
    test_invalid_elementi_before_element();
    test_invalid_mixed_icru_and_elements();
    return 0;
}
