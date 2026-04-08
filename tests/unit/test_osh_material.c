#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_rc.h"
#include "material/osh_material.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static int tmp_counter = 0;

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_material_test_%d.tmp", tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static double test_material_element_mass_da(struct material_element const *elem) {
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

static double test_material_bragg_mean_excitation_energy(struct material const *mat) {
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

static void test_material_fixture_icru_and_rho(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    snprintf(path, sizeof(path), "%s/test01/mat.dat", OSH_TEST_FIXTURES_DIR);

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    ASSERT_TRUE(wm->nmaterials == 3u);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_BLACKHOLE) != NULL);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_VACUUM) != NULL);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_BLACKHOLE)->name != NULL);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_VACUUM)->name != NULL);
    ASSERT_TRUE(osh_material_by_index(wm, OSH_MATERIAL_INDEX_BLACKHOLE)->rho < 0.0);
    ASSERT_TRUE(fabs(osh_material_by_index(wm, OSH_MATERIAL_INDEX_BLACKHOLE)->rgba[3] - 1.0f) < 1e-6);
    ASSERT_TRUE(fabs(osh_material_by_index(wm, OSH_MATERIAL_INDEX_VACUUM)->rgba[3]) < 1e-6);

    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(osh_material_by_name(wm, "water") == NULL);
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
    ASSERT_TRUE(mat->elements[0].mean_excitation_energy < 0.0);
    ASSERT_TRUE(mat->elements[1].z == 8u);
    ASSERT_TRUE(mat->elements[1].a == 0u);
    ASSERT_TRUE(fabs(mat->elements[1].mass_fraction - 0.888106) < 1e-7);
    ASSERT_TRUE(mat->elements[1].atom_count > 0.0);
    ASSERT_TRUE(mat->elements[1].mean_excitation_energy < 0.0);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count + mat->elements[1].atom_count - 1.0) < 1e-12);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_material_fixture_composition(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    snprintf(path, sizeof(path), "%s/tests/cases/01_simple_detect/mat.dat", OSH_PROJECT_SOURCE_DIR);

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    ASSERT_TRUE(wm->nmaterials == 3u);

    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->index == OSH_MATERIAL_INDEX_FIRST_USER);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_TRUE(fabs(mat->rho - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 78.02104963) < 1e-8);
    ASSERT_TRUE(mat->nelements == 2u);
    ASSERT_TRUE(mat->elements[0].z == 1u);
    ASSERT_TRUE(mat->elements[0].a == 0u);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count - 2.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[0].mass_fraction - 0.1118983441) < 1e-6);
    ASSERT_TRUE(fabs(mat->elements[0].mean_excitation_energy - 22.9) < 1e-12);
    ASSERT_TRUE(mat->elements[1].z == 8u);
    ASSERT_TRUE(mat->elements[1].a == 0u);
    ASSERT_TRUE(fabs(mat->elements[1].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[1].mass_fraction - 0.8881016559) < 1e-6);
    ASSERT_TRUE(fabs(mat->elements[1].mean_excitation_energy - 106.0) < 1e-12);
    ASSERT_TRUE(fabs(test_material_bragg_mean_excitation_energy(mat) - mat->mean_excitation_energy) < 1e-12);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
}

static void test_material_keys_and_state_are_case_insensitive(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MaTeRiAl Water\n"
                    "StAtE GaS\n"
                    "CoLoUr 0.1 0.2 0.3 0.4\n"
                    "rHo 0.001\n"
                    "nUcLiD 1 1\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->index == OSH_MATERIAL_INDEX_FIRST_USER);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_GAS);
    ASSERT_TRUE(fabs(mat->rgba[0] - 0.1f) < 1e-6);
    ASSERT_TRUE(fabs(mat->rgba[1] - 0.2f) < 1e-6);
    ASSERT_TRUE(fabs(mat->rgba[2] - 0.3f) < 1e-6);
    ASSERT_TRUE(fabs(mat->rgba[3] - 0.4f) < 1e-6);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_level_mean_excitation_and_dedx_table(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Water\n"
                    "MIVALUE 80.0\n"
                    "STATE gas\n"
                    "LOADDEDX tables/water_dedx.dat\n"
                    "ICRU 276\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 80.0) < 1e-12);
    ASSERT_TRUE(mat->elements[0].mean_excitation_energy < 0.0);
    ASSERT_TRUE(mat->elements[1].mean_excitation_energy < 0.0);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_GAS);
    ASSERT_TRUE(mat->dedx_table_path != NULL);
    ASSERT_TRUE(strcmp(mat->dedx_table_path, "tables/water_dedx.dat") == 0);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_level_mean_excitation_clears_compound_element_values(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Water\n"
                    "RHO 1.0\n"
                    "MATERIALI 75.0\n"
                    "ELEMENT 1 2\n"
                    "IVALUE 22.9\n"
                    "ELEMENT 8 1\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 75.0) < 1e-12);
    ASSERT_TRUE(mat->elements[0].mean_excitation_energy < 0.0);
    ASSERT_TRUE(mat->elements[1].mean_excitation_energy < 0.0);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_mass_fraction_input_completes_atom_counts(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL WaterByMass\n"
                    "RHO 1.0\n"
                    "ELEMENTBYMASS 1 0.1118983441\n"
                    "ELEMENTBYMASS 8 0.8881016559\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "WaterByMass");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->nelements == 2u);
    ASSERT_TRUE(fabs(mat->elements[0].mass_fraction - 0.1118983441) < 1e-10);
    ASSERT_TRUE(fabs(mat->elements[1].mass_fraction - 0.8881016559) < 1e-10);
    ASSERT_TRUE(mat->elements[0].atom_count > mat->elements[1].atom_count);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count + mat->elements[1].atom_count - 1.0) < 1e-12);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_explicit_isotope_uses_mass_number(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Li6\n"
                    "RHO 0.534\n"
                    "NUCLID 3 6 1.0\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "Li6");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->nelements == 1u);
    ASSERT_TRUE(mat->elements[0].z == 3u);
    ASSERT_TRUE(mat->elements[0].a == 6u);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[0].mass_fraction - 1.0) < 1e-12);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_defaults_state_to_condensed(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL DefaultState\n"
                    "RHO 1.0\n"
                    "NUCLID 1 1.0\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "DefaultState");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_CONDENSED);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_icru_element_expands_to_natural_element(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Lithium\n"
                    "ICRU 3\n"
                    "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "Lithium");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->icru_id == 3);
    ASSERT_TRUE(fabs(mat->rho - 0.534) < 1e-7);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 40.0) < 1e-12);
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_TRUE(mat->nelements == 1u);
    ASSERT_TRUE(mat->elements[0].z == 3u);
    ASSERT_TRUE(mat->elements[0].a == 0u);
    ASSERT_TRUE(fabs(mat->elements[0].mass_fraction - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat->elements[0].mean_excitation_energy - 40.0) < 1e-12);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_rejects_key_outside_material(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "RHO 1.0\n"
                    "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(wm == NULL);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_accepts_end_outside_material(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    ASSERT_TRUE(wm->nmaterials == 2u);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_rejects_end_arguments(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Water\n"
                    "ICRU 276\n"
                    "END trailing\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(wm == NULL);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_accepts_missing_end(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Water\n"
                    "ICRU 276\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->icru_id == 276);
    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_new_block_implicitly_closes_previous_one(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat_a;
    struct material const *mat_b;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL First\n"
                    "RHO 1.0\n"
                    "NUCLID 1 1.0\n"
                    "MATERIAL Second\n"
                    "RHO 2.0\n"
                    "NUCLID 8 1.0\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat_a = osh_material_by_name(wm, "First");
    mat_b = osh_material_by_name(wm, "Second");
    ASSERT_TRUE(mat_a != NULL);
    ASSERT_TRUE(mat_b != NULL);
    ASSERT_TRUE(fabs(mat_a->rho - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(mat_b->rho - 2.0) < 1e-12);
    ASSERT_TRUE(mat_a->nelements == 1u);
    ASSERT_TRUE(mat_b->nelements == 1u);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_rejects_element_i_before_element(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Water\n"
                    "RHO 1.0\n"
                    "IVALUE 22.9\n"
                    "NUCLID 1 1\n"
                    "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(wm == NULL);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_rejects_mixed_icru_and_elements(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Water\n"
                    "ICRU 276\n"
                    "RHO 1.0\n"
                    "NUCLID 1 2\n"
                    "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(wm == NULL);
    ASSERT_TRUE(remove(path) == 0);
}

static void test_material_numeric_name_is_still_a_name(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    struct material const *mat;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL 999\n"
                    "ICRU 276\n"
                    "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "999");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->index == OSH_MATERIAL_INDEX_FIRST_USER);

    ASSERT_TRUE(osh_material_workspace_free(wm) == OSH_OK);
    ASSERT_TRUE(remove(path) == 0);
}

int main(void) {
    test_material_fixture_icru_and_rho();
    test_material_fixture_composition();
    test_material_keys_and_state_are_case_insensitive();
    test_material_level_mean_excitation_and_dedx_table();
    test_material_level_mean_excitation_clears_compound_element_values();
    test_material_mass_fraction_input_completes_atom_counts();
    test_material_explicit_isotope_uses_mass_number();
    test_material_defaults_state_to_condensed();
    test_material_icru_element_expands_to_natural_element();
    test_material_rejects_key_outside_material();
    test_material_accepts_end_outside_material();
    test_material_rejects_end_arguments();
    test_material_accepts_missing_end();
    test_material_rejects_element_i_before_element();
    test_material_rejects_mixed_icru_and_elements();
    test_material_numeric_name_is_still_a_name();
    test_material_new_block_implicitly_closes_previous_one();
    return 0;
}
