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
    ASSERT_TRUE(fabs(osh_material_by_index(wm, OSH_MATERIAL_INDEX_BLACKHOLE)->rgba[3] - 1.0f) < 1e-6);
    ASSERT_TRUE(fabs(osh_material_by_index(wm, OSH_MATERIAL_INDEX_VACUUM)->rgba[3]) < 1e-6);

    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(mat->index == OSH_MATERIAL_INDEX_FIRST_USER);
    ASSERT_TRUE(mat->icru_id == 276);
    ASSERT_TRUE(fabs(mat->rho - 0.9987654321) < 1e-12);
    ASSERT_TRUE(mat->nelements == 0u);

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
    ASSERT_TRUE(mat->state == OSH_MATERIAL_STATE_GAS);
    ASSERT_TRUE(fabs(mat->rho - 1.0) < 1e-12);
    ASSERT_TRUE(mat->nelements == 2u);
    ASSERT_TRUE(mat->elements[0].z == 1u);
    ASSERT_TRUE(fabs(mat->elements[0].atom_count - 2.0) < 1e-12);
    ASSERT_TRUE(mat->elements[0].mass_fraction < 0.0);
    ASSERT_TRUE(fabs(mat->elements[0].mean_excitation_energy - 22.9) < 1e-12);
    ASSERT_TRUE(mat->elements[1].z == 8u);
    ASSERT_TRUE(fabs(mat->elements[1].atom_count - 1.0) < 1e-12);
    ASSERT_TRUE(mat->elements[1].mass_fraction < 0.0);
    ASSERT_TRUE(mat->elements[1].mean_excitation_energy < 0.0);

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
                    "nUcLiD 1 1\n"
                    "eNd\n");

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
                    "MIVALUE 75.0\n"
                    "LOADDEDX tables/water_dedx.dat\n"
                    "ICRU 276\n"
                    "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);
    mat = osh_material_by_name(wm, "Water");
    ASSERT_TRUE(mat != NULL);
    ASSERT_TRUE(fabs(mat->mean_excitation_energy - 75.0) < 1e-12);
    ASSERT_TRUE(mat->dedx_table_path != NULL);
    ASSERT_TRUE(strcmp(mat->dedx_table_path, "tables/water_dedx.dat") == 0);

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

static void test_material_rejects_end_outside_material(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), "END\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(wm == NULL);
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

static void test_material_rejects_missing_end(void) {
    char path[512];
    struct material_workspace *wm = NULL;
    enum osh_status rc;

    write_temp_file(path,
                    sizeof(path),
                    "MATERIAL Water\n"
                    "ICRU 276\n");

    rc = osh_material_setup_from_path(path, NULL, &wm);

    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(wm == NULL);
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
    test_material_rejects_key_outside_material();
    test_material_rejects_end_outside_material();
    test_material_rejects_end_arguments();
    test_material_rejects_missing_end();
    test_material_rejects_element_i_before_element();
    test_material_rejects_mixed_icru_and_elements();
    test_material_numeric_name_is_still_a_name();
    return 0;
}
