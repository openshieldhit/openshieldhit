#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/simulation.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

static int tmp_counter = 0;

static void write_temp_file(char *path_out, size_t cap, char const *content) {
    FILE *fp;
    snprintf(path_out, cap, "osh_sim_test_%d.tmp", tmp_counter++);
    fp = fopen(path_out, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    fclose(fp);
}

/* ---- Tests --------------------------------------------------------------- */

static void test_create_rejects_null_args(void) {
    struct osh_simulation *sim = NULL;

    ASSERT_TRUE(osh_simulation_create(NULL, NULL, NULL, NULL, NULL) == OSH_EINVAL);
    ASSERT_TRUE(osh_simulation_create(NULL, NULL, NULL, NULL, &sim) == OSH_EINVAL);
    ASSERT_TRUE(sim == NULL);
}

static void test_create_rejects_unprepared_geometry(void) {
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    enum osh_status rc;

    ASSERT_TRUE(osh_geometry_workspace_create(&geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_workspace_create(&beam) == OSH_OK);
    ASSERT_TRUE(osh_material_workspace_create(&mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_workspace_create(&scoring) == OSH_OK);

    /* geo->prepared is NULL: must be rejected before any expensive setup */
    rc = osh_simulation_create(beam, geo, mat, scoring, &sim);
    ASSERT_TRUE(rc == OSH_ESTATE);
    ASSERT_TRUE(sim == NULL);

    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
}

static void test_create_rejects_unknown_zone_material(void) {
    char geo_path[256];
    char mat_path[256];
    char beam_path[256];
    char scoring_path[256];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    enum osh_status rc;

    /* zone 001 is assigned "NoSuchMat" which is not defined in mat.dat */
    write_temp_file(geo_path,
                    sizeof(geo_path),
                    "  RCC 1  0.0 0.0 0.0  0.0 0.0 10.0\n"
                    "         5.0\n"
                    "  END\n"
                    "  001  +1\n"
                    "  END\n"
                    "  ASSIGNMAT NoSuchMat 001\n");
    write_temp_file(mat_path,
                    sizeof(mat_path),
                    "MATERIAL Water\n"
                    "  ICRU 276\n");
    write_temp_file(beam_path,
                    sizeof(beam_path),
                    "RNDSEED    1\n"
                    "PRIMARY    1 1\n"
                    "TMAX0    150.0 0.0\n"
                    "BEAMSIGMA  0.0 0.0\n"
                    "BEAMPOS    0.0 0.0 -5.0\n"
                    "NSTAT      1 -1\n"
                    "DELTAE   0.05\n"
                    "DEMIN    0.025\n");
    write_temp_file(scoring_path,
                    sizeof(scoring_path),
                    "Geometry Mesh\n"
                    "  Name G\n"
                    "  Z 0.0 10.0 1\n"
                    "Output\n"
                    "  Filename out.dat\n"
                    "  Fileformat TEXT\n"
                    "  Geo G\n"
                    "  Quantity Energy\n");

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    rc = osh_simulation_create(beam, geo, mat, scoring, &sim);
    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(sim == NULL);

    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);

    remove(geo_path);
    remove(mat_path);
    remove(beam_path);
    remove(scoring_path);
}

static void test_create_free_lifecycle(void) {
    char geo_path[512];
    char beam_path[512];
    char mat_path[512];
    char scoring_path[512];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/00_minimal/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/00_minimal/beam.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/00_minimal/mat.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(scoring_path, sizeof(scoring_path), "%s/tests/cases/00_minimal/detect.dat", OSH_PROJECT_SOURCE_DIR);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, &sim) == OSH_OK);
    ASSERT_TRUE(sim != NULL);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_free(NULL) == OSH_OK); /* NULL is a no-op */

    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
}

/* ---- Entry point --------------------------------------------------------- */

int main(void) {
    test_create_rejects_null_args();
    test_create_rejects_unprepared_geometry();
    test_create_rejects_unknown_zone_material();
    test_create_free_lifecycle();
    return 0;
}
