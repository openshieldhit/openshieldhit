#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
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
    struct osh_results const *results = NULL;

    ASSERT_TRUE(osh_simulation_create(NULL, NULL, NULL, NULL, NULL, NULL) == OSH_EINVAL);
    ASSERT_TRUE(osh_simulation_create(NULL, NULL, NULL, NULL, NULL, &sim) == OSH_EINVAL);
    ASSERT_TRUE(osh_simulation_get_results(NULL, NULL) == OSH_EINVAL);
    ASSERT_TRUE(osh_simulation_get_results(NULL, &results) == OSH_EINVAL);
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
    rc = osh_simulation_create(beam, geo, mat, scoring, NULL, &sim);
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

    rc = osh_simulation_create(beam, geo, mat, scoring, NULL, &sim);
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
    struct osh_results const *results = NULL;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/00_minimal/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/00_minimal/beam.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/00_minimal/mat.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(scoring_path, sizeof(scoring_path), "%s/tests/cases/00_minimal/detect.dat", OSH_PROJECT_SOURCE_DIR);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    ASSERT_TRUE(sim != NULL);
    ASSERT_TRUE(osh_simulation_get_results(sim, &results) == OSH_OK);
    ASSERT_TRUE(results != NULL);
    ASSERT_TRUE(osh_results_requested_nstat(results) == (unsigned long long) beam->nstat);
    ASSERT_TRUE(osh_results_completed_nstat(results) == 0ull);
    ASSERT_TRUE(osh_results_has_completed_run(results) == 0);
    ASSERT_TRUE(osh_simulation_save(sim) == OSH_ESTATE);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_free(NULL) == OSH_OK); /* NULL is a no-op */

    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
}

static void test_profiling_run_lifecycle(void) {
    char geo_path[512];
    char beam_path[512];
    char mat_path[512];
    char scoring_path[512];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    struct osh_simulation_profile prof;
    double phase_sum_s;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/00_minimal/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/00_minimal/beam.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/00_minimal/mat.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(scoring_path, sizeof(scoring_path), "%s/tests/cases/00_minimal/detect.dat", OSH_PROJECT_SOURCE_DIR);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    beam->nstat = 50u; /* keep the transport run short */

    ASSERT_TRUE(osh_simulation_set_profiling(NULL, 1) == OSH_EINVAL);

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);

    /* Profiling must be explicitly enabled before it can be queried. */
    ASSERT_TRUE(osh_simulation_get_profile(sim, &prof) == OSH_ESTATE);
    ASSERT_TRUE(osh_simulation_set_profiling(sim, 1) == OSH_OK);

    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);

    ASSERT_TRUE(osh_simulation_get_profile(sim, NULL) == OSH_EINVAL);
    ASSERT_TRUE(osh_simulation_get_profile(sim, &prof) == OSH_OK);
    ASSERT_TRUE(prof.steps > 0ull);
    ASSERT_TRUE(prof.iterations > 0ull);
    ASSERT_TRUE(prof.transport_s > 0.0);

    /* Phase timers must decompose the transport wall time: every phase is
     * non-negative and the sum cannot exceed the total (small slack for
     * progress reporting and clock granularity). */
    ASSERT_TRUE(prof.phase_fill_s >= 0.0);
    ASSERT_TRUE(prof.phase_zone_ref_s >= 0.0);
    ASSERT_TRUE(prof.phase_distance_s >= 0.0);
    ASSERT_TRUE(prof.phase_step_s > 0.0);
    ASSERT_TRUE(prof.phase_compact_s >= 0.0);
    phase_sum_s = prof.phase_fill_s + prof.phase_zone_ref_s + prof.phase_distance_s + prof.phase_step_s
                  + prof.phase_compact_s;
    ASSERT_TRUE(phase_sum_s <= prof.transport_s * 1.02 + 1.0e-6);

    /* Disabling profiling makes the profile unavailable again. */
    ASSERT_TRUE(osh_simulation_set_profiling(sim, 0) == OSH_OK);
    ASSERT_TRUE(osh_simulation_get_profile(sim, &prof) == OSH_ESTATE);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);

    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
}

static void test_create_rejects_voxel_geometry_without_hutable(void) {
    char beam_path[512];
    char mat_path[256];
    char scoring_path[512];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    struct osh_gemca_prepared *prepared = NULL;
    struct body *body = NULL;
    struct zone *zone = NULL;
    enum osh_status rc;

    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/00_minimal/beam.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(scoring_path, sizeof(scoring_path), "%s/tests/cases/00_minimal/detect.dat", OSH_PROJECT_SOURCE_DIR);
    write_temp_file(mat_path,
                    sizeof(mat_path),
                    "MATERIAL Water\n"
                    "  ICRU 276\n");

    ASSERT_TRUE(osh_geometry_workspace_create(&geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    prepared = (struct osh_gemca_prepared *) calloc(1, sizeof(*prepared));
    body = (struct body *) calloc(1, sizeof(*body));
    zone = (struct zone *) calloc(1, sizeof(*zone));
    ASSERT_TRUE(prepared != NULL);
    ASSERT_TRUE(body != NULL);
    ASSERT_TRUE(zone != NULL);

    prepared->bodies = (struct body **) calloc(1u, sizeof(struct body *));
    prepared->zones = (struct zone **) calloc(1u, sizeof(struct zone *));
    ASSERT_TRUE(prepared->bodies != NULL);
    ASSERT_TRUE(prepared->zones != NULL);

    body->type = OSH_GEMCA_BODY_VOX;
    zone->material_name = strdup("Water");
    ASSERT_TRUE(zone->material_name != NULL);

    prepared->bodies[0] = body;
    prepared->zones[0] = zone;
    prepared->nbodies = 1u;
    prepared->nzones = 1u;
    geo->prepared = prepared;

    rc = osh_simulation_create(beam, geo, mat, scoring, NULL, &sim);
    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(sim == NULL);

    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
    remove(mat_path);
}

/* ---- Entry point --------------------------------------------------------- */

int main(void) {
    test_create_rejects_null_args();
    test_create_rejects_unprepared_geometry();
    test_create_rejects_unknown_zone_material();
    test_create_free_lifecycle();
    test_profiling_run_lifecycle();
    test_create_rejects_voxel_geometry_without_hutable();
    return 0;
}
