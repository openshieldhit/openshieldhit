#include <math.h>
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
    phase_sum_s =
        prof.phase_fill_s + prof.phase_zone_ref_s + prof.phase_distance_s + prof.phase_step_s + prof.phase_compact_s;
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

/* Read the "# PRIMARIES: <n>" header value from an ASCII scoring output. */
static unsigned long long read_primaries_header(char const *path) {
    FILE *fp;
    char line[512];
    unsigned long long primaries = 0ull;
    int found = 0;

    fp = fopen(path, "r");
    ASSERT_TRUE(fp != NULL);
    while (fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "# PRIMARIES: %llu", &primaries) == 1) {
            found = 1;
            break;
        }
    }
    fclose(fp);
    ASSERT_TRUE(found == 1);
    return primaries;
}

/* Stop callback that always asks to stop (a constant-true should_stop). */
static int always_stop(void *user) {
    (void) user;
    return 1;
}

/*
 * Clean stop (issue #192 / #195): a should_stop callback that returns true
 * halts new-primary injection at the first safe point.  In-flight histories
 * drain, so the completed count is exactly one pool batch — deterministic, no
 * timing dependence — which makes this assertion robust on every platform/CI OS.
 * The partial result must save correctly, normalised by the *true* completed
 * count, and the ASCII "# PRIMARIES:" header must echo that same count.
 */
static void test_clean_stop_partial_result_is_exact(void) {
    char geo_path[512];
    char beam_path[512];
    char mat_path[512];
    char scoring_path[512];
    char scoring_text[1024];
    char out_path[256];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    struct osh_results const *results = NULL;
    size_t const pool_cap = 100u;
    unsigned long long const requested = 20000ull;
    unsigned long long completed;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/00_minimal/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/00_minimal/beam.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/00_minimal/mat.dat", OSH_PROJECT_SOURCE_DIR);

    /* A single TEXT output to a unique per-process filename, so there is exactly
     * one ASCII file to read the "# PRIMARIES:" header back from (the shared
     * 00_minimal detect.dat also emits a BDO, which would overwrite it). */
    snprintf(out_path, sizeof(out_path), "osh_clean_stop_%d.dat", tmp_counter++);
    snprintf(scoring_text,
             sizeof(scoring_text),
             "Geometry Mesh\n"
             "  Name M\n"
             "  X -5.0 5.0 1\n"
             "  Y -5.0 5.0 1\n"
             "  Z 0.0 20.0 200\n"
             "Output\n"
             "  Filename %s\n"
             "  Fileformat TEXT\n"
             "  Geo M\n"
             "  Quantity Energy\n",
             out_path);
    write_temp_file(scoring_path, sizeof(scoring_path), scoring_text);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    /* Request far more primaries than one pool batch so the stop is observable. */
    beam->nstat = (size_t) requested;

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);

    /* Small batch + always-stop callback → exactly one batch of primaries completes. */
    ASSERT_TRUE(osh_simulation_set_pool_capacity(sim, pool_cap) == OSH_OK);
    ASSERT_TRUE(osh_simulation_set_run_control(sim, 0.0, always_stop, NULL) == OSH_OK);

    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);

    ASSERT_TRUE(osh_simulation_get_results(sim, &results) == OSH_OK);
    completed = osh_results_completed_nstat(results);

    /* Exact: one pool batch finished, far short of the request. */
    ASSERT_TRUE(osh_results_requested_nstat(results) == requested);
    ASSERT_TRUE(completed == (unsigned long long) pool_cap);
    ASSERT_TRUE(completed < requested);
    ASSERT_TRUE(osh_results_has_completed_run(results) == 1);

    /* The partial result saves and its header echoes the true completed count. */
    ASSERT_TRUE(osh_simulation_save(sim) == OSH_OK);
    ASSERT_TRUE(read_primaries_header(out_path) == completed);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
    remove(out_path);
    remove(scoring_path);
}

/*
 * Run control left off (NULL callback, zero budget) must be bit-for-bit the old
 * behaviour: every requested primary completes.
 */
static void test_no_run_control_runs_to_completion(void) {
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

    beam->nstat = 50u;

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    /* A disabled policy (no budget, no callback) is a no-op: leaves the fast path. */
    ASSERT_TRUE(osh_simulation_set_run_control(sim, 0.0, NULL, NULL) == OSH_OK);
    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);

    ASSERT_TRUE(osh_simulation_get_results(sim, &results) == OSH_OK);
    ASSERT_TRUE(osh_results_completed_nstat(results) == 50ull);
    ASSERT_TRUE(osh_results_requested_nstat(results) == 50ull);

    ASSERT_TRUE(osh_simulation_set_run_control(NULL, 1.0, NULL, NULL) == OSH_EINVAL);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
}

/* Sum the last whitespace-separated column of every data (non-'#') line in an
 * ASCII scoring output.  For a single-quantity mesh output each data line is
 * "x y z value", so the last token is the per-primary scored value; summing them
 * gives an order-independent-ish total that a dropped batch would visibly change. */
static double sum_last_column(char const *path) {
    FILE *fp;
    char line[4096];
    double sum = 0.0;

    fp = fopen(path, "r");
    ASSERT_TRUE(fp != NULL);
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        char *tok;
        char *last = NULL;
        while (*p == ' ' || *p == '\t') {
            ++p;
        }
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') {
            continue;
        }
        for (tok = strtok(p, " \t\r\n"); tok != NULL; tok = strtok(NULL, " \t\r\n")) {
            last = tok;
        }
        if (last) {
            sum += strtod(last, NULL);
        }
    }
    fclose(fp);
    return sum;
}

/*
 * Drive one run of test case @p case_name at a given checkpoint cadence and
 * return the completed-primary count and total scored energy.  every_primaries
 * == 0 is FINAL-ONLY (one batch); > 0 runs family-complete batches of that size.
 * @p pool_capacity == 0 leaves the compiled default; > 0 sets the live-history
 * pool capacity (the perf knob) so callers can assert capacity-invariance.
 * When @p profile_out is non-NULL, profiling is enabled and the run profile is
 * copied out so callers can assert its counters accumulate across batches (and
 * read back ion_secondaries_dropped).  All bundled cases used here are water
 * cylinders spanning z ∈ [0, 20], so the same single-quantity Energy mesh reads
 * back the deposited energy for each.
 */
static void run_checkpoint_case(char const *case_name,
                                unsigned long long nstat,
                                unsigned long long every_primaries,
                                unsigned long long pool_capacity,
                                unsigned long long *completed_out,
                                double *energy_sum_out,
                                struct osh_simulation_profile *profile_out) {
    char geo_path[512];
    char beam_path[512];
    char mat_path[512];
    char scoring_path[512];
    char scoring_text[1024];
    char out_path[256];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    struct osh_results const *results = NULL;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/%s/geo.dat", OSH_PROJECT_SOURCE_DIR, case_name);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/%s/beam.dat", OSH_PROJECT_SOURCE_DIR, case_name);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/%s/mat.dat", OSH_PROJECT_SOURCE_DIR, case_name);

    /* One single-quantity TEXT mesh to a unique file, so the last column is the
     * scored energy and there is exactly one ASCII artifact to read back. */
    snprintf(out_path, sizeof(out_path), "osh_ckpt_%d.dat", tmp_counter++);
    snprintf(scoring_text,
             sizeof(scoring_text),
             "Geometry Mesh\n"
             "  Name M\n"
             "  X -5.0 5.0 1\n"
             "  Y -5.0 5.0 1\n"
             "  Z 0.0 20.0 200\n"
             "Output\n"
             "  Filename %s\n"
             "  Fileformat TEXT\n"
             "  Geo M\n"
             "  Quantity Energy\n",
             out_path);
    write_temp_file(scoring_path, sizeof(scoring_path), scoring_text);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    beam->nstat = (size_t) nstat; /* set before create so pools size to it */

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_set_checkpoint_policy(sim, every_primaries) == OSH_OK);
    if (pool_capacity != 0ull) {
        ASSERT_TRUE(osh_simulation_set_pool_capacity(sim, (size_t) pool_capacity) == OSH_OK);
    }
    if (profile_out) {
        ASSERT_TRUE(osh_simulation_set_profiling(sim, 1) == OSH_OK);
    }
    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);

    ASSERT_TRUE(osh_simulation_get_results(sim, &results) == OSH_OK);
    *completed_out = osh_results_completed_nstat(results);

    if (profile_out) {
        ASSERT_TRUE(osh_simulation_get_profile(sim, profile_out) == OSH_OK);
    }

    ASSERT_TRUE(osh_simulation_save(sim) == OSH_OK);
    *energy_sum_out = sum_last_column(out_path);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
    remove(out_path);
    remove(scoring_path);
}

/*
 * Like run_checkpoint_case, but drives periodic partial-result dumps (issue #193)
 * via osh_simulation_set_dump_control() instead of a bare checkpoint policy.  The
 * count cadence @p dump_every_primaries makes the run overwrite its output file
 * with the exact partial result at every checkpoint; the file read back at the
 * end therefore reflects the *final* save (the last dump is overwritten by it).
 * The point of the comparison in the caller is that those interleaved dumps —
 * non-destructive shadow snapshots — must not perturb the live accumulators, so
 * the final result is bit-identical to the same cadence run *without* dumps.
 */
static void run_dump_case(char const *case_name,
                          unsigned long long nstat,
                          unsigned long long dump_every_primaries,
                          unsigned long long *completed_out,
                          double *energy_sum_out) {
    char geo_path[512];
    char beam_path[512];
    char mat_path[512];
    char scoring_path[512];
    char scoring_text[1024];
    char out_path[256];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    struct osh_results const *results = NULL;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/%s/geo.dat", OSH_PROJECT_SOURCE_DIR, case_name);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/%s/beam.dat", OSH_PROJECT_SOURCE_DIR, case_name);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/%s/mat.dat", OSH_PROJECT_SOURCE_DIR, case_name);

    snprintf(out_path, sizeof(out_path), "osh_dump_%d.dat", tmp_counter++);
    snprintf(scoring_text,
             sizeof(scoring_text),
             "Geometry Mesh\n"
             "  Name M\n"
             "  X -5.0 5.0 1\n"
             "  Y -5.0 5.0 1\n"
             "  Z 0.0 20.0 200\n"
             "Output\n"
             "  Filename %s\n"
             "  Fileformat TEXT\n"
             "  Geo M\n"
             "  Quantity Energy\n",
             out_path);
    write_temp_file(scoring_path, sizeof(scoring_path), scoring_text);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    beam->nstat = (size_t) nstat;

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    /* Count cadence + no on-demand callback: deterministic periodic dumps. */
    ASSERT_TRUE(osh_simulation_set_dump_control(sim, 0.0, dump_every_primaries, NULL, NULL) == OSH_OK);
    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);

    ASSERT_TRUE(osh_simulation_get_results(sim, &results) == OSH_OK);
    *completed_out = osh_results_completed_nstat(results);

    ASSERT_TRUE(osh_simulation_save(sim) == OSH_OK);
    *energy_sum_out = sum_last_column(out_path);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
    remove(out_path);
    remove(scoring_path);
}

/*
 * A dump is a preview, never the run's product: a failed dump write must be
 * fail-soft — logged and skipped — so the run still completes.  Force every dump
 * to fail by pointing the scoring output at a non-existent subdirectory (fopen
 * cannot create the missing directory), then assert the run itself still succeeds
 * despite the periodic dumps at each checkpoint failing to write.
 */
static void test_dump_write_failure_is_fail_soft(void) {
    char geo_path[512];
    char beam_path[512];
    char mat_path[512];
    char scoring_path[512];
    char scoring_text[1024];
    char bad_out[256];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/00_minimal/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/00_minimal/beam.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/00_minimal/mat.dat", OSH_PROJECT_SOURCE_DIR);

    /* Output under a directory that does not exist → every fopen("wb") fails. */
    snprintf(bad_out, sizeof(bad_out), "osh_no_such_dir_%d/dose.dat", tmp_counter++);
    snprintf(scoring_text,
             sizeof(scoring_text),
             "Geometry Mesh\n"
             "  Name M\n"
             "  X -5.0 5.0 1\n"
             "  Y -5.0 5.0 1\n"
             "  Z 0.0 20.0 200\n"
             "Output\n"
             "  Filename %s\n"
             "  Fileformat TEXT\n"
             "  Geo M\n"
             "  Quantity Energy\n",
             bad_out);
    write_temp_file(scoring_path, sizeof(scoring_path), scoring_text);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    beam->nstat = 200u;

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    /* Cadence 40 over 200 → dumps attempted at 40/80/120/160, all failing. */
    ASSERT_TRUE(osh_simulation_set_dump_control(sim, 0.0, 40ull, NULL, NULL) == OSH_OK);
    /* The run must complete despite the failing dumps (fail-soft). */
    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
    remove(scoring_path);
}

/*
 * Periodic dumps (issue #193) must be non-destructive: interleaving mid-run
 * snapshots at a checkpoint cadence must leave the live accumulators untouched,
 * so the final result is *bit-identical* to the same-cadence run that takes no
 * dumps.  Comparing at a fixed cadence (not against final-only) isolates the dump
 * machinery from the FP-reduction-order difference that batching itself
 * introduces: the only variable between the two runs here is whether snapshots
 * are taken, so any discrepancy would be a snapshot perturbing live state.
 */
static void test_dump_control_is_non_destructive(void) {
    unsigned long long completed_batch = 0ull;
    unsigned long long completed_dump = 0ull;
    double energy_batch = 0.0;
    double energy_dump = 0.0;

    ASSERT_TRUE(osh_simulation_set_dump_control(NULL, 0.0, 0ull, NULL, NULL) == OSH_EINVAL);

    /* Same cadence (K = 40 over 200 → checkpoints at 40,80,120,160), one run with
     * periodic dumps and one without. */
    run_checkpoint_case("00_minimal", 200ull, 40ull, 0ull, &completed_batch, &energy_batch, NULL);
    run_dump_case("00_minimal", 200ull, 40ull, &completed_dump, &energy_dump);

    ASSERT_TRUE(completed_batch == 200ull);
    ASSERT_TRUE(completed_dump == 200ull);
    ASSERT_TRUE(energy_batch > 0.0);
    /* Bit-identical: the dumps changed nothing in the live run. */
    ASSERT_TRUE(energy_dump == energy_batch);
}

/*
 * Checkpoint policy (issue #195), default path: every_primaries == 0 is
 * FINAL-ONLY — the run must complete every requested primary in one batch, and
 * calling the setter with 0 must be indistinguishable from never calling it.
 */
static void test_checkpoint_final_only_default_completes(void) {
    unsigned long long completed = 0ull;
    double energy = 0.0;

    ASSERT_TRUE(osh_simulation_set_checkpoint_policy(NULL, 0ull) == OSH_EINVAL);

    run_checkpoint_case("00_minimal", 200ull, 0ull, 0ull, &completed, &energy, NULL);
    ASSERT_TRUE(completed == 200ull);
    ASSERT_TRUE(energy > 0.0);
}

/*
 * Checkpoint policy (issue #195), LIVE batching: a count cadence smaller than
 * nstat runs the whole request as several family-complete batches.  Every
 * primary must still finish (completed == nstat) and — because each batch drains
 * its families and no deposit is dropped at a checkpoint — the total scored
 * energy must match the single-batch result up to floating-point reduction order.
 * (00_minimal has NUCRE off, so the deposit multiset is identical; only the
 * per-voxel summation order differs between the two batchings.)
 */
static void test_checkpoint_live_batches_match_final_only(void) {
    unsigned long long completed_final = 0ull;
    unsigned long long completed_live = 0ull;
    double energy_final = 0.0;
    double energy_live = 0.0;
    double rel_diff;

    run_checkpoint_case("00_minimal", 200ull, 0ull, 0ull, &completed_final, &energy_final, NULL); /* one batch */
    run_checkpoint_case("00_minimal", 200ull, 64ull, 0ull, &completed_live, &energy_live, NULL);  /* 64,64,64,8 */

    ASSERT_TRUE(completed_final == 200ull);
    ASSERT_TRUE(completed_live == 200ull);
    ASSERT_TRUE(energy_final > 0.0);

    rel_diff = fabs(energy_live - energy_final) / energy_final;
    ASSERT_TRUE(rel_diff < 1.0e-6);
}

/*
 * Checkpoint policy (issue #195), LIVE batching WITH secondaries in flight.
 * 00_minimal has NUCRE off, so its equivalence test never actually drains a
 * neutron family across a checkpoint.  The idd_water_200mev_nucre1 reference deck
 * (NUCRE 1, 200 MeV p in water) enables nuclear reactions,
 * so each batch's ion pass banks neutrons that the family scheduler must fully
 * drain into scoring before the checkpoint boundary, and every nuclear event
 * appends charged secondaries (recoils, abrasion nucleons, break-up fragments)
 * to the ion pool.  The profile counters below assert secondaries really were
 * produced and drained (nuclear_events and neutrons_banked > 0), so this test
 * genuinely exercises the "family-complete, quiescent checkpoint" path rather
 * than a degenerate no-secondary run.
 *
 * The invariant (issue #213): scored energy is batch-size invariant.
 *   Final-only and multi-batch LIVE runs transport the *identical* set of
 *   histories — primaries and every nuclear secondary — because (a) the ion pool
 *   is sized with headroom beyond the primary wavefront, so no secondary is
 *   dropped (ion_secondaries_dropped == 0 below), and (b) each secondary's RNG
 *   stream is keyed by its lineage and ordinal, never drawn from the parent, so a
 *   secondary's presence or absence cannot shift its parent's stream.  The only
 *   remaining difference between the two batchings is the order deposits land in
 *   the shared scoring accumulators, i.e. floating-point summation reassociation,
 *   which is bounded far below any physical difference.  (Before #213 this fixture
 *   differed by ~2.6 %: a full primary wavefront filled the pool, its secondaries
 *   were silently dropped, and the drop count tracked the batch schedule — the
 *   "transport-kernel FP limitation" once documented here was that bug, not a
 *   kernel property.)
 */
static void test_checkpoint_live_batches_drain_families_with_nuclear(void) {
    unsigned long long completed_final = 0ull;
    unsigned long long completed_live_a = 0ull;
    unsigned long long completed_live_b = 0ull;
    double energy_final = 0.0;
    double energy_live_a = 0.0;
    double energy_live_b = 0.0;
    double rel_diff;
    struct osh_simulation_profile prof_final;
    struct osh_simulation_profile prof_live;

    /* Final-only baseline, with profiling so we can prove secondaries were in
     * flight (otherwise the fixture could silently regress to a no-nuclear run
     * and this test would prove nothing) and that none were dropped. */
    run_checkpoint_case(
        "../reference/idd_water_200mev_nucre1", 200ull, 0ull, 0ull, &completed_final, &energy_final, &prof_final);
    ASSERT_TRUE(completed_final == 200ull);
    ASSERT_TRUE(energy_final > 0.0);
    ASSERT_TRUE(prof_final.nuclear_events > 0ull);
    ASSERT_TRUE(prof_final.neutrons_banked > 0ull);
    /* Headroom fix (issue #213): the default configuration drops no ion secondary. */
    ASSERT_TRUE(prof_final.ion_secondaries_dropped == 0ull);

    /* LIVE batching drains a neutron family at every checkpoint, yet still
     * finishes every primary, still produces (and drains) neutrons, and still
     * drops nothing. */
    run_checkpoint_case(
        "../reference/idd_water_200mev_nucre1", 200ull, 64ull, 0ull, &completed_live_a, &energy_live_a, &prof_live);
    ASSERT_TRUE(completed_live_a == 200ull);
    ASSERT_TRUE(energy_live_a > 0.0);
    ASSERT_TRUE(prof_live.neutrons_banked > 0ull);
    ASSERT_TRUE(prof_live.ion_secondaries_dropped == 0ull);

    /* The real invariant: with no secondary dropped, final-only and LIVE agree to
     * floating-point summation order — NOT the ~2.6 % the pre-#213 code showed. */
    rel_diff = fabs(energy_live_a - energy_final) / energy_final;
    ASSERT_TRUE(rel_diff < 1.0e-9);

    /* Reproducibility of a fixed cadence: the same LIVE policy re-run reproduces
     * the scored energy exactly (bit-for-bit), the determinism guarantee the
     * count cadence is meant to provide. */
    run_checkpoint_case(
        "../reference/idd_water_200mev_nucre1", 200ull, 64ull, 0ull, &completed_live_b, &energy_live_b, NULL);
    ASSERT_TRUE(completed_live_b == 200ull);
    ASSERT_TRUE(energy_live_b == energy_live_a);
}

/*
 * Capacity/batch invariance with nuclear reactions ON (issue #213 regression
 * lock).  Scored energy for a fixed (seed, rndoffset) must not depend on:
 *   - the checkpoint batch size K (full / 64 / 1), nor
 *   - the live-history pool capacity (the cache/parallelism perf knob),
 * because every history's RNG stream is a pure function of its lineage and no
 * secondary is dropped (the counter is asserted 0 for each configuration).  All
 * runs therefore transport the identical set of histories; only the scorer's
 * floating-point summation order differs, which is bounded far below 1e-9
 * relative.  This is the guard the earlier equivalence test could not provide —
 * a bug whose signature was "the answer changes with K" would fail here.
 */
static void test_checkpoint_nuclear_invariant_across_batches_and_capacity(void) {
    unsigned long long const nstat = 200ull;

    /* K = 0 (final-only), 64, and 1 at the default capacity; then a small
     * explicit capacity at final-only.  Each row is (every_primaries, pool_cap). */
    struct {
        unsigned long long every_primaries;
        unsigned long long pool_capacity;
    } const configs[] = {
        {0ull, 0ull},   /* final-only, default capacity — the reference */
        {64ull, 0ull},  /* LIVE, batch 64 */
        {1ull, 0ull},   /* LIVE, batch 1 (checkpoint after every primary) */
        {0ull, 64ull},  /* final-only, capacity 64 (a different wavefront width) */
        {0ull, 256ull}, /* final-only, capacity 256 */
    };

    size_t const nconfigs = sizeof(configs) / sizeof(configs[0]);
    double energy_ref = 0.0;
    size_t i;

    for (i = 0u; i < nconfigs; ++i) {
        unsigned long long completed = 0ull;
        double energy = 0.0;
        struct osh_simulation_profile prof;

        run_checkpoint_case("../reference/idd_water_200mev_nucre1",
                            nstat,
                            configs[i].every_primaries,
                            configs[i].pool_capacity,
                            &completed,
                            &energy,
                            &prof);

        ASSERT_TRUE(completed == nstat);
        ASSERT_TRUE(energy > 0.0);
        /* Headroom holds at every tested capacity: no deposited energy is lost. */
        ASSERT_TRUE(prof.ion_secondaries_dropped == 0ull);

        if (i == 0u) {
            energy_ref = energy;
            ASSERT_TRUE(prof.nuclear_events > 0ull); /* the fixture really is nuclear */
        } else {
            ASSERT_TRUE(fabs(energy - energy_ref) / energy_ref < 1.0e-9);
        }
    }
}

/*
 * Overflow accounting (issue #213), the positive counterpart to the
 * ion_secondaries_dropped == 0 assertions above.  Those prove the default
 * configuration reserves enough headroom to drop nothing; this deliberately
 * starves the ion pool so a full wavefront's secondaries cannot all fit, and
 * asserts the overflow is COUNTED rather than silently discarded (the exact bug
 * #213 fixes).  It further checks that a drop never strands a primary (all
 * histories still complete), that the lost secondaries lower the scored energy
 * (their deposits are gone, not merely reordered), and that the starved run is
 * still bit-for-bit deterministic on re-run.  Guards the drop counter in
 * osh_transport_ion_step and the WARN in osh_simulation_run that surfaces it.
 *
 * Note this does NOT assert cross-batch invariance: once the pool overflows,
 * which secondaries are dropped depends on pool occupancy (and thus on the
 * checkpoint cadence), so the overflow regime trades exact invariance for a
 * counted, warned energy loss.  Batch/capacity invariance is a property of the
 * default, drop-free configuration and is locked by the test above.
 */
static void test_checkpoint_nuclear_overflow_is_counted_not_silent(void) {
    unsigned long long completed_ref = 0ull;
    unsigned long long completed_small = 0ull;
    unsigned long long completed_rerun = 0ull;
    double energy_ref = 0.0;
    double energy_small = 0.0;
    double energy_rerun = 0.0;
    struct osh_simulation_profile prof_ref;
    struct osh_simulation_profile prof_small;
    struct osh_simulation_profile prof_rerun;

    /* Reference: the default capacity has full secondary headroom, so nothing is
     * dropped and every nuclear deposit is scored. */
    run_checkpoint_case(
        "../reference/idd_water_200mev_nucre1", 200ull, 0ull, 0ull, &completed_ref, &energy_ref, &prof_ref);
    ASSERT_TRUE(completed_ref == 200ull);
    ASSERT_TRUE(prof_ref.nuclear_events > 0ull); /* the fixture really is nuclear */
    ASSERT_TRUE(prof_ref.ion_secondaries_dropped == 0ull);

    /* Starve the pool: capacity 2 leaves a full primary wavefront almost no room
     * for its recoils/abrasion nucleons/break-up fragments, so they overflow.
     * The overflow is counted in ion_secondaries_dropped, never a bare drop. */
    run_checkpoint_case(
        "../reference/idd_water_200mev_nucre1", 200ull, 0ull, 2ull, &completed_small, &energy_small, &prof_small);
    ASSERT_TRUE(prof_small.ion_secondaries_dropped > 0ull);
    /* A drop must never strand a primary: every history still completes. */
    ASSERT_TRUE(completed_small == 200ull);
    /* The dropped secondaries carried deposited energy that is now simply gone,
     * so the starved run scores strictly less than the headroom reference. */
    ASSERT_TRUE(energy_small > 0.0);
    ASSERT_TRUE(energy_small < energy_ref);

    /* Determinism survives overflow: re-running the identical (cadence, capacity)
     * reproduces the scored energy and the drop count exactly, because the drop
     * pattern is a pure function of the lineage-keyed histories and the fixed
     * pool occupancy — not of run-to-run chance. */
    run_checkpoint_case(
        "../reference/idd_water_200mev_nucre1", 200ull, 0ull, 2ull, &completed_rerun, &energy_rerun, &prof_rerun);
    ASSERT_TRUE(completed_rerun == 200ull);
    ASSERT_TRUE(energy_rerun == energy_small);
    ASSERT_TRUE(prof_rerun.ion_secondaries_dropped == prof_small.ion_secondaries_dropped);
}

/*
 * Checkpoint policy (issue #195), profiling accounting across batches.  The run
 * profile is a single master zeroed once at profiling-enable, and the batch loop
 * transports each checkpoint batch into it in turn.  Every counter must therefore
 * accumulate across batches: a regression that *assigns* (rather than sums) wall
 * time or step count would leave those two reflecting only the last batch, which
 * both undercounts the run and breaks the phase decomposition (phase_step_s could
 * exceed transport_s).  Because a count cadence is order-independent, the total
 * step count is a pure function of the primaries transported, so LIVE batching
 * must reproduce the single-batch step count exactly.
 */
static void test_checkpoint_profiling_accumulates_across_batches(void) {
    unsigned long long completed_final = 0ull;
    unsigned long long completed_live = 0ull;
    double energy_final = 0.0;
    double energy_live = 0.0;
    struct osh_simulation_profile prof_final;
    struct osh_simulation_profile prof_live;
    double phase_sum_s;

    run_checkpoint_case("00_minimal", 200ull, 0ull, 0ull, &completed_final, &energy_final, &prof_final);
    run_checkpoint_case("00_minimal", 200ull, 64ull, 0ull, &completed_live, &energy_live, &prof_live);

    ASSERT_TRUE(completed_final == 200ull);
    ASSERT_TRUE(completed_live == 200ull);

    /* Steps are order-independent under a count cadence: LIVE batching transports
     * exactly the same histories as the single batch, so the accumulated step
     * count must match to the last unit — the sharp check the old assignment bug
     * would fail (it reported only the final 8-primary batch).  (iterations is a
     * per-wavefront-pass counter, not per-history, so it legitimately differs
     * across batchings — each batch re-runs its own wavefront — and is not
     * compared here.) */
    ASSERT_TRUE(prof_final.steps > 0ull);
    ASSERT_TRUE(prof_live.steps == prof_final.steps);
    ASSERT_TRUE(prof_live.iterations > 0ull);

    /* Wall time and phase timers must stay a consistent decomposition after
     * accumulation: the summed phases cannot exceed the summed transport wall
     * time (small slack for clock granularity).  Under the old bug transport_s
     * carried only the last batch while phase_step_s summed all four, so this
     * would fail. */
    ASSERT_TRUE(prof_live.transport_s > 0.0);
    ASSERT_TRUE(prof_live.phase_step_s > 0.0);
    phase_sum_s = prof_live.phase_fill_s + prof_live.phase_zone_ref_s + prof_live.phase_distance_s
                  + prof_live.phase_step_s + prof_live.phase_compact_s;
    ASSERT_TRUE(phase_sum_s <= prof_live.transport_s * 1.02 + 1.0e-6);
}

/* ---- Entry point --------------------------------------------------------- */

int main(void) {
    test_create_rejects_null_args();
    test_create_rejects_unprepared_geometry();
    test_create_rejects_unknown_zone_material();
    test_create_free_lifecycle();
    test_profiling_run_lifecycle();
    test_create_rejects_voxel_geometry_without_hutable();
    test_clean_stop_partial_result_is_exact();
    test_no_run_control_runs_to_completion();
    test_checkpoint_final_only_default_completes();
    test_checkpoint_live_batches_match_final_only();
    test_checkpoint_live_batches_drain_families_with_nuclear();
    test_checkpoint_nuclear_invariant_across_batches_and_capacity();
    test_checkpoint_nuclear_overflow_is_counted_not_silent();
    test_checkpoint_profiling_accumulates_across_batches();
    test_dump_control_is_non_destructive();
    test_dump_write_failure_is_fail_soft();
    return 0;
}
