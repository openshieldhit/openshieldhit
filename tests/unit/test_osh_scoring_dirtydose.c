#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_step.h"
#include "material/runtime/osh_material_runtime.h"
#include "openshieldhit/const.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_estimator_common.h"
#include "scoring/runtime/osh_scoring_point.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/runtime/osh_scoring_step.h"

/*
 * DIRTYDOSE / DIRTYDOSEGY behave exactly like DOSE / DOSEGY but only accept the
 * contribution of charged particles whose mass stopping power (mass-LET) in the
 * scoring medium exceeds OSH_DIRTYDOSE_MASS_SP_THRESHOLD [MeV*cm^2/g].  These tests
 * exercise the LET gate on the step path, the point path, and under a Settings
 * medium override (dose-to-water, where the gate also uses mass-LET in water).
 */

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static int tmp_counter = 0;
static double const proton_mass_mev = OSH_PART_MASS_PROTON;

static void assert_close(double a, double b) {
    ASSERT_TRUE(fabs(a - b) < 1.0e-12);
}

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_scoring_dirtydose_test_%d.tmp", tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static struct osh_scoring_page_runtime *find_page_by_kind(struct osh_scoring_runtime *rt,
                                                          enum osh_scoring_score_kind kind) {
    size_t i;
    for (i = 0; i < rt->npages; ++i) {
        if (rt->pages[i].score_kind == kind) {
            return &rt->pages[i];
        }
    }
    return NULL;
}

/* Build a flat single-material SP table with the given mass stopping power (the
 * same value at both energy grid points, so the midpoint lookup returns it exactly).
 * Caller owns the backing arrays passed in and must keep them alive for the test. */
static void init_flat_mat_tables(struct osh_material_runtime *mat_rt,
                                 unsigned int *proj_z,
                                 unsigned int *proj_a,
                                 double *proj_mass,
                                 float *rho_arr,
                                 float *sp_values,
                                 float mass_sp) {
    proj_z[0] = 1u;
    proj_a[0] = 1u;
    proj_mass[0] = proton_mass_mev;
    rho_arr[0] = 1.0f;
    sp_values[0] = mass_sp;
    sp_values[1] = mass_sp;

    memset(mat_rt, 0, sizeof(*mat_rt));
    mat_rt->nprojectiles = 1u;
    mat_rt->nmaterials = 1u;
    mat_rt->nenergy = 2u;
    mat_rt->emin = 1.0;
    mat_rt->emax = 1000.0;
    mat_rt->log_emin = log(mat_rt->emin);
    mat_rt->inv_dlog = 1.0 / (log(mat_rt->emax) - log(mat_rt->emin));
    mat_rt->projectile_z = proj_z;
    mat_rt->projectile_a = proj_a;
    mat_rt->projectile_mass_mev = proj_mass;
    mat_rt->rho = rho_arr;
    mat_rt->mass_stopping_power = sp_values;
}

/* Score one Z-directed proton step through a 1x1x3 mesh (1 cm bins), with a flat
 * SP table of the given mass stopping power, and check DirtyDose against DOSE.
 * When above_threshold is expected, DirtyDose equals DOSE (0.5, 1.0, 0.5) and
 * DirtyDoseGy is that times OSH_MEVG2GY; otherwise DirtyDose is all zeros. */
static void run_step_gate_case(float mass_sp, int expect_scored) {
    char path[512];
    char const *detect = "Geometry Mesh\n"
                         "    Name G\n"
                         "    X -0.5 0.5 1\n"
                         "    Y -0.5 0.5 1\n"
                         "    Z  0.0 3.0 3\n"
                         "\n"
                         "Output\n"
                         "    Filename out.bdo\n"
                         "    Geo G\n"
                         "    Quantity Dose\n"
                         "    Quantity DirtyDose\n"
                         "    Quantity DirtyDoseGy\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dose_page;
    struct osh_scoring_page_runtime *dirty_page;
    struct osh_scoring_page_runtime *dirtygy_page;
    unsigned int proj_z[1];
    unsigned int proj_a[1];
    double proj_mass[1];
    float rho_arr[1];
    float sp_values[2];
    struct osh_material_runtime mat_rt;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), detect);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    init_flat_mat_tables(&mat_rt, proj_z, proj_a, proj_mass, rho_arr, sp_values, mass_sp);
    rt.mat_tables = &mat_rt;

    memset(&part, 0, sizeof(part));
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;
    part.mass = proton_mass_mev;

    memset(&st, 0, sizeof(st));
    st.p[2] = 0.5;
    st.p[3] = 150.0;
    st.q[2] = 2.5;
    st.q[3] = 150.0;
    st.v[2] = 1.0;
    st.ds = 2.0;
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    /* Postprocess so the DirtyDoseGy page gets its MeV/g -> Gy conversion; the
     * 1 cm^3 bins make the volume division a no-op, so DOSE stays 0.5, 1, 0.5. */
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    dose_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DOSE);
    dirty_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DIRTYDOSE);
    dirtygy_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DIRTYDOSEGY);
    ASSERT_TRUE(dose_page != NULL);
    ASSERT_TRUE(dirty_page != NULL);
    ASSERT_TRUE(dirtygy_page != NULL);

    /* DOSE is unaffected by the LET gate: base_scale = de/(score_len*rho) = 1.0. */
    assert_close(dose_page->acc.data[0], 0.5);
    assert_close(dose_page->acc.data[1], 1.0);
    assert_close(dose_page->acc.data[2], 0.5);

    if (expect_scored) {
        assert_close(dirty_page->acc.data[0], 0.5);
        assert_close(dirty_page->acc.data[1], 1.0);
        assert_close(dirty_page->acc.data[2], 0.5);
        assert_close(dirtygy_page->acc.data[0], 0.5 * OSH_MEVG2GY);
        assert_close(dirtygy_page->acc.data[1], 1.0 * OSH_MEVG2GY);
        assert_close(dirtygy_page->acc.data[2], 0.5 * OSH_MEVG2GY);
    } else {
        assert_close(dirty_page->acc.data[0], 0.0);
        assert_close(dirty_page->acc.data[1], 0.0);
        assert_close(dirty_page->acc.data[2], 0.0);
        assert_close(dirtygy_page->acc.data[0], 0.0);
        assert_close(dirtygy_page->acc.data[1], 0.0);
        assert_close(dirtygy_page->acc.data[2], 0.0);
    }

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* A proton whose mass-LET is above the threshold contributes exactly like DOSE. */
static void test_step_high_let_scores_like_dose(void) {
    run_step_gate_case((float) (OSH_DIRTYDOSE_MASS_SP_THRESHOLD + 10.0), 1);
}

/* A proton whose mass-LET is below the threshold contributes nothing. */
static void test_step_low_let_scores_zero(void) {
    run_step_gate_case((float) (OSH_DIRTYDOSE_MASS_SP_THRESHOLD - 10.0), 0);
}

/* A neutral particle has no mass stopping power, so it never scores dirty dose even
 * though it does deposit ordinary DOSE via its energy loss. */
static void test_step_neutral_scores_zero(void) {
    char path[512];
    char const *detect = "Geometry Mesh\n"
                         "    Name G\n"
                         "    X -0.5 0.5 1\n"
                         "    Y -0.5 0.5 1\n"
                         "    Z  0.0 3.0 3\n"
                         "\n"
                         "Output\n"
                         "    Filename out.bdo\n"
                         "    Geo G\n"
                         "    Quantity Dose\n"
                         "    Quantity DirtyDose\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dose_page;
    struct osh_scoring_page_runtime *dirty_page;
    unsigned int proj_z[1];
    unsigned int proj_a[1];
    double proj_mass[1];
    float rho_arr[1];
    float sp_values[2];
    struct osh_material_runtime mat_rt;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), detect);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    init_flat_mat_tables(
        &mat_rt, proj_z, proj_a, proj_mass, rho_arr, sp_values, (float) (OSH_DIRTYDOSE_MASS_SP_THRESHOLD + 10.0));
    rt.mat_tables = &mat_rt;

    /* Neutron: charge 0, z 0 — no SP-table column, so the gate value is 0. */
    memset(&part, 0, sizeof(part));
    part.charge = 0;
    part.z = 0u;
    part.a = 1u;

    memset(&st, 0, sizeof(st));
    st.p[2] = 0.5;
    st.p[3] = 150.0;
    st.q[2] = 2.5;
    st.q[3] = 150.0;
    st.v[2] = 1.0;
    st.ds = 2.0;
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    /* Assert on postprocessed output (1 cm^3 bins -> volume division is a no-op). */
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    dose_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DOSE);
    dirty_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DIRTYDOSE);
    ASSERT_TRUE(dose_page != NULL);
    ASSERT_TRUE(dirty_page != NULL);

    /* Ordinary DOSE still books the neutral's energy deposit. */
    assert_close(dose_page->acc.data[1], 1.0);
    /* DirtyDose rejects it entirely. */
    assert_close(dirty_page->acc.data[0], 0.0);
    assert_close(dirty_page->acc.data[1], 0.0);
    assert_close(dirty_page->acc.data[2], 0.0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* With a Settings medium override, both the scored dose and the LET gate refer to
 * the override medium.  Transport medium 0 has mass-LET below threshold (no dirty
 * dose), while override medium 1 (water) is above threshold: the overridden page
 * both passes the gate and scores dose-to-water (scaled by the S(water)/S(medium0)
 * stopping-power ratio). */
static void test_step_water_override_gates_and_scales(void) {
    char path[512];
    char const *detect = "Geometry Mesh\n"
                         "    Name G\n"
                         "    X -0.5 0.5 1\n"
                         "    Y -0.5 0.5 1\n"
                         "    Z  0.0 3.0 3\n"
                         "\n"
                         "Settings\n"
                         "    Name toWater\n"
                         "    Medium 1\n"
                         "\n"
                         "Output\n"
                         "    Filename out.bdo\n"
                         "    Geo G\n"
                         "    Quantity DirtyDose\n"
                         "    Quantity DirtyDose toWater\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dirty_local;
    struct osh_scoring_page_runtime *dirty_water;
    unsigned int proj_z[1];
    unsigned int proj_a[1];
    double proj_mass[1];
    float rho_arr[2];
    /* Two materials, one projectile, two energy points: layout is
     * [mat][proj][energy].  Medium 0 (transport) below threshold, medium 1 (water)
     * above, with S(water)/S(medium0) == 2. */
    float sp_values[4];
    struct osh_material_runtime mat_rt;
    float const sp_transport = (float) (OSH_DIRTYDOSE_MASS_SP_THRESHOLD - 10.0);
    float const sp_water = 2.0f * sp_transport;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), detect);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    proj_z[0] = 1u;
    proj_a[0] = 1u;
    proj_mass[0] = proton_mass_mev;
    rho_arr[0] = 1.0f;
    rho_arr[1] = 1.0f;
    sp_values[0] = sp_transport; /* mat 0 */
    sp_values[1] = sp_transport;
    sp_values[2] = sp_water; /* mat 1 */
    sp_values[3] = sp_water;

    memset(&mat_rt, 0, sizeof(mat_rt));
    mat_rt.nprojectiles = 1u;
    mat_rt.nmaterials = 2u;
    mat_rt.nenergy = 2u;
    mat_rt.emin = 1.0;
    mat_rt.emax = 1000.0;
    mat_rt.log_emin = log(mat_rt.emin);
    mat_rt.inv_dlog = 1.0 / (log(mat_rt.emax) - log(mat_rt.emin));
    mat_rt.projectile_z = proj_z;
    mat_rt.projectile_a = proj_a;
    mat_rt.projectile_mass_mev = proj_mass;
    mat_rt.rho = rho_arr;
    mat_rt.mass_stopping_power = sp_values;
    rt.mat_tables = &mat_rt;

    memset(&part, 0, sizeof(part));
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;
    part.mass = proton_mass_mev;

    memset(&st, 0, sizeof(st));
    st.p[2] = 0.5;
    st.p[3] = 150.0;
    st.q[2] = 2.5;
    st.q[3] = 150.0;
    st.v[2] = 1.0;
    st.ds = 2.0;
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    /* Assert on postprocessed output (1 cm^3 bins -> volume division is a no-op). */
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(rt.noutputs == 1u);
    ASSERT_TRUE(rt.outputs[0].npages == 2u);
    dirty_local = &rt.pages[rt.outputs[0].page_indices[0]];
    dirty_water = &rt.pages[rt.outputs[0].page_indices[1]];

    /* Local (transport) medium is below threshold: gated out entirely. */
    assert_close(dirty_local->acc.data[0], 0.0);
    assert_close(dirty_local->acc.data[1], 0.0);
    assert_close(dirty_local->acc.data[2], 0.0);

    /* Override to water: mass-LET in water is above threshold, so the page scores,
     * and the dose is dose-to-water = dose-to-medium * S(water)/S(medium0) = 2x. */
    assert_close(dirty_water->acc.data[0], 1.0);
    assert_close(dirty_water->acc.data[1], 2.0);
    assert_close(dirty_water->acc.data[2], 1.0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* Point-deposit path: a low-energy fragment born below transport threshold books
 * de/rho at one bin only when its mass-LET clears the gate. */
static void test_point_gate(void) {
    char path[512];
    char const *detect = "Geometry Mesh\n"
                         "    Name G\n"
                         "    X -0.5 0.5 1\n"
                         "    Y -0.5 0.5 1\n"
                         "    Z  0.0 3.0 3\n"
                         "\n"
                         "Output\n"
                         "    Filename out.bdo\n"
                         "    Geo G\n"
                         "    Quantity DirtyDose\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dirty_page;
    unsigned int proj_z[1];
    unsigned int proj_a[1];
    double proj_mass[1];
    float rho_arr[1];
    float sp_values[2];
    struct osh_material_runtime mat_rt;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), detect);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* Above threshold first: the point deposits de/rho = 4/2 = 2 at bin 1. */
    init_flat_mat_tables(
        &mat_rt, proj_z, proj_a, proj_mass, rho_arr, sp_values, (float) (OSH_DIRTYDOSE_MASS_SP_THRESHOLD + 10.0));
    rt.mat_tables = &mat_rt;

    memset(&part, 0, sizeof(part));
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;
    part.mass = proton_mass_mev;

    memset(&st, 0, sizeof(st));
    st.p[2] = 1.5; /* inside mesh bin 1 (Z in [1,2)) */
    st.p[3] = 150.0;
    st.q[2] = 1.5;
    st.q[3] = 150.0; /* point: q == p */
    st.de = 4.0;
    st.rho = 2.0;
    st.wt = 1.0;
    st.medium = 0;

    rc = osh_scoring_score_point(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    /* Assert on postprocessed output (1 cm^3 bins -> volume division is a no-op). */
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    dirty_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DIRTYDOSE);
    ASSERT_TRUE(dirty_page != NULL);
    assert_close(dirty_page->acc.data[0], 0.0);
    assert_close(dirty_page->acc.data[1], 2.0);
    assert_close(dirty_page->acc.data[2], 0.0);

    /* Now below threshold: the same point books nothing (raw check — the deposit is
     * zero before any postprocessing). */
    memset(dirty_page->acc.data, 0, dirty_page->len * sizeof(*dirty_page->acc.data));
    init_flat_mat_tables(
        &mat_rt, proj_z, proj_a, proj_mass, rho_arr, sp_values, (float) (OSH_DIRTYDOSE_MASS_SP_THRESHOLD - 10.0));

    rc = osh_scoring_score_point(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);
    assert_close(dirty_page->acc.data[1], 0.0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

int main(void) {
    test_step_high_let_scores_like_dose();
    test_step_low_let_scores_zero();
    test_step_neutral_scores_zero();
    test_step_water_override_gates_and_scales();
    test_point_gate();
    return 0;
}
