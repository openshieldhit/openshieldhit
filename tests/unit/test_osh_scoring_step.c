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
#include "particle/osh_particle_pdg.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/runtime/osh_scoring_step.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void assert_close(double a, double b) {
    ASSERT_TRUE(fabs(a - b) < 1.0e-12);
}

static void test_score_mesh_energy_and_fluence_with_filters(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    enum osh_status rc;
    size_t i;
    size_t energy0_idx;
    size_t energy1_idx;
    size_t fluence_idx;
    size_t filtered_idx;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test01/detect.dat", OSH_PROJECT_SOURCE_DIR);

    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&part, 0, sizeof(part));
    part.mass = 12.0 * 931.49410242;
    part.pdg = 1000060120;
    part.charge = 6;
    part.z = 6u;
    part.a = 12u;
    part.is_nucleus = 1u;

    memset(&st, 0, sizeof(st));
    st.p[0] = 0.0;
    st.p[1] = 0.0;
    st.p[2] = 0.5;
    st.p[3] = 100.0;
    st.q[0] = 0.0;
    st.q[1] = 0.0;
    st.q[2] = 8.5;
    st.q[3] = 98.0;
    st.v[0] = 0.0;
    st.v[1] = 0.0;
    st.v[2] = 1.0;
    st.w[0] = 0.0;
    st.w[1] = 0.0;
    st.w[2] = 1.0;
    st.ds = 8.0;
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 1;
    st.zone = 1;
    st.prim_idx = 7u;
    st.gen = 0u;

    rc = osh_scoring_score_step(&rt, &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    energy0_idx = rt.outputs[0].page_indices[0];
    energy1_idx = rt.outputs[1].page_indices[0];
    fluence_idx = rt.outputs[1].page_indices[1];
    filtered_idx = rt.outputs[1].page_indices[2];

    assert_close(rt.pages[energy0_idx].acc.data[0], 0.875);
    assert_close(rt.pages[energy0_idx].acc.data[1], 1.0);
    assert_close(rt.pages[energy0_idx].acc.data[2], 0.125);

    assert_close(rt.pages[energy1_idx].acc.data[0], 0.875);
    assert_close(rt.pages[energy1_idx].acc.data[1], 1.0);
    assert_close(rt.pages[energy1_idx].acc.data[2], 0.125);

    assert_close(rt.pages[fluence_idx].acc.data[0], 0.875);
    assert_close(rt.pages[fluence_idx].acc.data[1], 1.0);
    assert_close(rt.pages[fluence_idx].acc.data[2], 0.125);

    assert_close(rt.pages[filtered_idx].acc.data[0], 0.875);
    assert_close(rt.pages[filtered_idx].acc.data[1], 1.0);
    assert_close(rt.pages[filtered_idx].acc.data[2], 0.125);

    st.gen = 1u;
    rc = osh_scoring_score_step(&rt, &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    assert_close(rt.pages[energy0_idx].acc.data[0], 1.75);
    assert_close(rt.pages[energy0_idx].acc.data[1], 2.0);
    assert_close(rt.pages[energy0_idx].acc.data[2], 0.25);

    assert_close(rt.pages[energy1_idx].acc.data[0], 1.75);
    assert_close(rt.pages[energy1_idx].acc.data[1], 2.0);
    assert_close(rt.pages[energy1_idx].acc.data[2], 0.25);

    assert_close(rt.pages[fluence_idx].acc.data[0], 1.75);
    assert_close(rt.pages[fluence_idx].acc.data[1], 2.0);
    assert_close(rt.pages[fluence_idx].acc.data[2], 0.25);

    assert_close(rt.pages[filtered_idx].acc.data[0], 0.875);
    assert_close(rt.pages[filtered_idx].acc.data[1], 1.0);
    assert_close(rt.pages[filtered_idx].acc.data[2], 0.125);

    for (i = 3u; i < rt.pages[energy0_idx].len; ++i) {
        assert_close(rt.pages[energy0_idx].acc.data[i], 0.0);
        assert_close(rt.pages[energy1_idx].acc.data[i], 0.0);
        assert_close(rt.pages[fluence_idx].acc.data[i], 0.0);
        assert_close(rt.pages[filtered_idx].acc.data[i], 0.0);
    }

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

static void test_score_mesh_uses_step_chord_after_bending(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    enum osh_status rc;
    size_t i;
    size_t energy0_idx;
    size_t energy1_idx;
    size_t fluence_idx;
    size_t filtered_idx;
    double chord_len;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test01/detect.dat", OSH_PROJECT_SOURCE_DIR);

    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&part, 0, sizeof(part));
    part.mass = 12.0 * 931.49410242;
    part.pdg = 1000060120;
    part.charge = 6;
    part.z = 6u;
    part.a = 12u;
    part.is_nucleus = 1u;

    memset(&st, 0, sizeof(st));
    st.p[0] = 1.0;
    st.p[1] = 0.0;
    st.p[2] = 0.5;
    st.p[3] = 100.0;
    st.q[0] = 0.0;
    st.q[1] = 0.0;
    st.q[2] = 8.5;
    st.q[3] = 98.0;
    st.v[0] = 0.0;
    st.v[1] = 0.0;
    st.v[2] = 1.0;
    st.w[0] = -0.1;
    st.w[1] = 0.0;
    st.w[2] = 0.99498743710662;
    st.ds = sqrt(65.0);
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 1;
    st.zone = 1;
    st.prim_idx = 7u;
    st.gen = 0u;

    rc = osh_scoring_score_step(&rt, &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    energy0_idx = rt.outputs[0].page_indices[0];
    energy1_idx = rt.outputs[1].page_indices[0];
    fluence_idx = rt.outputs[1].page_indices[1];
    filtered_idx = rt.outputs[1].page_indices[2];
    chord_len = sqrt(65.0);

    assert_close(rt.pages[energy0_idx].acc.data[0], 0.0);
    assert_close(rt.pages[energy0_idx].acc.data[1], 0.875);
    assert_close(rt.pages[energy0_idx].acc.data[2], 0.125);

    assert_close(rt.pages[energy1_idx].acc.data[0], 0.0);
    assert_close(rt.pages[energy1_idx].acc.data[1], 0.875);
    assert_close(rt.pages[energy1_idx].acc.data[2], 0.125);

    assert_close(rt.pages[fluence_idx].acc.data[0], 0.0);
    assert_close(rt.pages[fluence_idx].acc.data[1], 3.5 * chord_len / 32.0);
    assert_close(rt.pages[fluence_idx].acc.data[2], 0.5 * chord_len / 32.0);

    assert_close(rt.pages[filtered_idx].acc.data[0], 0.0);
    assert_close(rt.pages[filtered_idx].acc.data[1], 3.5 * chord_len / 32.0);
    assert_close(rt.pages[filtered_idx].acc.data[2], 0.5 * chord_len / 32.0);

    for (i = 3u; i < rt.pages[energy0_idx].len; ++i) {
        assert_close(rt.pages[energy0_idx].acc.data[i], 0.0);
        assert_close(rt.pages[energy1_idx].acc.data[i], 0.0);
        assert_close(rt.pages[fluence_idx].acc.data[i], 0.0);
        assert_close(rt.pages[filtered_idx].acc.data[i], 0.0);
    }

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

static void test_score_mesh_neutron_id_filter(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle neutron;
    struct step st;
    enum osh_status rc;
    size_t fluence_idx;
    size_t neutron_idx;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test_neutron_filter/detect.dat", OSH_PROJECT_SOURCE_DIR);

    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&neutron, 0, sizeof(neutron));
    neutron.mass = OSH_PART_MASS_NEUTRON;
    neutron.pdg = OSH_PART_PDG_NEUTRON;
    neutron.charge = 0;
    neutron.z = 0u;
    neutron.a = 1u;
    neutron.is_nucleus = 0u;

    memset(&st, 0, sizeof(st));
    st.p[0] = 0.0;
    st.p[1] = 0.0;
    st.p[2] = 0.5;
    st.p[3] = 10.0;
    st.q[0] = 0.0;
    st.q[1] = 0.0;
    st.q[2] = 2.5;
    st.q[3] = 10.0;
    st.v[0] = 0.0;
    st.v[1] = 0.0;
    st.v[2] = 1.0;
    st.w[0] = 0.0;
    st.w[1] = 0.0;
    st.w[2] = 1.0;
    st.ds = 2.0;
    st.de = 0.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;
    st.zone = 0;
    st.prim_idx = 1u;
    st.gen = 1u;

    rc = osh_scoring_score_step(&rt, &neutron, &st);
    ASSERT_TRUE(rc == OSH_OK);

    fluence_idx = rt.outputs[0].page_indices[0];
    neutron_idx = rt.outputs[0].page_indices[1];

    assert_close(rt.pages[fluence_idx].acc.data[0], 0.5);
    assert_close(rt.pages[fluence_idx].acc.data[1], 1.0);
    assert_close(rt.pages[fluence_idx].acc.data[2], 0.5);
    assert_close(rt.pages[neutron_idx].acc.data[0], 0.5);
    assert_close(rt.pages[neutron_idx].acc.data[1], 1.0);
    assert_close(rt.pages[neutron_idx].acc.data[2], 0.5);

    neutron.pdg = OSH_PART_PDG_PROTON;
    rc = osh_scoring_score_step(&rt, &neutron, &st);
    ASSERT_TRUE(rc == OSH_OK);

    assert_close(rt.pages[fluence_idx].acc.data[0], 1.0);
    assert_close(rt.pages[fluence_idx].acc.data[1], 2.0);
    assert_close(rt.pages[fluence_idx].acc.data[2], 1.0);
    assert_close(rt.pages[neutron_idx].acc.data[0], 0.5);
    assert_close(rt.pages[neutron_idx].acc.data[1], 1.0);
    assert_close(rt.pages[neutron_idx].acc.data[2], 0.5);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
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

/* Step along Z through a 1x1x3 mesh (1cm bins).
 * p=(0,0,0.5) q=(0,0,2.5): crosses bins 0,1,2 with path_len 0.5, 1.0, 0.5 cm.
 * Verifies DOSE, DLET, TLET accumulation and postprocess without mat_tables. */
static void test_score_mesh_dose_and_let_geometric(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dose_page;
    struct osh_scoring_page_runtime *dlet_page;
    struct osh_scoring_page_runtime *tlet_page;
    struct osh_scoring_page_runtime *dqeff_page;
    struct osh_scoring_page_runtime *tqeff_page;
    enum osh_status rc;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test_let/detect.dat", OSH_PROJECT_SOURCE_DIR);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    /* mat_tables stays NULL → geometric LET fallback, DQEFF/TQEFF skipped */

    memset(&part, 0, sizeof(part));
    part.z = 1u;
    part.a = 1u;

    memset(&st, 0, sizeof(st));
    st.p[0] = 0.0;
    st.p[1] = 0.0;
    st.p[2] = 0.5;
    st.p[3] = 150.0;
    st.q[0] = 0.0;
    st.q[1] = 0.0;
    st.q[2] = 2.5;
    st.q[3] = 148.0;
    st.v[0] = 0.0;
    st.v[1] = 0.0;
    st.v[2] = 1.0;
    st.ds = 2.0;
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;

    rc = osh_scoring_score_step(&rt, &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* DOSE: base_scale = de*vol_inv/(score_len*rho) = 2.0*1.0/(2.0*1.0) = 1.0 MeV/g
     * No Gy conversion applied — DOSE is in MeV/g (SH12A-compatible). */
    dose_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DOSE);
    ASSERT_TRUE(dose_page != NULL);
    assert_close(dose_page->acc.data[0], 0.5);
    assert_close(dose_page->acc.data[1], 1.0);
    assert_close(dose_page->acc.data[2], 0.5);

    /* DOSEGY: same accumulated values, converted to Gy via OSH_MEVG2GY. */
    {
        struct osh_scoring_page_runtime *dosegy_page;
        dosegy_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DOSEGY);
        ASSERT_TRUE(dosegy_page != NULL);
        assert_close(dosegy_page->acc.data[0], 0.5 * OSH_MEVG2GY);
        assert_close(dosegy_page->acc.data[1], 1.0 * OSH_MEVG2GY);
        assert_close(dosegy_page->acc.data[2], 0.5 * OSH_MEVG2GY);
    }

    /* DLET geometric: let_step = de/score_len = 1.0 MeV/cm in all traversed bins */
    dlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DLET);
    ASSERT_TRUE(dlet_page != NULL);
    assert_close(dlet_page->acc.data[0], 1.0);
    assert_close(dlet_page->acc.data[1], 1.0);
    assert_close(dlet_page->acc.data[2], 1.0);

    /* TLET geometric: same value as DLET for geometric fallback */
    tlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_TLET);
    ASSERT_TRUE(tlet_page != NULL);
    assert_close(tlet_page->acc.data[0], 1.0);
    assert_close(tlet_page->acc.data[1], 1.0);
    assert_close(tlet_page->acc.data[2], 1.0);

    /* DQEFF/TQEFF: mat_tables NULL → no contribution, data stays zero */
    dqeff_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DQEFF);
    ASSERT_TRUE(dqeff_page != NULL);
    assert_close(dqeff_page->acc.data[0], 0.0);

    tqeff_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_TQEFF);
    ASSERT_TRUE(tqeff_page != NULL);
    assert_close(tqeff_page->acc.data[0], 0.0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* Verify DQEFF and TQEFF produce (z_eff/β)² averaged over the step.
 * For a monoenergetic proton DQEFF == TQEFF == (z_eff/β)². */
static void test_score_mesh_dqeff_tqeff(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dqeff_page;
    struct osh_scoring_page_runtime *tqeff_page;
    /* Minimal mat_tables with proton rest mass only. */
    double proton_mass_mev = 938.272046;
    unsigned int proj_z[1];
    unsigned int proj_a[1];
    double proj_mass[1];
    float rho_arr[1];
    float sp_dummy[1];
    struct osh_material_runtime mat_rt;
    double mean_energy;
    double gamma_inv;
    double beta;
    double expected_qeff;
    enum osh_status rc;

    proj_z[0] = 1u;
    proj_a[0] = 1u;
    proj_mass[0] = proton_mass_mev;
    rho_arr[0] = 1.0f;
    sp_dummy[0] = 0.0f;

    memset(&mat_rt, 0, sizeof(mat_rt));
    mat_rt.nprojectiles = 1u;
    mat_rt.nmaterials = 1u;
    mat_rt.projectile_z = proj_z;
    mat_rt.projectile_a = proj_a;
    mat_rt.projectile_mass_mev = proj_mass;
    mat_rt.rho = rho_arr;
    mat_rt.mass_stopping_power = sp_dummy;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test_let/detect.dat", OSH_PROJECT_SOURCE_DIR);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    rt.mat_tables = &mat_rt;

    memset(&part, 0, sizeof(part));
    part.z = 1u;
    part.a = 1u;

    memset(&st, 0, sizeof(st));
    st.p[0] = 0.0;
    st.p[1] = 0.0;
    st.p[2] = 0.5;
    st.p[3] = 150.0;
    st.q[0] = 0.0;
    st.q[1] = 0.0;
    st.q[2] = 2.5;
    st.q[3] = 150.0; /* monoenergetic */
    st.v[0] = 0.0;
    st.v[1] = 0.0;
    st.v[2] = 1.0;
    st.ds = 2.0;
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;

    rc = osh_scoring_score_step(&rt, &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* Expected qeff = (z_eff/β)² for proton at 150 MeV */
    mean_energy = 150.0;
    gamma_inv = proton_mass_mev / (mean_energy + proton_mass_mev);
    beta = sqrt(1.0 - gamma_inv * gamma_inv);
    /* z_eff ≈ 1.0 for proton at β≈0.507 (125·β·1^(-2/3) ≈ 63 → exp term ≈ 0) */
    expected_qeff = 1.0 / (beta * beta);

    dqeff_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DQEFF);
    ASSERT_TRUE(dqeff_page != NULL);
    /* All traversed bins share the same qeff (monoenergetic) */
    assert_close(dqeff_page->acc.data[0], expected_qeff);
    assert_close(dqeff_page->acc.data[1], expected_qeff);
    assert_close(dqeff_page->acc.data[2], expected_qeff);

    tqeff_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_TQEFF);
    ASSERT_TRUE(tqeff_page != NULL);
    assert_close(tqeff_page->acc.data[0], expected_qeff);
    assert_close(tqeff_page->acc.data[1], expected_qeff);
    assert_close(tqeff_page->acc.data[2], expected_qeff);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

int main(void) {
    test_score_mesh_energy_and_fluence_with_filters();
    test_score_mesh_uses_step_chord_after_bending();
    test_score_mesh_neutron_id_filter();
    test_score_mesh_dose_and_let_geometric();
    test_score_mesh_dqeff_tqeff();
    return 0;
}
