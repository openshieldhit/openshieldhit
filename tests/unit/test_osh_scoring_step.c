#include <math.h>
#include <stdint.h>
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

static int tmp_counter = 0;

static void assert_close(double a, double b) {
    ASSERT_TRUE(fabs(a - b) < 1.0e-12);
}

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_scoring_step_test_%d.tmp", tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
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

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
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

    /* FLUENCE now deposits the raw track length (cm) per bin; the ÷volume happens
     * in postprocess.  Voxel volume is 4 cm³ (dz=4), so the pre-postprocess raw is
     * the old fluence value / vol_inv (× 4). */
    assert_close(rt.pages[fluence_idx].acc.data[0], 3.5);
    assert_close(rt.pages[fluence_idx].acc.data[1], 4.0);
    assert_close(rt.pages[fluence_idx].acc.data[2], 0.5);

    assert_close(rt.pages[filtered_idx].acc.data[0], 3.5);
    assert_close(rt.pages[filtered_idx].acc.data[1], 4.0);
    assert_close(rt.pages[filtered_idx].acc.data[2], 0.5);

    st.gen = 1u;
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
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

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
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

    /* FLUENCE deposits the raw track length (÷volume is deferred to postprocess);
     * the old values had vol_inv = 1/4 folded in (the /32 below was chord/8 · 1/4). */
    assert_close(rt.pages[fluence_idx].acc.data[0], 0.0);
    assert_close(rt.pages[fluence_idx].acc.data[1], 3.5 * chord_len / 8.0);
    assert_close(rt.pages[fluence_idx].acc.data[2], 0.5 * chord_len / 8.0);

    assert_close(rt.pages[filtered_idx].acc.data[0], 0.0);
    assert_close(rt.pages[filtered_idx].acc.data[1], 3.5 * chord_len / 8.0);
    assert_close(rt.pages[filtered_idx].acc.data[2], 0.5 * chord_len / 8.0);

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

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &neutron, &st);
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
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &neutron, &st);
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
    /* mat_tables stays NULL -> geometric LET fallback.  part.mass stays unset
     * here, so DQEFF/TQEFF are skipped because beta cannot be computed. */

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

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
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
    /* Two-point log energy grid: enough for osh_material_runtime_sp_lookup(),
     * which interpolates between grid points idx and idx+1 (so nenergy >= 2 and
     * the table must hold two entries per [material][projectile]). */
    float sp_dummy[2];
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
    /* Flat, non-zero stopping power. Its magnitude does not affect the Qeff
     * pages exercised here (they carry no dose-to-medium override, so the
     * looked-up value only feeds sp_tr, which stays unused); it only has to be a
     * valid in-bounds read. */
    sp_dummy[0] = 1.0f;
    sp_dummy[1] = 1.0f;

    memset(&mat_rt, 0, sizeof(mat_rt));
    mat_rt.nprojectiles = 1u;
    mat_rt.nmaterials = 1u;
    mat_rt.nenergy = 2u;
    mat_rt.emin = OSH_MATERIAL_RUNTIME_EMIN;
    mat_rt.emax = OSH_MATERIAL_RUNTIME_EMAX;
    mat_rt.log_emin = log(OSH_MATERIAL_RUNTIME_EMIN);
    mat_rt.inv_dlog = (double) (mat_rt.nenergy - 1u) / log(OSH_MATERIAL_RUNTIME_EMAX / OSH_MATERIAL_RUNTIME_EMIN);
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
    part.mass = proton_mass_mev;

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

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* Expected qeff = (z_eff/β)² for proton at 150 MeV */
    mean_energy = 150.0;
    gamma_inv = proton_mass_mev / (mean_energy + proton_mass_mev);
    beta = sqrt(1.0 - (gamma_inv * gamma_inv));
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

/* LET and QEFF can be differential-axis values for ordinary extensive scores
 * such as FLUENCE.  They are not the same thing as the averaged quantities DLET,
 * TLET, DQEFF, and TQEFF.  This test verifies that a FLUENCE page can bin track
 * length by LET and by QEFF. */
static void test_score_mesh_fluence_diff_let_qeff(void) {
    char path[512];
    char const *detect = "Geometry Mesh\n"
                         "    Name G\n"
                         "    X -0.5 0.5 1\n"
                         "    Y -0.5 0.5 1\n"
                         "    Z  0.0 2.0 1\n"
                         "\n"
                         "Output\n"
                         "    Filename out.bdo\n"
                         "    Geo G\n"
                         "    Quantity Fluence\n"
                         "    Diff1 0.0 10.0 5\n"
                         "    Diff1Type LET\n"
                         "    Quantity Fluence\n"
                         "    Diff1 0.0 8.0 8\n"
                         "    Diff1Type QEFF\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *let_page;
    struct osh_scoring_page_runtime *qeff_page;
    double proton_mass_mev;
    unsigned int proj_z[1];
    unsigned int proj_a[1];
    double proj_mass[1];
    float rho_arr[1];
    float sp_values[2];
    struct osh_material_runtime mat_rt;
    double mean_energy;
    double gamma_inv;
    double beta;
    double qeff;
    size_t qeff_bin;
    enum osh_status rc;
    size_t i;

    write_temp_file(path, sizeof(path), detect);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    proton_mass_mev = 938.272046;
    proj_z[0] = 1u;
    proj_a[0] = 1u;
    proj_mass[0] = proton_mass_mev;
    rho_arr[0] = 1.0f;
    sp_values[0] = 4.0f;
    sp_values[1] = 4.0f;

    memset(&mat_rt, 0, sizeof(mat_rt));
    mat_rt.nprojectiles = 1u;
    mat_rt.nmaterials = 1u;
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
    part.z = 1u;
    part.a = 1u;
    part.mass = proton_mass_mev;

    memset(&st, 0, sizeof(st));
    st.p[0] = 0.0;
    st.p[1] = 0.0;
    st.p[2] = 0.5;
    st.p[3] = 150.0;
    st.q[0] = 0.0;
    st.q[1] = 0.0;
    st.q[2] = 1.5;
    st.q[3] = 150.0;
    st.v[2] = 1.0;
    st.ds = 1.0;
    st.de = 1.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(rt.noutputs == 1u);
    ASSERT_TRUE(rt.outputs[0].npages == 2u);
    let_page = &rt.pages[rt.outputs[0].page_indices[0]];
    qeff_page = &rt.pages[rt.outputs[0].page_indices[1]];

    /* LET axis: SP is constant 4 MeV*cm2/g and rho is 1 g/cm3, so LET is
     * 4 MeV/cm.  Diff1 0..10 with 5 bins has width 2, so LET=4 lands in bin 2.
     * The raw FLUENCE deposit is the crossed track length, here 1 cm. */
    ASSERT_TRUE(let_page->len == 5u);
    for (i = 0u; i < let_page->len; ++i) {
        assert_close(let_page->acc.data[i], (i == 2u) ? 1.0 : 0.0);
    }

    /* QEFF axis: qeff is computed at the step midpoint and used only as the
     * differential-axis coordinate.  The scored quantity is still raw FLUENCE,
     * so the selected spectrum bin receives 1 cm of track length. */
    mean_energy = 150.0;
    gamma_inv = proton_mass_mev / (mean_energy + proton_mass_mev);
    beta = sqrt(1.0 - (gamma_inv * gamma_inv));
    qeff = 1.0 / (beta * beta);
    qeff_bin = (size_t) ((qeff - 0.0) / (8.0 - 0.0) * 8.0);
    ASSERT_TRUE(qeff_bin < qeff_page->len);
    for (i = 0u; i < qeff_page->len; ++i) {
        assert_close(qeff_page->acc.data[i], (i == qeff_bin) ? 1.0 : 0.0);
    }

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

static void test_score_zone_energy_fluence_dose(void) {
    char path[512];
    char const *detect = "Geometry Zone\n"
                         "    Name Z\n"
                         "    Zone Entrance\n"
                         "    Volume 2.0\n"
                         "    Zone Target\n"
                         "    Volume 4.0\n"
                         "\n"
                         "Output\n"
                         "    Filename zone.bdo\n"
                         "    Geo Z\n"
                         "    Quantity Energy\n"
                         "    Quantity Fluence\n"
                         "    Quantity Dose\n"
                         "    Quantity DoseGy\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *energy_page;
    struct osh_scoring_page_runtime *fluence_page;
    struct osh_scoring_page_runtime *dose_page;
    struct osh_scoring_page_runtime *dosegy_page;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), detect);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    ASSERT_TRUE(ws->ngeometries == 1u);
    ASSERT_TRUE(ws->geometries[0].nzone_indices == 2u);

    /* App-level parsing stores Zone selectors by name.  The OSH app resolves
     * those names against geo.dat before compile; this unit test fills the same
     * resolved 0-based transport zone ids directly to exercise the library path. */
    ws->geometries[0].zone_indices = (size_t *) calloc(2u, sizeof(*ws->geometries[0].zone_indices));
    ASSERT_TRUE(ws->geometries[0].zone_indices != NULL);
    ws->geometries[0].zone_indices[0] = 3u;
    ws->geometries[0].zone_indices[1] = 7u;

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&part, 0, sizeof(part));
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;

    memset(&st, 0, sizeof(st));
    st.p[2] = 0.0;
    st.q[2] = 2.0;
    st.v[2] = 1.0;
    st.ds = 2.0;
    st.de = 8.0;
    st.rho = 2.0;
    st.wt = 1.0;
    st.medium = 0;
    st.zone = 7; /* selected Zone entry #2, dense scorer bin 1 */

    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    energy_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_ENERGY);
    fluence_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_FLUENCE);
    dose_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DOSE);
    dosegy_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DOSEGY);
    ASSERT_TRUE(energy_page != NULL);
    ASSERT_TRUE(fluence_page != NULL);
    ASSERT_TRUE(dose_page != NULL);
    ASSERT_TRUE(dosegy_page != NULL);

    /* Before postprocess, Zone scoring has deposited into the dense zone bin
     * selected by st.zone.  The first selected zone is untouched. */
    assert_close(energy_page->acc.data[0], 0.0);
    assert_close(energy_page->acc.data[1], 8.0);
    assert_close(fluence_page->acc.data[0], 0.0);
    assert_close(fluence_page->acc.data[1], 2.0);
    assert_close(dose_page->acc.data[0], 0.0);
    assert_close(dose_page->acc.data[1], 4.0);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* Target has Volume 4 cm3.  FLUENCE becomes track length / volume = 2/4,
     * and DOSE becomes energy / (rho*volume) = 8/(2*4). */
    assert_close(energy_page->acc.data[0], 0.0);
    assert_close(energy_page->acc.data[1], 8.0);
    assert_close(fluence_page->acc.data[0], 0.0);
    assert_close(fluence_page->acc.data[1], 0.5);
    assert_close(dose_page->acc.data[0], 0.0);
    assert_close(dose_page->acc.data[1], 1.0);
    assert_close(dosegy_page->acc.data[0], 0.0);
    assert_close(dosegy_page->acc.data[1], OSH_MEVG2GY);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* Allocate a private accumulator set cloning the shape of every page in rt
 * (same len and data2 presence), zero-initialised.  Mirrors what a future
 * worker_accumulators_alloc() does, but inline for the test. */
static struct osh_scoring_accumulator *alloc_private_acc_set(struct osh_scoring_runtime const *rt) {
    struct osh_scoring_accumulator *set;
    size_t i;

    set = (struct osh_scoring_accumulator *) calloc(rt->npages, sizeof(*set));
    /* calloc(0, ...) may legitimately return NULL; an empty set is valid. */
    ASSERT_TRUE(set != NULL || rt->npages == 0u);
    for (i = 0; i < rt->npages; ++i) {
        ASSERT_TRUE(osh_scoring_accumulator_alloc(&set[i], rt->pages[i].len, rt->pages[i].has_data2) == OSH_OK);
    }
    return set;
}

static void free_acc_set(struct osh_scoring_accumulator *set, size_t n) {
    size_t i;
    for (i = 0; i < n; ++i) {
        osh_scoring_accumulator_free(&set[i]);
    }
    free(set);
}

static void assert_arrays_bit_equal(double const *a, double const *b, size_t n) {
    size_t i;
    if (!a || !b) {
        ASSERT_TRUE(a == b);
        return;
    }
    /* Compare the IEEE-754 representations as integers: this is a true bit-for-bit
     * check (distinguishes +0.0 from -0.0, matches identical NaN payloads) rather
     * than the value comparison == would do.  memcpy into uint64_t avoids both the
     * value comparison and the bugprone-suspicious-memory-comparison clang-tidy
     * diagnostic that memcmp on a double (no unique object representation) trips. */
    for (i = 0; i < n; ++i) {
        uint64_t ua;
        uint64_t ub;
        memcpy(&ua, &a[i], sizeof(ua));
        memcpy(&ub, &b[i], sizeof(ub));
        ASSERT_TRUE(ua == ub);
    }
}

/* Acceptance criterion: scoring a step into a private acc_set and then merging
 * that set into a zeroed master must equal scoring the same step directly into
 * the master.  Both paths run identical floating-point ops and the merge is a
 * plain add onto zero, so the two results are bit-for-bit identical. */
static void test_score_private_then_merge_equals_direct(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_accumulator *priv;
    struct osh_scoring_accumulator *merged;
    enum osh_status rc;
    size_t i;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test_let/detect.dat", OSH_PROJECT_SOURCE_DIR);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&part, 0, sizeof(part));
    part.z = 1u;
    part.a = 1u;

    memset(&st, 0, sizeof(st));
    st.p[2] = 0.5;
    st.p[3] = 150.0;
    st.q[2] = 2.5;
    st.q[3] = 148.0;
    st.v[2] = 1.0;
    st.ds = 2.0;
    st.de = 2.0;
    st.rho = 1.0;
    st.wt = 1.0;
    st.medium = 0;

    /* Direct: deposit straight into the master view (rt->pages[*].acc). */
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    /* Private: deposit into a private set, then fold into a zeroed master. */
    priv = alloc_private_acc_set(&rt);
    merged = alloc_private_acc_set(&rt); /* starts zeroed */
    rc = osh_scoring_score_step(&rt, priv, osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);
    for (i = 0; i < rt.npages; ++i) {
        ASSERT_TRUE(osh_scoring_accumulator_merge(&merged[i], &priv[i]) == OSH_OK);
    }

    /* merged (zero + private) must equal the direct master deposit, bit-for-bit. */
    for (i = 0; i < rt.npages; ++i) {
        assert_arrays_bit_equal(merged[i].data, rt.pages[i].acc.data, rt.pages[i].len);
        assert_arrays_bit_equal(merged[i].data2, rt.pages[i].acc.data2, rt.pages[i].len);
    }

    free_acc_set(priv, rt.npages);
    free_acc_set(merged, rt.npages);
    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

int main(void) {
    test_score_mesh_energy_and_fluence_with_filters();
    test_score_mesh_uses_step_chord_after_bending();
    test_score_mesh_neutron_id_filter();
    test_score_mesh_dose_and_let_geometric();
    test_score_mesh_dqeff_tqeff();
    test_score_mesh_fluence_diff_let_qeff();
    test_score_zone_energy_fluence_dose();
    test_score_private_then_merge_equals_direct();
    return 0;
}
