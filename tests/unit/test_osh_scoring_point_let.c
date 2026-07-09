/*
 * Unit tests for the point-deposit LET scorers (osh_scoring_estimator_point_dlet /
 * _tlet, issue #227).  A heavy recoil or fragment born below the transport
 * threshold deposits its whole birth energy at one bin with no track length.  The
 * point LET scorers give it a representative dE/dx — the stopping power at its
 * birth energy S(medium, E_birth) * rho — so that the single highest-LET component
 * of the field becomes visible in the DLET/TLET scorers instead of contributing
 * energy/dose only.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_step.h"
#include "material/runtime/osh_material_runtime.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_point.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static int tmp_counter = 0;

/* LET ratios come from a few divides; compare with a relative tolerance that is
 * immune to floating-point noise but far tighter than the physics differences the
 * tests assert (DLET vs TLET differ by tens of MeV/cm). */
static void assert_near(double a, double b) {
    ASSERT_TRUE(fabs(a - b) < 1.0e-9 * (1.0 + fabs(b)));
}

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_scoring_point_let_test_%d.tmp", tmp_counter++);
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

/* detect.dat with a 1x1x3 mesh (1 cm bins) scoring both DLET and TLET.  The point
 * deposits land in Z bin 1; the volume cancels in the LET ratio, so the mesh sizing
 * is irrelevant beyond locating the bin. */
static char const *const DETECT_TEXT = "Geometry Mesh\n"
                                       "    Name G\n"
                                       "    X -0.5 0.5 1\n"
                                       "    Y -0.5 0.5 1\n"
                                       "    Z  0.0 3.0 3\n"
                                       "\n"
                                       "Output\n"
                                       "    Filename out.bdo\n"
                                       "    Geo G\n"
                                       "    Quantity DLET\n"
                                       "    Quantity TLET\n";

/* Build a flat SP table: two projectiles (Z=2 alpha, Z=6 carbon), one material,
 * two energy grid points holding the same mass stopping power per projectile so the
 * midpoint lookup returns it exactly regardless of the birth energy.  rho is 1
 * g/cm^3, so LET [MeV/cm] equals the mass stopping power [MeV*cm^2/g].  Caller owns
 * the backing arrays and must keep them alive for the test. */
static void init_two_proj_tables(struct osh_material_runtime *mat_rt,
                                 unsigned int proj_z[2],
                                 unsigned int proj_a[2],
                                 double proj_mass[2],
                                 float rho_arr[1],
                                 float sp_values[4],
                                 float sp_alpha,
                                 float sp_carbon) {
    proj_z[0] = 2u;
    proj_a[0] = 4u;
    proj_mass[0] = 4.0 * OSH_PART_MASS_PROTON;
    proj_z[1] = 6u;
    proj_a[1] = 12u;
    proj_mass[1] = 12.0 * OSH_PART_MASS_PROTON;
    rho_arr[0] = 1.0f;
    /* Layout [material][projectile][energy] with nenergy == 2. */
    sp_values[0] = sp_alpha;  /* proj 0, energy 0 */
    sp_values[1] = sp_alpha;  /* proj 0, energy 1 */
    sp_values[2] = sp_carbon; /* proj 1, energy 0 */
    sp_values[3] = sp_carbon; /* proj 1, energy 1 */

    memset(mat_rt, 0, sizeof(*mat_rt));
    mat_rt->nprojectiles = 2u;
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

/* Fill a point-deposit step at Z-bin 1 (Z in [1,2)) for species with mass number a. */
static void point_step(struct step *st, double energy, double de) {
    memset(st, 0, sizeof(*st));
    st->p[2] = 1.5;
    st->p[3] = energy;
    st->q[2] = 1.5;
    st->q[3] = energy; /* point: q == p, so mean energy == birth energy */
    st->ds = 0.0;      /* zero track length: the case issue #227 handles */
    st->de = de;
    st->rho = 1.0;
    st->wt = 1.0;
    st->medium = 0;
}

/* An isolated recoil deposit makes both DLET and TLET equal its birth-energy LET
 * (S * rho), and books it at the located bin only. */
static void test_point_single_recoil_let(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dlet_page;
    struct osh_scoring_page_runtime *tlet_page;
    unsigned int proj_z[2];
    unsigned int proj_a[2];
    double proj_mass[2];
    float rho_arr[1];
    float sp_values[4];
    struct osh_material_runtime mat_rt;
    enum osh_status rc;
    float const sp_alpha = 50.0f; /* LET = 50 MeV/cm at rho = 1 */

    write_temp_file(path, sizeof(path), DETECT_TEXT);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    init_two_proj_tables(&mat_rt, proj_z, proj_a, proj_mass, rho_arr, sp_values, sp_alpha, 400.0f);
    rt.mat_tables = &mat_rt;

    /* An alpha recoil (Z=2, A=4) born at 40 MeV (10 MeV/u, inside the table). */
    memset(&part, 0, sizeof(part));
    part.charge = 2;
    part.z = 2u;
    part.a = 4u;
    part.mass = proj_mass[0];

    point_step(&st, 40.0, 4.0);
    rc = osh_scoring_score_point(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    /* Check the raw two-pass accumulators before postprocess collapses the ratio. */
    dlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DLET);
    tlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_TLET);
    ASSERT_TRUE(dlet_page != NULL);
    ASSERT_TRUE(tlet_page != NULL);
    /* DLET numerator/denominator: (LET*de, de) = (200, 4) at bin 1, zero elsewhere. */
    assert_near(dlet_page->acc.data[1], 200.0);
    assert_near(dlet_page->acc.data2[1], 4.0);
    assert_near(dlet_page->acc.data[0], 0.0);
    assert_near(dlet_page->acc.data[2], 0.0);
    /* TLET numerator/denominator: (de, de/LET) = (4, 0.08). */
    assert_near(tlet_page->acc.data[1], 4.0);
    assert_near(tlet_page->acc.data2[1], 4.0 / 50.0);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* Both averages collapse to the representative birth-energy LET. */
    assert_near(dlet_page->acc.data[1], 50.0);
    assert_near(tlet_page->acc.data[1], 50.0);
    assert_near(dlet_page->acc.data[0], 0.0);
    assert_near(tlet_page->acc.data[2], 0.0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* Two recoils of different LET sharing one bin: DLET is dose-weighted and TLET is
 * track-weighted, so the two averages differ in the physically expected direction
 * (the high-LET, short-range recoil dominates the dose average more than the track
 * average). */
static void test_point_mixed_recoils_dlet_vs_tlet(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle alpha;
    struct particle carbon;
    struct step st;
    struct osh_scoring_page_runtime *dlet_page;
    struct osh_scoring_page_runtime *tlet_page;
    unsigned int proj_z[2];
    unsigned int proj_a[2];
    double proj_mass[2];
    float rho_arr[1];
    float sp_values[4];
    struct osh_material_runtime mat_rt;
    enum osh_status rc;
    double const let_a = 100.0; /* alpha  recoil LET [MeV/cm] */
    double const let_c = 300.0; /* carbon recoil LET [MeV/cm] */
    double const de_a = 3.0;
    double const de_c = 6.0;
    double expected_dlet;
    double expected_tlet;
    double tw_a;
    double tw_c;

    write_temp_file(path, sizeof(path), DETECT_TEXT);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    init_two_proj_tables(&mat_rt, proj_z, proj_a, proj_mass, rho_arr, sp_values, (float) let_a, (float) let_c);
    rt.mat_tables = &mat_rt;

    memset(&alpha, 0, sizeof(alpha));
    alpha.charge = 2;
    alpha.z = 2u;
    alpha.a = 4u;
    alpha.mass = proj_mass[0];

    memset(&carbon, 0, sizeof(carbon));
    carbon.charge = 6;
    carbon.z = 6u;
    carbon.a = 12u;
    carbon.mass = proj_mass[1];

    /* Both recoils deposit into the same bin (Z bin 1). */
    point_step(&st, 40.0, de_a);
    rc = osh_scoring_score_point(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &alpha, &st);
    ASSERT_TRUE(rc == OSH_OK);

    point_step(&st, 120.0, de_c);
    rc = osh_scoring_score_point(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &carbon, &st);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    dlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DLET);
    tlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_TLET);
    ASSERT_TRUE(dlet_page != NULL);
    ASSERT_TRUE(tlet_page != NULL);

    expected_dlet = (let_a * de_a + let_c * de_c) / (de_a + de_c);
    tw_a = de_a / let_a;
    tw_c = de_c / let_c;
    expected_tlet = (let_a * tw_a + let_c * tw_c) / (tw_a + tw_c);

    assert_near(dlet_page->acc.data[1], expected_dlet);
    assert_near(tlet_page->acc.data[1], expected_tlet);
    /* The dose average is pulled higher than the track average by the high-LET
     * carbon recoil (sanity check on the weighting distinction). */
    ASSERT_TRUE(dlet_page->acc.data[1] > tlet_page->acc.data[1]);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* A recoil species with no SP-table column has no representable dE/dx, so it feeds
 * neither DLET nor TLET (the two-pass denominator stays zero -> ratio 0). */
static void test_point_species_not_in_table(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *dlet_page;
    struct osh_scoring_page_runtime *tlet_page;
    unsigned int proj_z[2];
    unsigned int proj_a[2];
    double proj_mass[2];
    float rho_arr[1];
    float sp_values[4];
    struct osh_material_runtime mat_rt;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), DETECT_TEXT);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* Table has Z=2 and Z=6 only. */
    init_two_proj_tables(&mat_rt, proj_z, proj_a, proj_mass, rho_arr, sp_values, 50.0f, 400.0f);
    rt.mat_tables = &mat_rt;

    /* An oxygen recoil (Z=8) is absent from the projectile columns. */
    memset(&part, 0, sizeof(part));
    part.charge = 8;
    part.z = 8u;
    part.a = 16u;
    part.mass = 16.0 * OSH_PART_MASS_PROTON;

    point_step(&st, 40.0, 4.0);
    rc = osh_scoring_score_point(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    dlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_DLET);
    tlet_page = find_page_by_kind(&rt, OSH_SCORING_SCORE_TLET);
    ASSERT_TRUE(dlet_page != NULL);
    ASSERT_TRUE(tlet_page != NULL);
    /* Nothing booked: both numerator and denominator stay zero. */
    assert_near(dlet_page->acc.data[1], 0.0);
    assert_near(dlet_page->acc.data2[1], 0.0);
    assert_near(tlet_page->acc.data[1], 0.0);
    assert_near(tlet_page->acc.data2[1], 0.0);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);
    assert_near(dlet_page->acc.data[1], 0.0);
    assert_near(tlet_page->acc.data[1], 0.0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

int main(void) {
    test_point_single_recoil_let();
    test_point_mixed_recoils_dlet_vs_tlet();
    test_point_species_not_in_table();
    printf("All osh_scoring_point_let tests passed.\n");
    return 0;
}
