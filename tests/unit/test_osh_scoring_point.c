/*
 * Unit tests for the point scorer (osh_scoring_score_point, issue #179) on
 * cylindrical (R,Z) geometry.  A sub-threshold recoil deposits its energy at a
 * single (x,y,z) with no track length, and must land in the correct CYL voxel
 * (r_bin, z_bin) with the radius-dependent 1/V — the path previously skipped.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_step.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_pdg.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_point.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#define DETECT_PATH "osh_scoring_point_cyl_detect.tmp"

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

/*
 * R 0..2 in 2 bins (dr=1): r_bin0=[0,1), r_bin1=[1,2).
 * Z 0..2 in 2 bins (dz=1): z_bin0=[0,1), z_bin1=[1,2).
 * nr=2 -> flat idx = z_bin*2 + r_bin.
 * V(r_bin0)=pi*(1-0)*1=pi ; V(r_bin1)=pi*(4-1)*1=3pi  ->  dose(r0) = 3*dose(r1).
 */
static char const *const DETECT_TEXT = "Geometry Cyl\n"
                                       "    Name C\n"
                                       "    R 0.0 2.0 2\n"
                                       "    Z 0.0 2.0 2\n"
                                       "\n"
                                       "Output\n"
                                       "    Filename point_cyl_e.txt\n"
                                       "    Fileformat TEXT\n"
                                       "    Geo C\n"
                                       "    Quantity Energy\n"
                                       "\n"
                                       "Output\n"
                                       "    Filename point_cyl_d.txt\n"
                                       "    Fileformat TEXT\n"
                                       "    Geo C\n"
                                       "    Quantity Dose\n";

static void write_detect(void) {
    FILE *fp = fopen(DETECT_PATH, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(DETECT_TEXT, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static void build(struct osh_scoring_workspace **ws, struct osh_scoring_runtime *rt) {
    enum osh_status rc;

    write_detect();
    *ws = NULL;
    memset(rt, 0, sizeof(*rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(*ws, NULL, rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt->noutputs == 2u);
}

static void teardown(struct osh_scoring_workspace *ws, struct osh_scoring_runtime *rt) {
    osh_scoring_runtime_free(rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
}

static void point_step(struct step *st, double x, double y, double z, double de) {
    memset(st, 0, sizeof(*st));
    st->p[0] = x;
    st->p[1] = y;
    st->p[2] = z;
    st->de = de;
    st->rho = 1.0;
    st->wt = 1.0;
    st->medium = 1;
    st->zone = 1;
}

static enum osh_status deposit(struct osh_scoring_runtime *rt, struct particle const *part, struct step const *st) {
    return osh_scoring_score_point(
        rt, osh_scoring_runtime_master_accumulators(rt), osh_scoring_runtime_master_scratch(rt), part, st);
}

/* Energy lands in the located CYL voxel; out-of-bounds points deposit nothing;
 * dose reflects the radius-dependent voxel volume (inner bin 3x the outer). */
static void test_cyl_point_locate(void) {
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    struct particle proton;
    struct step st;
    size_t e_idx;
    size_t d_idx;

    build(&ws, &rt);
    memset(&proton, 0, sizeof(proton));
    proton.mass = 938.272;
    proton.pdg = OSH_PART_PDG_PROTON;
    proton.charge = 1;
    proton.z = 1u;
    proton.a = 1u;

    e_idx = rt.outputs[0].page_indices[0];
    d_idx = rt.outputs[1].page_indices[0];

    /* A: r=0.5 (bin0), z=0.5 (bin0) -> flat idx 0. */
    point_step(&st, 0.5, 0.0, 0.5, 2.0);
    ASSERT_TRUE(deposit(&rt, &proton, &st) == OSH_OK);

    /* B: r=1.5 (bin1), z=1.5 (bin1) -> flat idx 3.  y-component exercises r=sqrt(x^2+y^2). */
    point_step(&st, 0.9, 1.2, 1.5, 2.0); /* r = sqrt(0.81+1.44)=1.5 */
    ASSERT_TRUE(deposit(&rt, &proton, &st) == OSH_OK);

    /* Outside R (r=2.5 >= 2.0): no deposit. */
    point_step(&st, 2.5, 0.0, 0.5, 2.0);
    ASSERT_TRUE(deposit(&rt, &proton, &st) == OSH_OK);

    /* Below the axial stack (z=-0.5): no deposit. */
    point_step(&st, 0.5, 0.0, -0.5, 2.0);
    ASSERT_TRUE(deposit(&rt, &proton, &st) == OSH_OK);

    /* Energy: exactly de in idx 0 and idx 3, zero elsewhere (no misattribution). */
    assert_close(rt.pages[e_idx].acc.data[0], 2.0);
    assert_close(rt.pages[e_idx].acc.data[1], 0.0);
    assert_close(rt.pages[e_idx].acc.data[2], 0.0);
    assert_close(rt.pages[e_idx].acc.data[3], 2.0);

    /* Dose ~ 1/V: inner bin (V=pi) is 3x the outer bin (V=3pi), constant-free. */
    ASSERT_TRUE(rt.pages[d_idx].acc.data[0] > 0.0);
    ASSERT_TRUE(rt.pages[d_idx].acc.data[3] > 0.0);
    assert_close(rt.pages[d_idx].acc.data[0], 3.0 * rt.pages[d_idx].acc.data[3]);
    assert_close(rt.pages[d_idx].acc.data[1], 0.0);
    assert_close(rt.pages[d_idx].acc.data[2], 0.0);

    teardown(ws, &rt);
}

/* A neutral point deposit (de-excitation gamma) books energy but no dose. */
static void test_cyl_point_neutral_no_dose(void) {
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    struct particle gamma;
    struct step st;
    size_t e_idx;
    size_t d_idx;

    build(&ws, &rt);
    memset(&gamma, 0, sizeof(gamma));
    gamma.pdg = OSH_PART_PDG_GAMMA;
    gamma.charge = 0;

    e_idx = rt.outputs[0].page_indices[0];
    d_idx = rt.outputs[1].page_indices[0];

    point_step(&st, 0.5, 0.0, 0.5, 2.0);
    ASSERT_TRUE(deposit(&rt, &gamma, &st) == OSH_OK);

    assert_close(rt.pages[e_idx].acc.data[0], 2.0);
    assert_close(rt.pages[d_idx].acc.data[0], 0.0);

    teardown(ws, &rt);
}

int main(void) {
    test_cyl_point_locate();
    test_cyl_point_neutral_no_dose();
    printf("All osh_scoring_point tests passed.\n");
    return 0;
}
