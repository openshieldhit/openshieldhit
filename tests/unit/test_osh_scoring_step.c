#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_step.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_compile.h"
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

    assert_close(rt.pages[energy0_idx].data[0], 0.875);
    assert_close(rt.pages[energy0_idx].data[1], 1.0);
    assert_close(rt.pages[energy0_idx].data[2], 0.125);

    assert_close(rt.pages[energy1_idx].data[0], 0.875);
    assert_close(rt.pages[energy1_idx].data[1], 1.0);
    assert_close(rt.pages[energy1_idx].data[2], 0.125);

    assert_close(rt.pages[fluence_idx].data[0], 0.875);
    assert_close(rt.pages[fluence_idx].data[1], 1.0);
    assert_close(rt.pages[fluence_idx].data[2], 0.125);

    assert_close(rt.pages[filtered_idx].data[0], 0.875);
    assert_close(rt.pages[filtered_idx].data[1], 1.0);
    assert_close(rt.pages[filtered_idx].data[2], 0.125);

    st.gen = 1u;
    rc = osh_scoring_score_step(&rt, &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    assert_close(rt.pages[energy0_idx].data[0], 1.75);
    assert_close(rt.pages[energy0_idx].data[1], 2.0);
    assert_close(rt.pages[energy0_idx].data[2], 0.25);

    assert_close(rt.pages[energy1_idx].data[0], 1.75);
    assert_close(rt.pages[energy1_idx].data[1], 2.0);
    assert_close(rt.pages[energy1_idx].data[2], 0.25);

    assert_close(rt.pages[fluence_idx].data[0], 1.75);
    assert_close(rt.pages[fluence_idx].data[1], 2.0);
    assert_close(rt.pages[fluence_idx].data[2], 0.25);

    assert_close(rt.pages[filtered_idx].data[0], 0.875);
    assert_close(rt.pages[filtered_idx].data[1], 1.0);
    assert_close(rt.pages[filtered_idx].data[2], 0.125);

    for (i = 3u; i < rt.pages[energy0_idx].len; ++i) {
        assert_close(rt.pages[energy0_idx].data[i], 0.0);
        assert_close(rt.pages[energy1_idx].data[i], 0.0);
        assert_close(rt.pages[fluence_idx].data[i], 0.0);
        assert_close(rt.pages[filtered_idx].data[i], 0.0);
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

    assert_close(rt.pages[energy0_idx].data[0], 0.0);
    assert_close(rt.pages[energy0_idx].data[1], 0.875);
    assert_close(rt.pages[energy0_idx].data[2], 0.125);

    assert_close(rt.pages[energy1_idx].data[0], 0.0);
    assert_close(rt.pages[energy1_idx].data[1], 0.875);
    assert_close(rt.pages[energy1_idx].data[2], 0.125);

    assert_close(rt.pages[fluence_idx].data[0], 0.0);
    assert_close(rt.pages[fluence_idx].data[1], 3.5 * chord_len / 32.0);
    assert_close(rt.pages[fluence_idx].data[2], 0.5 * chord_len / 32.0);

    assert_close(rt.pages[filtered_idx].data[0], 0.0);
    assert_close(rt.pages[filtered_idx].data[1], 3.5 * chord_len / 32.0);
    assert_close(rt.pages[filtered_idx].data[2], 0.5 * chord_len / 32.0);

    for (i = 3u; i < rt.pages[energy0_idx].len; ++i) {
        assert_close(rt.pages[energy0_idx].data[i], 0.0);
        assert_close(rt.pages[energy1_idx].data[i], 0.0);
        assert_close(rt.pages[fluence_idx].data[i], 0.0);
        assert_close(rt.pages[filtered_idx].data[i], 0.0);
    }

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

int main(void) {
    test_score_mesh_energy_and_fluence_with_filters();
    test_score_mesh_uses_step_chord_after_bending();
    return 0;
}
