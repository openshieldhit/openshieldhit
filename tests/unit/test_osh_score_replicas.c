#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/simulation.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "test_assert.h"

/*
 * Unit coverage for the sequential score-replica harness (issue #230): the
 * contiguous integer partition, the private-set clone helpers, the merge reduce
 * that folds them into the master, and the set_score_replicas API guard.  The
 * full end-to-end reproducibility contract (serial == N=1 byte-for-byte, etc.) is
 * exercised by tests/cases/run_score_replicas.cmake.
 */

/* Mirror of the contiguous integer partition used by run_score_replicas() in
 * src/transport/osh_transport.c: sub-range r is [nstat*r/N, nstat*(r+1)/N). */
static size_t part_lo(size_t nstat, size_t n, size_t r) {
    return (size_t) ((unsigned long long) nstat * r / n);
}

static void test_partition_tiles_range(void) {
    size_t const nstats[] = {1u, 2u, 7u, 100u, 999u, 1000u};
    size_t const ns[] = {1u, 2u, 3u, 4u, 7u, 100u};
    size_t si;
    size_t ni;

    for (si = 0u; si < sizeof(nstats) / sizeof(nstats[0]); ++si) {
        for (ni = 0u; ni < sizeof(ns) / sizeof(ns[0]); ++ni) {
            size_t const nstat = nstats[si];
            size_t const n = ns[ni];
            size_t r;

            if (n > nstat) {
                continue; /* N > nstat is rejected upstream, never partitioned */
            }
            /* First sub-range starts at 0; last ends exactly at nstat; every
             * sub-range is non-empty and consecutive ones share no history. */
            ASSERT_TRUE(part_lo(nstat, n, 0u) == 0u);
            ASSERT_TRUE(part_lo(nstat, n, n) == nstat);
            for (r = 0u; r < n; ++r) {
                size_t const lo = part_lo(nstat, n, r);
                size_t const hi = part_lo(nstat, n, r + 1u);
                ASSERT_TRUE(hi > lo);     /* non-empty */
                ASSERT_TRUE(hi <= nstat); /* in range */
            }
        }
    }
}

/* Compile the shared test01 scoring runtime; caller frees ws + rt. */
static void compile_runtime(struct osh_scoring_workspace **ws_out, struct osh_scoring_runtime *rt) {
    char path[512];

    snprintf(path, sizeof(path), "%s/test01/detect.dat", OSH_TEST_FIXTURES_DIR);
    *ws_out = NULL;
    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, ws_out) == OSH_OK);
    ASSERT_TRUE(*ws_out != NULL);
    memset(rt, 0, sizeof(*rt));
    ASSERT_TRUE(osh_scoring_compile(*ws_out, NULL, rt) == OSH_OK);
    ASSERT_TRUE(rt->npages > 0u);
}

static void test_clone_accumulators_and_scratch_match_master(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct osh_scoring_accumulator *set = NULL;
    struct osh_scoring_scratch scratch;
    size_t p;

    compile_runtime(&ws, &rt);

    /* A shaped set is npages long and matches each master page exactly. */
    set = (struct osh_scoring_accumulator *) calloc(rt.npages, sizeof(*set));
    ASSERT_TRUE(set != NULL);
    ASSERT_TRUE(osh_scoring_runtime_alloc_accumulator_set(&rt, set) == OSH_OK);
    for (p = 0u; p < rt.npages; ++p) {
        ASSERT_TRUE(set[p].data != NULL);
        ASSERT_TRUE(set[p].len == rt.pages[p].acc.len);
        ASSERT_TRUE((set[p].data2 != NULL) == (rt.pages[p].acc.data2 != NULL));
        ASSERT_TRUE(set[p].weight == 0.0);
        ASSERT_TRUE(set[p].nbatch == 0u);
    }

    /* The cloned scratch mirrors the master's pre-sized crossing buffer. */
    memset(&scratch, 0, sizeof(scratch));
    ASSERT_TRUE(osh_scoring_runtime_clone_scratch(&rt, &scratch) == OSH_OK);
    ASSERT_TRUE(scratch.crossing_cap == rt.master_scratch.crossing_cap);
    if (rt.master_scratch.crossing_cap > 0u) {
        ASSERT_TRUE(scratch.crossing_buf != NULL);
    }

    osh_scoring_runtime_free_scratch(&scratch);
    osh_scoring_runtime_free_accumulator_set(set, rt.npages);
    free(set);
    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/*
 * The reduce at the heart of the harness: two private sets folded into the empty
 * master must equal their element-wise sum (data and data2), matching what a
 * single depositor would have accumulated.  This is the invariant N==1 relies on
 * (merge into empty master == identity) and that N>1 relies on (order-independent
 * additive combine).
 */
static void test_private_sets_merge_into_master(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct osh_scoring_accumulator *s0 = NULL;
    struct osh_scoring_accumulator *s1 = NULL;
    struct osh_scoring_accumulator *master;
    size_t p;
    size_t i;

    compile_runtime(&ws, &rt);
    s0 = (struct osh_scoring_accumulator *) calloc(rt.npages, sizeof(*s0));
    s1 = (struct osh_scoring_accumulator *) calloc(rt.npages, sizeof(*s1));
    ASSERT_TRUE(s0 != NULL && s1 != NULL);
    ASSERT_TRUE(osh_scoring_runtime_alloc_accumulator_set(&rt, s0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_runtime_alloc_accumulator_set(&rt, s1) == OSH_OK);

    /* Seed each private set with distinct, index-dependent values. */
    for (p = 0u; p < rt.npages; ++p) {
        for (i = 0u; i < s0[p].len; ++i) {
            s0[p].data[i] = (double) (p + 1u) + 0.5 * (double) i;
            s1[p].data[i] = 100.0 * (double) (p + 1u) - 0.25 * (double) i;
            if (s0[p].data2) {
                s0[p].data2[i] = 2.0 * (double) i + 1.0;
                s1[p].data2[i] = 3.0 * (double) i + 7.0;
            }
        }
    }

    /* Merge both into the master (zeroed by compile): 0 + s0 + s1. */
    master = osh_scoring_runtime_master_accumulators(&rt);
    ASSERT_TRUE(master != NULL);
    for (p = 0u; p < rt.npages; ++p) {
        ASSERT_TRUE(osh_scoring_accumulator_merge(&master[p], &s0[p]) == OSH_OK);
        ASSERT_TRUE(osh_scoring_accumulator_merge(&master[p], &s1[p]) == OSH_OK);
    }

    for (p = 0u; p < rt.npages; ++p) {
        for (i = 0u; i < master[p].len; ++i) {
            double const expect =
                ((double) (p + 1u) + 0.5 * (double) i) + (100.0 * (double) (p + 1u) - 0.25 * (double) i);
            ASSERT_TRUE(master[p].data[i] == expect); /* exact: additive of the two seeded sets */
            if (master[p].data2) {
                double const expect2 = (2.0 * (double) i + 1.0) + (3.0 * (double) i + 7.0);
                ASSERT_TRUE(master[p].data2[i] == expect2);
            }
        }
    }

    osh_scoring_runtime_free_accumulator_set(s0, rt.npages);
    osh_scoring_runtime_free_accumulator_set(s1, rt.npages);
    free(s0);
    free(s1);
    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* set_score_replicas accepts 0..nstat and rejects N > nstat and a NULL sim. */
static void test_set_score_replicas_validates(void) {
    char geo_path[512];
    char beam_path[512];
    char mat_path[512];
    char detect_path[512];
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/00_minimal/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(beam_path, sizeof(beam_path), "%s/tests/cases/00_minimal/beam.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/00_minimal/mat.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(detect_path, sizeof(detect_path), "%s/tests/cases/00_minimal/detect.dat", OSH_PROJECT_SOURCE_DIR);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(detect_path, NULL, &scoring) == OSH_OK);

    beam->nstat = 10u;
    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);

    ASSERT_TRUE(osh_simulation_set_score_replicas(NULL, 1u) == OSH_EINVAL);
    ASSERT_TRUE(osh_simulation_set_score_replicas(sim, 11u) == OSH_EINVAL); /* > nstat */
    ASSERT_TRUE(osh_simulation_set_score_replicas(sim, 10u) == OSH_OK);     /* == nstat (boundary) */
    ASSERT_TRUE(osh_simulation_set_score_replicas(sim, 1u) == OSH_OK);
    ASSERT_TRUE(osh_simulation_set_score_replicas(sim, 0u) == OSH_OK); /* 0 disables the harness */

    osh_simulation_free(sim);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
}

int main(void) {
    test_partition_tiles_range();
    test_clone_accumulators_and_scratch_match_master();
    test_private_sets_merge_into_master();
    test_set_score_replicas_validates();
    return 0;
}
