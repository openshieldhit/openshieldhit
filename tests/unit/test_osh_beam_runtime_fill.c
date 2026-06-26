/*
 * Unit tests for osh_beam_runtime_fill_pool_at() — the cursor-free primary
 * generation that makes a run's history range partitionable across workers.
 *
 * The contract under test: a primary's phase space and both of its RNG streams
 * (the transient BEAM stream and the persistent PHYSICS stream carried on the
 * pool slot) are a pure function of its GLOBAL history index, supplied by the
 * caller as global_prim_base + i.  Therefore filling [0, n) in one call must be
 * byte-for-byte identical to filling [0, k) and then [k, n) in two calls — and
 * identical regardless of the order the two slices are produced.  That identity
 * is exactly what lets future workers own disjoint slices of the run.
 *
 * Coverage:
 *   test_split_matches_single   — [0,k)+[k,n) reproduces [0,n) exactly
 *   test_order_independent      — producing the high slice first changes nothing
 *   test_wrapper_advances_cursor— the serial wrapper bases off / advances the
 *                                 cursor; _at() leaves it untouched
 *   test_at_does_not_touch_rt   — _at() never mutates primaries_generated
 */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "apps/osh/osh_app_osh.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_particle_pool.h"
#include "openshieldhit/status.h"
#include "random/osh_rng.h"
#include "test_assert.h"

#define POOL_CAP 16u
#define NPRIM 10u
#define SPLIT 4u /* fill [0,SPLIT) then [SPLIT,NPRIM) */

static int _tmp_counter = 0;

static void write_temp_beam(char *path, size_t path_cap) {
    FILE *fp;
    char const *content = "PRIMARY proton\n"
                          "TMAX0 150.0 1.0\n"
                          "BEAMPOS 0.5 -0.5 -20.0\n"
                          "BEAMSIGMA 0.3 0.4\n"
                          "BEAMDIV 2.0 3.0 0.0\n";

    snprintf(path, path_cap, "osh_test_beamfill_%d.tmp", _tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

/* Build a compiled beam runtime from a freshly-written beam file. The workspace
 * is returned too so the caller can free it after the runtime. */
static void make_runtime(struct osh_beam_workspace **wb_out, struct osh_beam_runtime **rt_out) {
    char beam_path[512];
    int rc;

    write_temp_beam(beam_path, sizeof(beam_path));
    rc = (int) osh_beam_setup_from_path(beam_path, NULL, wb_out);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(*wb_out != NULL);
    ASSERT_TRUE(osh_beam_compile(*wb_out, rt_out) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void seeding_init(struct osh_rng_seeding *s) {
    s->type = OSH_RNG_TYPE_PCG32;
    s->seed = 20260626u;  /* arbitrary fixed run seed */
    s->hist_base = 1000u; /* arbitrary non-zero RNDOFFSET to catch base bugs */
}

/* Two carried PHYSICS streams are equal iff their engine state matches.  We
 * compare the live engine fields rather than memcmp the whole struct: seeding
 * leaves the gauss cache and the unused union tail uninitialised, so a raw byte
 * compare would test garbage.  The run seeds PCG32, so state+inc is the stream. */
static int rng_state_equal(struct osh_rng const *a, struct osh_rng const *b) {
    return a->type == b->type && a->u.pcg32.state == b->u.pcg32.state && a->u.pcg32.inc == b->u.pcg32.inc;
}

/* Compare every per-primary field the transport loop relies on, slot by slot. */
static void assert_slots_identical(struct osh_particle_pool const *a, struct osh_particle_pool const *b, size_t n) {
    size_t i;
    ASSERT_TRUE(a->n == b->n);
    ASSERT_TRUE(a->n >= n);
    for (i = 0; i < n; ++i) {
        ASSERT_TRUE(a->prim_idx[i] == b->prim_idx[i]);
        ASSERT_TRUE(a->gen[i] == b->gen[i]);
        ASSERT_TRUE(a->x[i] == b->x[i]);
        ASSERT_TRUE(a->y[i] == b->y[i]);
        ASSERT_TRUE(a->z[i] == b->z[i]);
        ASSERT_TRUE(a->ux[i] == b->ux[i]);
        ASSERT_TRUE(a->uy[i] == b->uy[i]);
        ASSERT_TRUE(a->uz[i] == b->uz[i]);
        ASSERT_TRUE(a->e[i] == b->e[i]);
        ASSERT_TRUE(a->wt[i] == b->wt[i]);
        /* Same global history index => same persistent PHYSICS stream. */
        ASSERT_TRUE(rng_state_equal(&a->rng[i], &b->rng[i]));
    }
}

/* ---- [0,k)+[k,n) reproduces [0,n) exactly -------------------------------- */

static void test_split_matches_single(void) {
    struct osh_beam_workspace *wb = NULL;
    struct osh_beam_runtime *rt = NULL;
    struct osh_particle_pool whole;
    struct osh_particle_pool split;
    struct osh_rng_seeding seeding;

    make_runtime(&wb, &rt);
    seeding_init(&seeding);
    ASSERT_TRUE(osh_particle_pool_init(&whole, POOL_CAP) == OSH_OK);
    ASSERT_TRUE(osh_particle_pool_init(&split, POOL_CAP) == OSH_OK);

    /* One shot over the whole range. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &whole, NPRIM, 0u) == OSH_OK);

    /* Two contiguous slices, low then high; bases carry the global index. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &split, SPLIT, 0u) == OSH_OK);
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &split, NPRIM - SPLIT, SPLIT) == OSH_OK);

    ASSERT_TRUE(whole.n == NPRIM);
    ASSERT_TRUE(split.n == NPRIM);
    assert_slots_identical(&whole, &split, NPRIM);

    osh_particle_pool_free(&whole);
    osh_particle_pool_free(&split);
    osh_beam_runtime_free(&rt);
    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
}

/* ---- producing the high slice first changes nothing ---------------------- */

static void test_order_independent(void) {
    struct osh_beam_workspace *wb = NULL;
    struct osh_beam_runtime *rt = NULL;
    struct osh_particle_pool lo_hi; /* fill [0,SPLIT) then [SPLIT,NPRIM) */
    struct osh_particle_pool hi_only;
    struct osh_rng_seeding seeding;
    size_t i;

    make_runtime(&wb, &rt);
    seeding_init(&seeding);
    ASSERT_TRUE(osh_particle_pool_init(&lo_hi, POOL_CAP) == OSH_OK);
    ASSERT_TRUE(osh_particle_pool_init(&hi_only, POOL_CAP) == OSH_OK);

    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &lo_hi, SPLIT, 0u) == OSH_OK);
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &lo_hi, NPRIM - SPLIT, SPLIT) == OSH_OK);

    /* Produce only the high slice, into a fresh pool starting at slot 0. The
     * high-slice histories must match regardless of whether the low slice was
     * ever generated or in what order — index is the only thing that matters. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &hi_only, NPRIM - SPLIT, SPLIT) == OSH_OK);

    ASSERT_TRUE(hi_only.n == NPRIM - SPLIT);
    for (i = 0; i < NPRIM - SPLIT; ++i) {
        ASSERT_TRUE(lo_hi.prim_idx[SPLIT + i] == hi_only.prim_idx[i]);
        ASSERT_TRUE(lo_hi.x[SPLIT + i] == hi_only.x[i]);
        ASSERT_TRUE(lo_hi.e[SPLIT + i] == hi_only.e[i]);
        ASSERT_TRUE(rng_state_equal(&lo_hi.rng[SPLIT + i], &hi_only.rng[i]));
    }

    osh_particle_pool_free(&lo_hi);
    osh_particle_pool_free(&hi_only);
    osh_beam_runtime_free(&rt);
    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
}

/* ---- serial wrapper bases off and advances the cursor; _at() does not ----- */

static void test_wrapper_advances_cursor(void) {
    struct osh_beam_workspace *wb = NULL;
    struct osh_beam_runtime *rt = NULL;
    struct osh_particle_pool via_wrapper;
    struct osh_particle_pool via_at;
    struct osh_rng_seeding seeding;

    make_runtime(&wb, &rt);
    seeding_init(&seeding);
    ASSERT_TRUE(osh_particle_pool_init(&via_wrapper, POOL_CAP) == OSH_OK);
    ASSERT_TRUE(osh_particle_pool_init(&via_at, POOL_CAP) == OSH_OK);

    /* Two successive wrapper calls chain through the internal cursor: the second
     * fill continues where the first left off. */
    ASSERT_TRUE(rt->primaries_generated == 0u);
    ASSERT_TRUE(osh_beam_runtime_fill_pool(rt, &seeding, &via_wrapper, SPLIT) == OSH_OK);
    ASSERT_TRUE(rt->primaries_generated == SPLIT);
    ASSERT_TRUE(osh_beam_runtime_fill_pool(rt, &seeding, &via_wrapper, NPRIM - SPLIT) == OSH_OK);
    ASSERT_TRUE(rt->primaries_generated == NPRIM);

    /* The same sequence built with explicit bases must be byte-identical, and
     * must NOT have moved the cursor. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &via_at, SPLIT, 0u) == OSH_OK);
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &via_at, NPRIM - SPLIT, SPLIT) == OSH_OK);
    ASSERT_TRUE(rt->primaries_generated == NPRIM); /* unchanged by _at() */

    assert_slots_identical(&via_wrapper, &via_at, NPRIM);

    osh_particle_pool_free(&via_wrapper);
    osh_particle_pool_free(&via_at);
    osh_beam_runtime_free(&rt);
    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
}

/* ---- _at() never mutates the cursor, even across many calls -------------- */

static void test_at_does_not_touch_rt(void) {
    struct osh_beam_workspace *wb = NULL;
    struct osh_beam_runtime *rt = NULL;
    struct osh_particle_pool pool;
    struct osh_rng_seeding seeding;

    make_runtime(&wb, &rt);
    seeding_init(&seeding);
    ASSERT_TRUE(osh_particle_pool_init(&pool, POOL_CAP) == OSH_OK);

    ASSERT_TRUE(rt->primaries_generated == 0u);
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &pool, NPRIM, 500u) == OSH_OK);
    ASSERT_TRUE(rt->primaries_generated == 0u); /* explicit base, cursor untouched */

    /* n == 0 is a well-defined no-op; NULL operands are rejected. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &pool, 0u, 0u) == OSH_OK);
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(NULL, &seeding, &pool, 1u, 0u) == OSH_EINVAL);
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, NULL, &pool, 1u, 0u) == OSH_EINVAL);
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, NULL, 1u, 0u) == OSH_EINVAL);

    osh_particle_pool_free(&pool);
    osh_beam_runtime_free(&rt);
    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
}

/* ---- a base whose range would wrap uint64_t is rejected ------------------ */

static void test_overflow_guard(void) {
    struct osh_beam_workspace *wb = NULL;
    struct osh_beam_runtime *rt = NULL;
    struct osh_particle_pool pool;
    struct osh_rng_seeding seeding;

    make_runtime(&wb, &rt);
    seeding_init(&seeding);
    ASSERT_TRUE(osh_particle_pool_init(&pool, POOL_CAP) == OSH_OK);

    /* [base, base + NPRIM) would wrap prim_idx: rejected, pool left untouched —
     * fail loudly rather than reuse history indices / RNG streams. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &pool, NPRIM, UINT64_MAX - 3u) == OSH_EINVAL);
    ASSERT_TRUE(pool.n == 0u);

    /* Even when prim_idx itself would not wrap, RNDOFFSET + prim_idx must not
     * wrap the final RNG history key. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &pool, NPRIM, UINT64_MAX - 1000u) == OSH_EINVAL);
    ASSERT_TRUE(pool.n == 0u);

    /* A large but non-wrapping final history-index range is accepted. With
     * hist_base = 1000 and NPRIM = 10, this lands exactly on UINT64_MAX. */
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &pool, NPRIM, UINT64_MAX - 1009u) == OSH_OK);
    ASSERT_TRUE(pool.n == NPRIM);
    ASSERT_TRUE(pool.prim_idx[0] == UINT64_MAX - 1009u);
    ASSERT_TRUE(pool.prim_idx[NPRIM - 1u] == UINT64_MAX - 1000u);

    /* Off-by-one boundary: a single primary at UINT64_MAX is valid when
     * RNDOFFSET is zero, because the last index is base + (n - 1). */
    seeding.hist_base = 0u;
    ASSERT_TRUE(osh_beam_runtime_fill_pool_at(rt, &seeding, &pool, 1u, UINT64_MAX) == OSH_OK);
    ASSERT_TRUE(pool.n == NPRIM + 1u);
    ASSERT_TRUE(pool.prim_idx[NPRIM] == UINT64_MAX);

    osh_particle_pool_free(&pool);
    osh_beam_runtime_free(&rt);
    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
}

int main(void) {
    test_split_matches_single();
    test_order_independent();
    test_wrapper_advances_cursor();
    test_at_does_not_touch_rt();
    test_overflow_guard();
    printf("All osh_beam_runtime_fill tests passed.\n");
    return 0;
}
