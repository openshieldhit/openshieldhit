#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "random/osh_rng.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_pcg32_known_sequence(void) {
    struct osh_rng r;
    osh_rng_init(&r, OSH_RNG_TYPE_PCG32, 42u, 54u);

    /* Expected sequence for (seed=42, stream=54) with current implementation */
    uint32_t const exp[5] = {
        2707161783u,
        2068313097u,
        3122475824u,
        2211639955u,
        3215226955u,
    };

    for (int i = 0; i < 5; ++i) {
        uint32_t x = osh_rng_u32(&r);
        ASSERT_TRUE(x == exp[i]);
    }
}

static void test_xoshiro256ss_known_sequence(void) {
    struct osh_rng r;
    osh_rng_init(&r, OSH_RNG_TYPE_XOSHIRO256SS, 42u, 54u);

    uint64_t const exp[5] = {
        9619421891339311063ull,
        17143628181114060176ull,
        17740981343507171333ull,
        7781542089684599863ull,
        309587622474537745ull,
    };

    for (int i = 0; i < 5; ++i) {
        uint64_t x = osh_rng_u64(&r);
        ASSERT_TRUE(x == exp[i]);
    }
}

static void test_uniform_ranges(void) {
    struct osh_rng r;
    osh_rng_init(&r, OSH_RNG_TYPE_PCG32, 1u, 2u);

    for (int i = 0; i < 10000; ++i) {
        double d = osh_rng_double(&r);
        float f = osh_rng_float(&r);

        ASSERT_TRUE(d >= 0.0 && d < 1.0);
        ASSERT_TRUE(f >= 0.0f && f < 1.0f);
    }
}

static void test_gauss01_known_values(void) {
    struct osh_rng r;
    osh_rng_init(&r, OSH_RNG_TYPE_PCG32, 42u, 54u);

    /* Expected first values (seed=42, stream=54) */
    double const exp[4] = {
        0.8010234011838121,
        1.3954298012510542,
        0.83712788683267014,
        0.83888355725713559,
    };

    for (int i = 0; i < 4; ++i) {
        double x = osh_rng_gauss01(&r);
        /* tight tolerance is fine; these are deterministic */
        ASSERT_TRUE(fabs(x - exp[i]) < 1e-15);
    }
}

/* ---- Per-history seeding (osh_rng_seed_history / osh_rng_split) ---------- */

/*
 * The seeding primitives back the transport parallelization guarantee: a
 * history's stream depends only on its key, never on execution order, so the
 * same history sees the same draws regardless of pool capacity, thread, or
 * rank.  These tests lock the properties the transport loop relies on.
 */

#define HIST_N 256
#define DRAWS 4

/* Stream is a pure function of (seed, hist_index, purpose): seeding histories
 * in any order yields the exact same per-history sequence.  This is the
 * order-independence the wavefront loop depends on. */
static void test_seed_history_order_independence(void) {
    uint64_t ref[HIST_N][DRAWS];
    uint64_t h;
    int k;

    for (h = 0u; h < HIST_N; ++h) {
        struct osh_rng r;
        osh_rng_seed_history(&r, OSH_RNG_TYPE_PCG32, 12345u, h, OSH_RNG_PURPOSE_PHYSICS);
        for (k = 0; k < DRAWS; ++k) {
            ref[h][k] = osh_rng_u64(&r);
        }
    }

    /* Re-seed in reverse order — every history must reproduce its sequence. */
    for (h = HIST_N; h-- > 0u;) {
        struct osh_rng r;
        osh_rng_seed_history(&r, OSH_RNG_TYPE_PCG32, 12345u, h, OSH_RNG_PURPOSE_PHYSICS);
        for (k = 0; k < DRAWS; ++k) {
            ASSERT_TRUE(osh_rng_u64(&r) == ref[h][k]);
        }
    }
}

/* The same history index drives independent streams under different purposes,
 * so beam sampling stays decoupled from in-transport physics. */
static void test_seed_history_purpose_independence(void) {
    struct osh_rng beam;
    struct osh_rng phys;
    int differ = 0;
    int k;

    osh_rng_seed_history(&beam, OSH_RNG_TYPE_PCG32, 777u, 42u, OSH_RNG_PURPOSE_BEAM);
    osh_rng_seed_history(&phys, OSH_RNG_TYPE_PCG32, 777u, 42u, OSH_RNG_PURPOSE_PHYSICS);

    for (k = 0; k < 8; ++k) {
        if (osh_rng_u64(&beam) != osh_rng_u64(&phys)) {
            differ = 1;
        }
    }
    ASSERT_TRUE(differ);
}

/* Disjoint history-index ranges (e.g. one per MPI rank via hist_base) yield
 * disjoint streams: no first-draw collisions across a large index span. */
static void test_seed_history_disjoint_ranges(void) {
    uint64_t first[HIST_N];
    uint64_t h;
    uint64_t i;
    uint64_t j;

    for (h = 0u; h < HIST_N; ++h) {
        struct osh_rng r;
        osh_rng_seed_history(&r, OSH_RNG_TYPE_PCG32, 9u, 1000000u + h, OSH_RNG_PURPOSE_PHYSICS);
        first[h] = osh_rng_u64(&r);
    }
    for (i = 0u; i < HIST_N; ++i) {
        for (j = i + 1u; j < HIST_N; ++j) {
            ASSERT_TRUE(first[i] != first[j]);
        }
    }
}

/* osh_rng_split is deterministic along a lineage: identical parents at the
 * same draw position produce identical, reproducible child streams.  This is
 * what makes nuclear secondaries reproducible irrespective of scheduling. */
static void test_split_deterministic(void) {
    struct osh_rng pa;
    struct osh_rng pb;
    struct osh_rng ca;
    struct osh_rng cb;
    int k;

    osh_rng_seed_history(&pa, OSH_RNG_TYPE_PCG32, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);
    osh_rng_seed_history(&pb, OSH_RNG_TYPE_PCG32, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);

    /* Advance both parents identically, then split at the same lineage point. */
    for (k = 0; k < 5; ++k) {
        (void) osh_rng_u64(&pa);
        (void) osh_rng_u64(&pb);
    }
    osh_rng_split(&ca, &pa);
    osh_rng_split(&cb, &pb);

    for (k = 0; k < 8; ++k) {
        ASSERT_TRUE(osh_rng_u64(&ca) == osh_rng_u64(&cb));
    }

    /* The child stream is distinct from the parent's continuation. */
    {
        struct osh_rng cc;
        int differ = 0;
        osh_rng_seed_history(&pa, OSH_RNG_TYPE_PCG32, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);
        osh_rng_split(&cc, &pa);
        for (k = 0; k < 8; ++k) {
            if (osh_rng_u64(&cc) != osh_rng_u64(&pa)) {
                differ = 1;
            }
        }
        ASSERT_TRUE(differ);
    }
}

int main(void) {
    test_pcg32_known_sequence();
    test_xoshiro256ss_known_sequence();
    test_uniform_ranges();
    test_gauss01_known_values();
    test_seed_history_order_independence();
    test_seed_history_purpose_independence();
    test_seed_history_disjoint_ranges();
    test_split_deterministic();

    return 0;
}
