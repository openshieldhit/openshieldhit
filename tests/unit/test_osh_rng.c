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

/* Disjoint history-index ranges (e.g. one per MPI rank via a disjoint
 * global_prim_base) must not alias.  Two streams alias only if they are
 * seeded to the *same generator state*, so we compare the seeded PCG32 state
 * (state, inc) pairs rather than a single output: distinct streams can
 * legitimately share a first draw, but they must never share full state. */
static void test_seed_history_disjoint_ranges(void) {
    uint64_t state[HIST_N];
    uint64_t inc[HIST_N];
    uint64_t h;
    uint64_t i;
    uint64_t j;

    for (h = 0u; h < HIST_N; ++h) {
        struct osh_rng r;
        osh_rng_seed_history(&r, OSH_RNG_TYPE_PCG32, 9u, 1000000u + h, OSH_RNG_PURPOSE_PHYSICS);
        state[h] = r.u.pcg32.state;
        inc[h] = r.u.pcg32.inc;
    }
    for (i = 0u; i < HIST_N; ++i) {
        for (j = i + 1u; j < HIST_N; ++j) {
            ASSERT_TRUE(state[i] != state[j] || inc[i] != inc[j]);
        }
    }
}

/* Regression lock for issue #317: two runs sharing an RNDSEED but using
 * distinct RNDOFFSET values must decorrelate streams across the *entire*
 * shared [0, HIST_N) history-index range -- not just a non-overlapping tail.
 * Before the fix, RNDOFFSET shifted the history index instead of the seed, so
 * two offsets one apart (e.g. -N 0 and -N 1) produced byte-identical streams
 * for every index they both cover, and only the disjoint tail differed.
 *
 * Goes through osh_rng_seeding_init() -- the real production fold -- rather
 * than hand-computing the two seeds, so a regression in the fold itself (not
 * just in rng_mix_stream) would be caught here too. */
static void test_seed_history_independent_seedoffsets(void) {
    struct osh_rng_seeding seeding_a;
    struct osh_rng_seeding seeding_b;
    uint64_t h;

    osh_rng_seeding_init(&seeding_a, OSH_RNG_TYPE_PCG32, 42u, 0u);
    osh_rng_seeding_init(&seeding_b, OSH_RNG_TYPE_PCG32, 42u, 1u);

    for (h = 0u; h < HIST_N; ++h) {
        struct osh_rng ra;
        struct osh_rng rb;

        osh_rng_seed_history(&ra, OSH_RNG_TYPE_PCG32, seeding_a.seed, h, OSH_RNG_PURPOSE_PHYSICS);
        osh_rng_seed_history(&rb, OSH_RNG_TYPE_PCG32, seeding_b.seed, h, OSH_RNG_PURPOSE_PHYSICS);
        ASSERT_TRUE(ra.u.pcg32.state != rb.u.pcg32.state || ra.u.pcg32.inc != rb.u.pcg32.inc);
    }
}

/* Regression lock for the RNDSEED/RNDOFFSET aliasing pitfall found reviewing
 * the #317 fix: folding RNDOFFSET into the seed by plain addition made
 * (RNDSEED=S, RNDOFFSET=k) indistinguishable from (RNDSEED=S+k,
 * RNDOFFSET=0) -- two different, independently-chosen configurations
 * silently sharing one stream family, which double-counts histories if their
 * outputs are ever merged as independent replicas. osh_rng_seeding_init()
 * must decorrelate these via a hash-mix instead of letting the sums collide. */
static void test_seedoffset_does_not_alias_rndseed(void) {
    struct osh_rng_seeding seeding_a;
    struct osh_rng_seeding seeding_b;

    osh_rng_seeding_init(&seeding_a, OSH_RNG_TYPE_PCG32, 100u, 1u);
    osh_rng_seeding_init(&seeding_b, OSH_RNG_TYPE_PCG32, 101u, 0u);

    ASSERT_TRUE(seeding_a.seed != seeding_b.seed);
}

/* osh_rng_split is deterministic along a lineage: identical parents at the
 * same draw position produce identical, reproducible child streams for a given
 * ordinal.  This is what makes nuclear secondaries reproducible irrespective of
 * scheduling. */
static void test_split_deterministic(void) {
    struct osh_rng pa;
    struct osh_rng pb;
    struct osh_rng ca;
    struct osh_rng cb;
    int k;

    osh_rng_seed_history(&pa, OSH_RNG_TYPE_PCG32, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);
    osh_rng_seed_history(&pb, OSH_RNG_TYPE_PCG32, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);

    /* Advance both parents identically, then split at the same lineage point
     * and ordinal. */
    for (k = 0; k < 5; ++k) {
        (void) osh_rng_u64(&pa);
        (void) osh_rng_u64(&pb);
    }
    osh_rng_split(&ca, &pa, 0u);
    osh_rng_split(&cb, &pb, 0u);

    for (k = 0; k < 8; ++k) {
        ASSERT_TRUE(osh_rng_u64(&ca) == osh_rng_u64(&cb));
    }

    /* The child stream is distinct from the parent's own continuation. */
    {
        struct osh_rng cc;
        int differ = 0;
        osh_rng_seed_history(&pa, OSH_RNG_TYPE_PCG32, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);
        osh_rng_split(&cc, &pa, 0u);
        for (k = 0; k < 8; ++k) {
            if (osh_rng_u64(&cc) != osh_rng_u64(&pa)) {
                differ = 1;
            }
        }
        ASSERT_TRUE(differ);
    }
}

/*
 * Splitting must not consume from the parent, and each ordinal must key an
 * independent child.  Together these make secondary seeding drop-, reorder-,
 * and overflow-proof (issue #213): whether a sibling secondary is injected or
 * silently dropped cannot shift the parent's stream or any other sibling's.
 */
static void test_split_nonconsuming_and_drop_independent(void) {
    struct osh_rng parent;
    struct osh_rng snapshot;
    struct osh_rng pristine;
    struct osh_rng c0;
    struct osh_rng c1;
    struct osh_rng c2;
    int k;

    osh_rng_seed_history(&parent, OSH_RNG_TYPE_PCG32, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);
    for (k = 0; k < 5; ++k) {
        (void) osh_rng_u64(&parent);
    }
    snapshot = parent; /* consumed by step 1's lockstep comparison  */
    pristine = parent; /* untouched copy of the pre-split state for step 3 */

    osh_rng_split(&c0, &parent, 0u);
    osh_rng_split(&c1, &parent, 1u);
    osh_rng_split(&c2, &parent, 2u);

    /* (1) Non-consuming: the parent's own draws are exactly what they would have
     *     been had it never been split. */
    for (k = 0; k < 16; ++k) {
        ASSERT_TRUE(osh_rng_u64(&parent) == osh_rng_u64(&snapshot));
    }

    /* (2) Distinct ordinals key distinct streams (c0 and c2 are preserved for
     *     step 3 by drawing from copies here). */
    {
        struct osh_rng a = c0;
        struct osh_rng b = c1;
        struct osh_rng bb = c1;
        struct osh_rng cc = c2;
        int differ01 = 0;
        int differ12 = 0;
        for (k = 0; k < 8; ++k) {
            if (osh_rng_u64(&a) != osh_rng_u64(&b)) {
                differ01 = 1;
            }
            if (osh_rng_u64(&bb) != osh_rng_u64(&cc)) {
                differ12 = 1;
            }
        }
        ASSERT_TRUE(differ01);
        ASSERT_TRUE(differ12);
    }

    /* (3) Drop-independence: the ordinal-2 child is identical whether or not
     *     ordinals 0 and 1 were ever split.  A dropped/absent sibling cannot
     *     change another secondary's stream. */
    {
        struct osh_rng c2_alone;
        osh_rng_split(&c2_alone, &pristine, 2u);
        for (k = 0; k < 8; ++k) {
            ASSERT_TRUE(osh_rng_u64(&c2_alone) == osh_rng_u64(&c2));
        }
    }
}

/*
 * P-6 / issue #299 regression: a split child must be seeded from the parent's
 * internal *state*, not from the parent's forthcoming output window.  The old
 * design fed the parent's very next two u64 draws into the child as
 * (seed, stream) for ordinal 0; reconstruct that window-scan child and assert
 * the current split produces a different stream, so the child seed is no longer
 * structurally derived from the parent's next draws.  Runs for both engines,
 * whose lineage-key derivation differs, and re-checks the non-advancing
 * guarantee.
 */
static void test_split_not_seeded_from_parent_output(void) {
    enum osh_rng_type const engines[2] = {OSH_RNG_TYPE_PCG32, OSH_RNG_TYPE_XOSHIRO256SS};

    int e;

    for (e = 0; e < 2; ++e) {
        struct osh_rng parent;
        struct osh_rng before;
        struct osh_rng scan;
        struct osh_rng old_style;
        struct osh_rng child;
        uint64_t window_seed;
        uint64_t window_stream;
        int differ = 0;
        int k;

        osh_rng_seed_history(&parent, engines[e], 2025u, 7u, OSH_RNG_PURPOSE_PHYSICS);
        before = parent; /* snapshot: the split must never advance the parent */

        /* The parent's next two u64 draws are exactly what the old window-scan
         * split fed into osh_rng_init as (seed, stream) for ordinal 0. */
        scan = parent;
        window_seed = osh_rng_u64(&scan);
        window_stream = osh_rng_u64(&scan);
        osh_rng_init(&old_style, engines[e], window_seed, window_stream);

        osh_rng_split(&child, &parent, 0u);

        /* Non-advancing across both engines. */
        ASSERT_TRUE(osh_rng_u64(&parent) == osh_rng_u64(&before));

        /* The child no longer reuses the parent's forthcoming output window. */
        for (k = 0; k < 8; ++k) {
            if (osh_rng_u64(&child) != osh_rng_u64(&old_style)) {
                differ = 1;
            }
        }
        ASSERT_TRUE(differ);
    }
}

/*
 * The split's core guarantees must hold for every engine, not only PCG32 —
 * osh_rng_split derives its lineage key per engine type.  This mirrors
 * test_split_deterministic / _drop_independent for xoshiro256**.
 */
static void test_split_xoshiro_properties(void) {
    struct osh_rng parent;
    struct osh_rng pristine;
    struct osh_rng c0;
    struct osh_rng c1;
    struct osh_rng c2;
    int k;

    osh_rng_seed_history(&parent, OSH_RNG_TYPE_XOSHIRO256SS, 55u, 7u, OSH_RNG_PURPOSE_PHYSICS);
    for (k = 0; k < 5; ++k) {
        (void) osh_rng_u64(&parent);
    }
    pristine = parent; /* untouched copy of the pre-split state */

    osh_rng_split(&c0, &parent, 0u);
    osh_rng_split(&c1, &parent, 1u);
    osh_rng_split(&c2, &parent, 2u);

    /* Distinct ordinals key distinct streams. */
    {
        struct osh_rng a = c0;
        struct osh_rng b = c1;
        int differ = 0;
        for (k = 0; k < 8; ++k) {
            if (osh_rng_u64(&a) != osh_rng_u64(&b)) {
                differ = 1;
            }
        }
        ASSERT_TRUE(differ);
    }

    /* Drop-independence: the ordinal-2 child is identical whether or not
     * ordinals 0 and 1 were ever split (split reads but never advances). */
    {
        struct osh_rng c2_alone;
        osh_rng_split(&c2_alone, &pristine, 2u);
        for (k = 0; k < 8; ++k) {
            ASSERT_TRUE(osh_rng_u64(&c2_alone) == osh_rng_u64(&c2));
        }
    }
}

/* ---- Vector helpers (osh_rng_*_vec) ------------------------------------- */

/*
 * The *_vec helpers must produce exactly the same sequence as repeated scalar
 * calls — they are just unrolled batch wrappers.  Seeding two RNGs identically
 * and comparing a batch fill against element-by-element scalar draws locks
 * that contract (and exercises the non-multiple-of-4 tail of each loop).
 */
#define VEC_N 10 /* not a multiple of 4 — covers the unrolled body and the tail */

static void test_double_vec_matches_scalar(void) {
    struct osh_rng a;
    struct osh_rng b;
    double vec[VEC_N];
    int i;

    osh_rng_init(&a, OSH_RNG_TYPE_PCG32, 7u, 3u);
    osh_rng_init(&b, OSH_RNG_TYPE_PCG32, 7u, 3u);

    osh_rng_double_vec(&a, vec, VEC_N);
    for (i = 0; i < VEC_N; ++i) {
        ASSERT_TRUE(vec[i] == osh_rng_double(&b));
        ASSERT_TRUE(vec[i] >= 0.0 && vec[i] < 1.0);
    }
}

static void test_float_vec_matches_scalar(void) {
    struct osh_rng a;
    struct osh_rng b;
    float vec[VEC_N];
    int i;

    osh_rng_init(&a, OSH_RNG_TYPE_PCG32, 11u, 5u);
    osh_rng_init(&b, OSH_RNG_TYPE_PCG32, 11u, 5u);

    osh_rng_float_vec(&a, vec, VEC_N);
    for (i = 0; i < VEC_N; ++i) {
        ASSERT_TRUE(vec[i] == osh_rng_float(&b));
        ASSERT_TRUE(vec[i] >= 0.0f && vec[i] < 1.0f);
    }
}

static void test_u32_vec_matches_scalar(void) {
    struct osh_rng a;
    struct osh_rng b;
    uint32_t vec[VEC_N];
    int i;

    osh_rng_init(&a, OSH_RNG_TYPE_XOSHIRO256SS, 13u, 1u);
    osh_rng_init(&b, OSH_RNG_TYPE_XOSHIRO256SS, 13u, 1u);

    osh_rng_u32_vec(&a, vec, VEC_N);
    for (i = 0; i < VEC_N; ++i) {
        ASSERT_TRUE(vec[i] == osh_rng_u32(&b));
    }
}

static void test_gauss_vecs_match_scalar(void) {
    struct osh_rng a;
    struct osh_rng b;
    double g01[VEC_N];
    double g[VEC_N];
    double const mu = 2.5;
    double const sigma = 0.75;
    int i;

    osh_rng_init(&a, OSH_RNG_TYPE_PCG32, 21u, 2u);
    osh_rng_init(&b, OSH_RNG_TYPE_PCG32, 21u, 2u);
    osh_rng_gauss01_vec(&a, g01, VEC_N);
    for (i = 0; i < VEC_N; ++i) {
        ASSERT_TRUE(g01[i] == osh_rng_gauss01(&b));
    }

    /* gauss_vec is gauss01_vec scaled by (mu, sigma). */
    osh_rng_init(&a, OSH_RNG_TYPE_PCG32, 21u, 2u);
    osh_rng_init(&b, OSH_RNG_TYPE_PCG32, 21u, 2u);
    osh_rng_gauss_vec(&a, mu, sigma, g, VEC_N);
    for (i = 0; i < VEC_N; ++i) {
        ASSERT_TRUE(fabs(g[i] - (mu + sigma * osh_rng_gauss01(&b))) < 1e-15);
    }
}

static void test_poisson(void) {
    struct osh_rng r;
    long count;
    int i;
    double mean;
    int const samples = 20000;
    double const lambda = 3.0;

    osh_rng_init(&r, OSH_RNG_TYPE_PCG32, 99u, 1u);

    /* Non-positive lambda is defined to return 0. */
    ASSERT_TRUE(osh_rng_poisson(&r, 0.0) == 0);
    ASSERT_TRUE(osh_rng_poisson(&r, -1.0) == 0);

    count = 0;
    for (i = 0; i < samples; ++i) {
        int k = osh_rng_poisson(&r, lambda);
        ASSERT_TRUE(k >= 0);
        count += k;
    }
    mean = (double) count / (double) samples;
    /* Sample mean should be close to lambda (loose bound, ~5 sigma). */
    ASSERT_TRUE(fabs(mean - lambda) < 0.1);
}

int main(void) {
    test_pcg32_known_sequence();
    test_xoshiro256ss_known_sequence();
    test_uniform_ranges();
    test_gauss01_known_values();
    test_seed_history_order_independence();
    test_seed_history_purpose_independence();
    test_seed_history_disjoint_ranges();
    test_seed_history_independent_seedoffsets();
    test_seedoffset_does_not_alias_rndseed();
    test_split_deterministic();
    test_split_nonconsuming_and_drop_independent();
    test_split_not_seeded_from_parent_output();
    test_split_xoshiro_properties();
    test_double_vec_matches_scalar();
    test_float_vec_matches_scalar();
    test_u32_vec_matches_scalar();
    test_gauss_vecs_match_scalar();
    test_poisson();

    return 0;
}
