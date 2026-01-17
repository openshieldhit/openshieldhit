#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "osh_rng.h"

#define ASSERT_TRUE(cond)                                                      \
    do {                                                                       \
        if (!(cond)) {                                                         \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__,    \
                    __LINE__);                                                 \
            exit(1);                                                           \
        }                                                                      \
    } while (0)

static void test_pcg32_known_sequence(void) {
    struct osh_rng r;
    osh_rng_init(&r, OSH_RNG_TYPE_PCG32, 42u, 54u);

    /* Expected sequence for (seed=42, stream=54) with current implementation */
    const uint32_t exp[5] = {
        2707161783u, 2068313097u, 3122475824u, 2211639955u, 3215226955u,
    };

    for (int i = 0; i < 5; ++i) {
        uint32_t x = osh_rng_u32(&r);
        ASSERT_TRUE(x == exp[i]);
    }
}

static void test_xoshiro256ss_known_sequence(void) {
    struct osh_rng r;
    osh_rng_init(&r, OSH_RNG_TYPE_XOSHIRO256SS, 42u, 54u);

    const uint64_t exp[5] = {
        9619421891339311063ull,  17143628181114060176ull,
        17740981343507171333ull, 7781542089684599863ull,
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
    const double exp[4] = {
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

int main(void) {
    test_pcg32_known_sequence();
    test_xoshiro256ss_known_sequence();
    test_uniform_ranges();
    test_gauss01_known_values();

    return 0;
}
