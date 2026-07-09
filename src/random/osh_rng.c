#include "random/osh_rng.h"

#include <math.h>

void osh_rng_init(struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t stream) {
    rng->type = type;
    rng->gauss_has_spare = 0;

    switch (type) {
    case OSH_RNG_TYPE_PCG32:
        osh_rng_pcg32_init(rng, seed, stream);
        break;

    case OSH_RNG_TYPE_XOSHIRO256SS:
        osh_rng_xoshiro256ss_init(rng, seed, stream);
        break;

    default:
        osh_rng_pcg32_init(rng, seed, stream);
        rng->type = OSH_RNG_TYPE_PCG32;
        break;
    }
}

/*
 * rng_mix_stream() — derive a 64-bit stream id from (seed, history index,
 * purpose).  It hash-combines the three inputs, then runs the SplitMix64
 * finaliser.  Being a pure function of its inputs is the whole point: a
 * history maps to the same stream id no matter when, on which thread, or in
 * which pool capacity it is seeded.
 *
 * Why each constant — these are standard mixing constants with known good
 * diffusion, not arbitrary values:
 *
 *   Step 1 — spread each input across all 64 bits, then XOR-combine.
 *   Each field is multiplied by a large ODD constant.  An odd constant is a
 *   unit in the ring Z/2^64, so the multiply is a bijection (it cannot lose
 *   information) and it scatters the entropy of low-magnitude inputs — a small
 *   history index, a 0/1 purpose — across the full word so they survive the
 *   XOR instead of sitting in the low bits.  The constants only need to be odd
 *   and distinct so the index and purpose axes do not alias:
 *     - 0x9E3779B97F4A7C15 = floor(2^64 / phi), the 64-bit golden-ratio
 *       constant ("Fibonacci hashing").  Consecutive indices land far apart —
 *       the classic Weyl-sequence spacing.  Refs: Knuth, TAOCP Vol. 3 sec. 6.4
 *       (multiplicative hashing); Steele, Lea & Flood, "Fast Splittable
 *       Pseudorandom Number Generators", OOPSLA 2014 (SplitMix's GOLDEN_GAMMA).
 *     - 0xD1B54A32D192ED03 spreads the (tiny) purpose enum so BEAM and PHYSICS
 *       land on well-separated, full-width offsets.  It is a recognised 64-bit
 *       mixing multiplier — the same constant used in widely-copied hash-table
 *       de-randomisation templates (the SplitMix64-based custom hashes from the
 *       competitive-programming "blowing up unordered_map" idiom).  Note it is
 *       odd but NOT prime: primality is irrelevant here; what matters is that an
 *       odd multiplier is a unit in Z/2^64 (so the multiply is invertible) and
 *       that its bit pattern diffuses well.  Any distinct large odd constant
 *       would satisfy the seeding contract; this one is a known-good choice.
 *
 *   Step 2 — the SplitMix64 finaliser (David Stafford's "Mix13" variant), an
 *   improvement on MurmurHash3's fmix64 with stronger avalanche: flipping one
 *   input bit flips ~half the output bits with low bias.  The multipliers
 *   (0xBF58476D1CE4E5B9, 0x94D049BB133111EB) and shift amounts (30, 27, 31)
 *   are taken verbatim.  Each xorshift-multiply step is invertible, so the
 *   finaliser is a bijection.  Refs: D. Stafford, "Better Bit Mixing —
 *   Improving on MurmurHash3's 64-bit Finalizer" (2011),
 *   http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html ;
 *   public-domain reference implementation https://prng.di.unimi.it/splitmix64.c
 *   This is the same finaliser used to seed xoshiro256** in
 *   osh_rng_xoshiro256ss_init(), so both engines get identical lane
 *   separation.  Disjointness over the index ranges we actually use is also
 *   checked empirically by test_seed_history_disjoint_ranges.
 */
static uint64_t rng_mix_stream(uint64_t seed, uint64_t hist_index, uint64_t purpose) {
    /* Golden-ratio constant on the index axis; a second odd constant on the
     * purpose axis (see the constant rationale above). */
    uint64_t x = seed ^ (hist_index * 0x9E3779B97F4A7C15ULL) ^ (purpose * 0xD1B54A32D192ED03ULL);

    /* Stafford Mix13 finaliser: xorshift / multiply / xorshift / multiply /
     * xorshift, giving full 64-bit avalanche. */
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

void osh_rng_seed_history(
    struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t hist_index, enum osh_rng_purpose purpose) {
    uint64_t const stream = rng_mix_stream(seed, hist_index, (uint64_t) purpose);

    osh_rng_init(rng, type, seed, stream);
}

/*
 * rng_lineage_key() — fold a parent RNG's *internal state* into one 64-bit
 * lineage key with rng_mix_stream().
 *
 * The key is derived from the engine's raw state words, never from its output.
 * That distinction is what removes the structural seed reuse P-6 flagged
 * (issue #299): an engine's output function permutes the state before emitting
 * it (PCG32's XSH-RR permutation of the pre-advance state; xoshiro256**'s
 * rotate/multiply scrambler), so the state words a split hashes here are not
 * the values the parent emits.  The child seed is therefore decoupled from the
 * parent's output window rather than — as before — taken straight from it.  A
 * chance 64-bit coincidence with some later parent draw is of course still
 * possible, but it no longer happens by construction and carries no
 * information.
 */
static uint64_t rng_lineage_key(struct osh_rng const *parent) {
    switch (parent->type) {
    case OSH_RNG_TYPE_XOSHIRO256SS: {
        /* Fold all four 64-bit state words in, each fully mixed so no bits are
         * lost: two lineages alias only if their whole 256-bit states hash
         * together, which is negligibly unlikely. */
        uint64_t key = parent->u.xoshiro256ss.s[0];

        key = rng_mix_stream(key, parent->u.xoshiro256ss.s[1], 0u);
        key = rng_mix_stream(key, parent->u.xoshiro256ss.s[2], 0u);
        key = rng_mix_stream(key, parent->u.xoshiro256ss.s[3], 0u);
        return key;
    }

    case OSH_RNG_TYPE_PCG32:
    default:
        /* Both state words matter: two PCG streams sharing state but not inc
         * are different streams, so inc must enter the key. */
        return rng_mix_stream(parent->u.pcg32.state, parent->u.pcg32.inc, 0u);
    }
}

void osh_rng_split(struct osh_rng *child, struct osh_rng const *parent, uint64_t ordinal) {
    /*
     * Seed the child by hashing the parent's *current internal state*, keyed by
     * the child's ordinal.  The parent's own stream is deliberately never
     * consumed: splitting reads parent state but does not advance it.
     *
     * This is what makes secondary seeding drop-, reorder-, and overflow-proof
     * (issue #213; design in #148).  Because splitting never draws from the
     * parent, whether a sibling secondary is injected, reordered, or silently
     * dropped when a pool is full can never shift the parent's — or any other
     * sibling's — subsequent draws.  Each child stream is a pure function of
     * its lineage key (a hash of the parent's current state, itself a pure
     * function of the parent's own draws) and its ordinal, so reproducibility
     * no longer depends on pool occupancy or wavefront scheduling.
     *
     * Unlike an earlier design that scanned a copy of the parent's *output*
     * window, the lineage key here comes from the parent's raw state words,
     * which the engine permutes before emitting — so a child seed is no longer
     * taken from, or structurally correlated with, the values the parent draws
     * next (P-6, issue #299).  Ordinal k owns the disjoint slot pair [2k, 2k+1]
     * on the hash's index axis, mirroring the old two-slot window layout, and
     * the derivation is O(1) rather than O(ordinal).
     */
    uint64_t const key = rng_lineage_key(parent);
    uint64_t const child_seed = rng_mix_stream(key, 2u * ordinal, 0u);
    uint64_t const child_stream = rng_mix_stream(key, (2u * ordinal) + 1u, 0u);

    osh_rng_init(child, parent->type, child_seed, child_stream);
}

uint32_t osh_rng_u32(struct osh_rng *rng) {
    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        return osh_rng_pcg32_u32(rng);

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return (uint32_t) (osh_rng_xoshiro256ss_u64(rng) >> 32);

    default:
        return osh_rng_pcg32_u32(rng);
    }
}

uint64_t osh_rng_u64(struct osh_rng *rng) {
    uint64_t hi;
    uint64_t lo;

    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        hi = (uint64_t) osh_rng_pcg32_u32(rng);
        lo = (uint64_t) osh_rng_pcg32_u32(rng);
        return (hi << 32) | lo;

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return osh_rng_xoshiro256ss_u64(rng);

    default:
        return ((uint64_t) osh_rng_u32(rng) << 32) | (uint64_t) osh_rng_u32(rng);
    }
}

float osh_rng_float(struct osh_rng *rng) {
    uint32_t r;
    uint32_t mant;

    /* Use 24 bits for float mantissa */
    r = osh_rng_u32(rng);
    mant = r >> 8; /* top 24 bits */

    return (float) mant * (1.0f / 16777216.0f); /* 2^24 */
}

double osh_rng_double(struct osh_rng *rng) {
    uint64_t r;
    uint64_t mant;

    /* Use 53 bits for double mantissa */
    r = osh_rng_u64(rng);
    mant = r >> 11; /* top 53 bits */

    return (double) mant * (1.0 / 9007199254740992.0); /* 2^53 */
}

double osh_rng_gauss01(struct osh_rng *rng) {
    double u;
    double v;
    double s;
    double m;

    if (rng->gauss_has_spare) {
        rng->gauss_has_spare = 0;
        return rng->gauss_spare;
    }

    /* Rejection sample in the unit disk */
    do {
        u = 2.0 * osh_rng_double(rng) - 1.0;
        v = 2.0 * osh_rng_double(rng) - 1.0;
        s = u * u + v * v;
    } while (s <= 0.0 || s >= 1.0);

    m = sqrt(-2.0 * log(s) / s);

    rng->gauss_spare = v * m;
    rng->gauss_has_spare = 1;

    return u * m;
}

double osh_rng_gauss(struct osh_rng *rng, double mu, double sigma) {
    return mu + (sigma * osh_rng_gauss01(rng));
}

/* vectors */

void osh_rng_double_vec(struct osh_rng *rng, double *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_double(rng);
        x[i + 1] = osh_rng_double(rng);
        x[i + 2] = osh_rng_double(rng);
        x[i + 3] = osh_rng_double(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_double(rng);
        i++;
    }
}

void osh_rng_float_vec(struct osh_rng *rng, float *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_float(rng);
        x[i + 1] = osh_rng_float(rng);
        x[i + 2] = osh_rng_float(rng);
        x[i + 3] = osh_rng_float(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_float(rng);
        i++;
    }
}

void osh_rng_gauss01_vec(struct osh_rng *rng, double *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_gauss01(rng);
        x[i + 1] = osh_rng_gauss01(rng);
        x[i + 2] = osh_rng_gauss01(rng);
        x[i + 3] = osh_rng_gauss01(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_gauss01(rng);
        i++;
    }
}

void osh_rng_gauss_vec(struct osh_rng *rng, double mu, double sigma, double *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = mu + sigma * osh_rng_gauss01(rng);
        x[i + 1] = mu + sigma * osh_rng_gauss01(rng);
        x[i + 2] = mu + sigma * osh_rng_gauss01(rng);
        x[i + 3] = mu + sigma * osh_rng_gauss01(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = mu + sigma * osh_rng_gauss01(rng);
        i++;
    }
}

void osh_rng_u32_vec(struct osh_rng *rng, uint32_t *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_u32(rng);
        x[i + 1] = osh_rng_u32(rng);
        x[i + 2] = osh_rng_u32(rng);
        x[i + 3] = osh_rng_u32(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_u32(rng);
        i++;
    }
}

int osh_rng_poisson(struct osh_rng *rng, double lambda) {
    double L;
    double p;
    int k;

    if (lambda <= 0.0) {
        return 0;
    }
    L = exp(-lambda);
    p = 1.0;
    k = 0;
    do {
        ++k;
        p *= osh_rng_double(rng);
    } while (p > L && k <= 64);
    return k - 1;
}
