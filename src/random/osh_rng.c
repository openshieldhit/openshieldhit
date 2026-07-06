#include "random/osh_rng.h"
#include "random/osh_rng_hd.h"

#include <math.h>

/* Body lives in osh_rng_hd.h so device kernels compile the same lines. */
void osh_rng_init(struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t stream) {
    _osh_rng_init_hd(rng, type, seed, stream);
}

void osh_rng_seed_history(
    struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t hist_index, enum osh_rng_purpose purpose) {
    _osh_rng_seed_history_hd(rng, type, seed, hist_index, purpose);
}

void osh_rng_split(struct osh_rng *child, struct osh_rng const *parent, uint64_t ordinal) {
    /*
     * Seed the child from a *private copy* of the parent advanced to the
     * child's ordinal slot.  The parent's own stream is deliberately never
     * consumed: splitting reads parent state but does not advance it.
     *
     * This is what makes secondary seeding drop-, reorder-, and overflow-proof
     * (issue #213; design in #148).  Because splitting no longer draws from the
     * parent, whether a sibling secondary is injected, reordered, or silently
     * dropped when a pool is full can never shift the parent's — or any other
     * sibling's — subsequent draws.  Each child stream is a pure function of
     * its lineage key (the parent's current state, itself a pure function of
     * the parent's own draws) and its ordinal, so reproducibility no longer
     * depends on pool occupancy or wavefront scheduling.
     *
     * Ordinal k owns the disjoint two-draw window [2k, 2k+1] of the copied
     * stream — the exact layout the old sequential split produced — so an
     * injected child keeps the identical stream it had before, while the parent
     * is left untouched.  The scan cost is O(ordinal); secondary counts per
     * event are bounded (OSH_NUCLEAR_MAX_SECONDARIES), so this is negligible.
     */
    struct osh_rng scan = *parent; /* copy: parent is const and stays put */
    uint64_t child_seed;
    uint64_t child_stream;
    uint64_t i;

    for (i = 0u; i < ordinal; ++i) {
        (void) osh_rng_u64(&scan);
        (void) osh_rng_u64(&scan);
    }
    child_seed = osh_rng_u64(&scan);
    child_stream = osh_rng_u64(&scan);

    osh_rng_init(child, parent->type, child_seed, child_stream);
}

uint32_t osh_rng_u32(struct osh_rng *rng) {
    return _osh_rng_u32_hd(rng);
}

uint64_t osh_rng_u64(struct osh_rng *rng) {
    return _osh_rng_u64_hd(rng);
}

float osh_rng_float(struct osh_rng *rng) {
    return _osh_rng_float_hd(rng);
}

double osh_rng_double(struct osh_rng *rng) {
    return _osh_rng_double_hd(rng);
}

double osh_rng_gauss01(struct osh_rng *rng) {
    return _osh_rng_gauss01_hd(rng);
}

double osh_rng_gauss(struct osh_rng *rng, double mu, double sigma) {
    return _osh_rng_gauss_hd(rng, mu, sigma);
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
