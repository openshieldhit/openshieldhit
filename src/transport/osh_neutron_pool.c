#include "transport/osh_neutron_pool.h"

#include <stdlib.h>
#include <string.h>

#include "random/osh_rng.h"

/*
 * Single-slab allocation layout (mirrors osh_particle_pool, no species array).
 *
 *   [ x | y | z | ux | uy | uz | e | wt ]  — 8 * capacity doubles
 *   [ prim_idx ]                             — capacity uint64_t
 *   [ rng ]                                  — capacity struct osh_rng
 *   [ gen ]                                  — capacity uint8_t
 *
 * 8-byte-aligned fields (prim_idx, rng) are placed before the 1-byte gen field
 * so no per-segment padding is needed beyond the conservative round-up below.
 */

#define TAIL_BYTES_PER_ENTRY (sizeof(uint64_t) + sizeof(struct osh_rng) + sizeof(uint8_t))
#define TAIL_DOUBLES_PER_ENTRY ((TAIL_BYTES_PER_ENTRY + sizeof(double) - 1u) / sizeof(double))

/* 8 phase-space+weight doubles + conservative tail */
#define DOUBLES_PER_ENTRY (8u + TAIL_DOUBLES_PER_ENTRY)

enum osh_status osh_neutron_pool_init(struct osh_neutron_pool *pool, size_t capacity) {
    double *slab;
    size_t slab_doubles;
    size_t off;

    if (!pool || capacity == 0u) {
        return OSH_EINVAL;
    }

    /* Guard against size_t overflow before multiplying. */
    if (capacity > ((size_t) -1) / (DOUBLES_PER_ENTRY * sizeof(double))) {
        return OSH_ENOMEM;
    }

    memset(pool, 0, sizeof(*pool));

    slab_doubles = DOUBLES_PER_ENTRY * capacity;
    slab = (double *) malloc(slab_doubles * sizeof(double));
    if (!slab) {
        return OSH_ENOMEM;
    }

    /* Partition the double region. */
    off = 0u;
    pool->x = slab + off;
    off += capacity;
    pool->y = slab + off;
    off += capacity;
    pool->z = slab + off;
    off += capacity;
    pool->ux = slab + off;
    off += capacity;
    pool->uy = slab + off;
    off += capacity;
    pool->uz = slab + off;
    off += capacity;
    pool->e = slab + off;
    off += capacity;
    pool->wt = slab + off;
    off += capacity;

    /* Tail fields carved from the remaining slab bytes. */
    {
        unsigned char *tail;
        size_t tail_off;

        tail = (unsigned char *) (slab + off);
        tail_off = 0u;

        pool->prim_idx = (uint64_t *) (tail + tail_off);
        tail_off += capacity * sizeof(uint64_t);

        tail_off = (tail_off + _Alignof(struct osh_rng) - 1u) & ~(_Alignof(struct osh_rng) - 1u);
        pool->rng = (struct osh_rng *) (tail + tail_off);
        tail_off += capacity * sizeof(struct osh_rng);

        /* gen has 1-byte alignment — no padding needed. */
        pool->gen = (uint8_t *) (tail + tail_off);
    }

    pool->n = 0u;
    pool->capacity = capacity;
    pool->n_created = 0u;
    pool->n_dropped = 0u;

    return OSH_OK;
}

void osh_neutron_pool_free(struct osh_neutron_pool *pool) {
    if (!pool) {
        return;
    }
    free(pool->x);
    memset(pool, 0, sizeof(*pool));
}

void osh_neutron_pool_reset(struct osh_neutron_pool *pool) {
    if (!pool) {
        return;
    }
    pool->n = 0u;
    pool->n_created = 0u;
    pool->n_dropped = 0u;
}

void osh_neutron_pool_compact(struct osh_neutron_pool *pool) {
    size_t dst;
    size_t src;

    if (!pool || pool->n == 0u) {
        return;
    }

    dst = 0u;
    for (src = 0u; src < pool->n; ++src) {
        if (pool->e[src] <= 0.0) {
            continue; /* dead — skip */
        }
        if (dst != src) {
            pool->x[dst] = pool->x[src];
            pool->y[dst] = pool->y[src];
            pool->z[dst] = pool->z[src];
            pool->ux[dst] = pool->ux[src];
            pool->uy[dst] = pool->uy[src];
            pool->uz[dst] = pool->uz[src];
            pool->e[dst] = pool->e[src];
            pool->wt[dst] = pool->wt[src];
            pool->prim_idx[dst] = pool->prim_idx[src];
            pool->gen[dst] = pool->gen[src];
            pool->rng[dst] = pool->rng[src];
        }
        ++dst;
    }
    pool->n = dst;
}
