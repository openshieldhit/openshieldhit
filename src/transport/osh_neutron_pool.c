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

enum osh_status osh_neutron_pool_alloc(size_t capacity, struct osh_neutron_pool **out) {
    struct osh_neutron_pool *pool;
    double *slab;
    size_t slab_doubles;
    size_t off;

    if (!out || capacity == 0u) {
        return OSH_EINVAL;
    }

    /* Guard against size_t overflow before multiplying. */
    if (capacity > ((size_t) -1) / (DOUBLES_PER_ENTRY * sizeof(double))) {
        return OSH_ENOMEM;
    }

    pool = (struct osh_neutron_pool *) malloc(sizeof(*pool));
    if (!pool) {
        return OSH_ENOMEM;
    }

    slab_doubles = DOUBLES_PER_ENTRY * capacity;
    slab = (double *) malloc(slab_doubles * sizeof(double));
    if (!slab) {
        free(pool);
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

    *out = pool;
    return OSH_OK;
}

void osh_neutron_pool_free(struct osh_neutron_pool *pool) {
    if (!pool) {
        return;
    }
    /* x points to the slab base; all other arrays are offsets within it. */
    free(pool->x);
    free(pool);
}

void osh_neutron_pool_reset(struct osh_neutron_pool *pool) {
    if (!pool) {
        return;
    }
    pool->n = 0u;
    pool->n_created = 0u;
    pool->n_dropped = 0u;
}
