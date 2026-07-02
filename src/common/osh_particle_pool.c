#include "common/osh_particle_pool.h"

#include <stdlib.h>
#include <string.h>

#include "random/osh_rng.h"

/*
 * Single-slab allocation layout.
 *
 * All arrays are carved from one malloc() call so that pool_free() is a
 * single free().  The slab is laid out in field order:
 *
 *   [ x | y | z | ux | uy | uz | e ]  — 7 * capacity doubles
 *   [ wt ]                             — 1 * capacity doubles
 *   [ prim_idx ]                       — capacity uint64_t
 *   [ rng ]                            — capacity struct osh_rng
 *   [ species ]                        — capacity pointers
 *   [ gen ]                            — capacity uint8_t
 *
 * Each segment is aligned to the natural alignment of its element type by
 * padding the preceding segment to a multiple of that alignment.  The 8-byte
 * fields (prim_idx, rng, species) are placed first so each starts on an
 * 8-aligned offset (capacity * 8k bytes); the 1-byte gen field comes last.
 */

/*
 * Bytes occupied by the non-double tail fields per entry, then rounded up to
 * whole doubles so the slab is conservatively sized:
 *   prim_idx : sizeof(uint64_t)            = 8 bytes
 *   rng      : sizeof(struct osh_rng)      — per-slot RNG state
 *   species  : sizeof(struct particle **)  = 8 bytes
 *   gen      : sizeof(uint8_t)             = 1 byte
 *
 * The round-up to whole doubles also absorbs any per-segment alignment
 * padding, so the partition below never runs past the allocated slab.
 */
#define TAIL_BYTES_PER_ENTRY                                                                                           \
    (sizeof(uint64_t) + sizeof(struct osh_rng) + sizeof(struct particle const *) + sizeof(uint8_t))
#define TAIL_DOUBLES_PER_ENTRY ((TAIL_BYTES_PER_ENTRY + sizeof(double) - 1u) / sizeof(double))

/* Total double-sized slots needed per entry: 8 phase-space+weight + tail. */
#define DOUBLES_PER_ENTRY (8u + TAIL_DOUBLES_PER_ENTRY)

/* Shared slab-partition logic used by both _init and _alloc. */
static enum osh_status particle_pool_alloc_slab(struct osh_particle_pool *pool, size_t capacity) {
    double *slab;
    size_t slab_doubles;
    size_t off;

    /* Guard the slab size product against size_t overflow. */
    if (capacity > ((size_t) -1) / (DOUBLES_PER_ENTRY * sizeof(double))) {
        return OSH_ENOMEM;
    }

    slab_doubles = DOUBLES_PER_ENTRY * capacity;
    slab = (double *) malloc(slab_doubles * sizeof(double));
    if (!slab) {
        return OSH_ENOMEM;
    }

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

    /*
     * Tail fields after the 8 double arrays.  8-byte-aligned fields first so
     * offsets (multiples of capacity * 8) need no extra padding; gen (1-byte)
     * is placed last.  The slab base is double-aligned, satisfying uint64_t,
     * struct osh_rng, and void* alignment.
     */
    {
        unsigned char *tail = (unsigned char *) (slab + off);
        size_t tail_off = 0u;

        pool->prim_idx = (uint64_t *) (tail + tail_off);
        tail_off += capacity * sizeof(uint64_t);

        tail_off = (tail_off + _Alignof(struct osh_rng) - 1u) & ~(_Alignof(struct osh_rng) - 1u);
        pool->rng = (struct osh_rng *) (tail + tail_off);
        tail_off += capacity * sizeof(struct osh_rng);

        tail_off = (tail_off + sizeof(void *) - 1u) & ~(sizeof(void *) - 1u);
        pool->species = (struct particle const **) (tail + tail_off);
        tail_off += capacity * sizeof(struct particle const *);

        pool->gen = (uint8_t *) (tail + tail_off);
    }

    pool->n = 0u;
    pool->capacity = capacity;
    pool->n_dropped = 0u;
    return OSH_OK;
}

enum osh_status osh_particle_pool_init(struct osh_particle_pool *pool, size_t capacity) {
    if (!pool || capacity == 0u) {
        return OSH_EINVAL;
    }
    memset(pool, 0, sizeof(*pool));
    return particle_pool_alloc_slab(pool, capacity);
}

enum osh_status osh_particle_pool_alloc(size_t capacity, struct osh_particle_pool **pool_out) {
    struct osh_particle_pool *pool;
    enum osh_status rc;

    if (!pool_out || capacity == 0u) {
        return OSH_EINVAL;
    }

    pool = (struct osh_particle_pool *) malloc(sizeof(*pool));
    if (!pool) {
        return OSH_ENOMEM;
    }

    rc = particle_pool_alloc_slab(pool, capacity);
    if (rc != OSH_OK) {
        free(pool);
        return rc;
    }

    *pool_out = pool;
    return OSH_OK;
}

void osh_particle_pool_free(struct osh_particle_pool *pool) {
    if (!pool) {
        return;
    }
    /* x points to the slab base; free it.  The pool struct itself is not freed
     * here — it is owned by the caller (typically embedded in osh_simulation). */
    free(pool->x);
    memset(pool, 0, sizeof(*pool));
}

void osh_particle_pool_compact(struct osh_particle_pool *pool) {
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
            pool->rng[dst] = pool->rng[src]; /* carry the slot's RNG stream */
            pool->species[dst] = pool->species[src];
        }
        ++dst;
    }
    pool->n = dst;
}
