#include "common/osh_particle_pool.h"

#include <stdlib.h>
#include <string.h>

/*
 * Single-slab allocation layout.
 *
 * All arrays are carved from one malloc() call so that pool_free() is a
 * single free().  The slab is laid out in field order:
 *
 *   [ x | y | z | ux | uy | uz | e ]  — 7 * capacity doubles
 *   [ wt ]                             — 1 * capacity doubles
 *   [ prim_idx ]                       — capacity uint32_t
 *   [ gen ]                            — capacity uint8_t
 *   [ species ]                        — capacity pointers
 *
 * Each segment is aligned to the natural alignment of its element type by
 * padding the preceding segment to a multiple of that alignment.  For
 * simplicity the slab is sized conservatively: all segments are treated as
 * double-sized (8 bytes), which satisfies every element's alignment.
 */

/*
 * Number of double-sized slots occupied by the non-double tail fields per
 * entry, rounded up for alignment:
 *   prim_idx : sizeof(uint32_t)            = 4 bytes → 1 double covers 2 entries
 *   gen      : sizeof(uint8_t)             = 1 byte  → 1 double covers 8 entries
 *   species  : sizeof(struct particle **)  = 8 bytes → 1 double per entry
 *
 * Conservatively: per entry, tail needs sizeof(uint32_t) + sizeof(uint8_t)
 * + sizeof(void*) bytes.  We round the whole tail up to doubles.
 */
#define TAIL_BYTES_PER_ENTRY (sizeof(uint32_t) + sizeof(uint8_t) + sizeof(struct particle const *))
#define TAIL_DOUBLES_PER_ENTRY ((TAIL_BYTES_PER_ENTRY + sizeof(double) - 1u) / sizeof(double))

/* Total double-sized slots needed per entry: 8 phase-space+weight + tail. */
#define DOUBLES_PER_ENTRY (8u + TAIL_DOUBLES_PER_ENTRY)

enum osh_status osh_particle_pool_alloc(size_t capacity, struct osh_particle_pool **pool_out) {
    struct osh_particle_pool *pool;
    double *slab;
    size_t slab_doubles;
    size_t off;

    if (!pool_out || capacity == 0u) {
        return OSH_EINVAL;
    }

    pool = (struct osh_particle_pool *) malloc(sizeof(*pool));
    if (!pool) {
        return OSH_ENOMEM;
    }

    slab_doubles = DOUBLES_PER_ENTRY * capacity;
    slab = (double *) malloc(slab_doubles * sizeof(double));
    if (!slab) {
        free(pool);
        return OSH_ENOMEM;
    }

    /* Partition the slab. */
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
     * Tail fields (prim_idx, gen, species) are placed in the remaining slab
     * space after the 8 double arrays.  Each pointer is derived from the
     * slab base plus a byte offset, cast to the appropriate type.  The slab
     * allocation above ensures sufficient space and natural alignment for
     * double*, so uint32_t* and uint8_t* are always within bounds.
     */
    {
        unsigned char *tail = (unsigned char *) (slab + off);
        size_t tail_off = 0u;

        pool->prim_idx = (uint32_t *) (tail + tail_off);
        tail_off += capacity * sizeof(uint32_t);

        /* Align gen to its natural alignment (1 byte — no padding needed). */
        pool->gen = (uint8_t *) (tail + tail_off);
        tail_off += capacity * sizeof(uint8_t);

        /* Align species pointer array to pointer alignment. */
        tail_off = (tail_off + sizeof(void *) - 1u) & ~(sizeof(void *) - 1u);
        pool->species = (struct particle const **) (tail + tail_off);
    }

    pool->n = 0u;
    pool->capacity = capacity;

    *pool_out = pool;
    return OSH_OK;
}

void osh_particle_pool_free(struct osh_particle_pool *pool) {
    if (!pool) {
        return;
    }
    /* x points to the start of the slab; all other arrays are offsets within it. */
    free(pool->x);
    free(pool);
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
            pool->species[dst] = pool->species[src];
        }
        ++dst;
    }
    pool->n = dst;
}
