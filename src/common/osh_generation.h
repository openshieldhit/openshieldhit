#ifndef OSH_GENERATION_H
#define OSH_GENERATION_H

#include <stdint.h>

/*
 * Particle generation is stored as uint8_t in transport pools.  Children use
 * parent + 1 until the storage limit is reached, then saturate at this value
 * instead of wrapping to zero.
 */
#define OSH_GENERATION_MAX UINT8_MAX

static inline uint8_t osh_generation_child(uint8_t parent_gen) {
    uint8_t child_gen;

    if (parent_gen < OSH_GENERATION_MAX) {
        child_gen = (uint8_t) (parent_gen + 1u);
    } else {
        child_gen = OSH_GENERATION_MAX;
    }
    return child_gen;
}

#endif /* OSH_GENERATION_H */
