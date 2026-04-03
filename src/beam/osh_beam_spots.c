#include "beam/osh_beam_spots.h"

#include <stdlib.h>
#include <string.h>

#include "common/osh_rc.h"

int osh_beam_spots_init(struct beam_spot **sl, size_t nspots) {
    if (!sl || nspots == 0) {
        return OSH_EINVAL;
    }

    *sl = calloc(1, sizeof(struct beam_spot));
    if (!*sl) {
        return OSH_ENOMEM;
    }

    return OSH_OK;
}

int osh_beam_spots_free(struct beam_spot *sl) {
    if (!sl) {
        return OSH_EINVAL;
    }
    free(sl);
    return OSH_OK;
}

int osh_beam_shared_init(struct beam_shared *shared) {
    if (!shared) {
        return OSH_EINVAL;
    }
    memset(shared, 0, sizeof *shared);

    /* emax/pmax are derived after parse — leave at 0.0 until computed */
    shared->emax = 0.0;
    shared->pmax = 0.0;
    shared->sad[0] = 0.0;
    shared->sad[1] = 0.0;
    shared->focus = 0.0;
    shared->theta = 0.0; /* rad — along +Z by default */
    shared->phi = 0.0;   /* rad */
    shared->use_div = 0;
    shared->use_sad = 0;

    return OSH_OK;
}
