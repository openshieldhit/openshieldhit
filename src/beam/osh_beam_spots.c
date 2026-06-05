#include "beam/osh_beam_spots.h"

#include <stdlib.h>
#include <string.h>

#include "openshieldhit/status.h"

enum osh_status osh_beam_spots_init(struct beam_spot **sl, size_t nspots) {
    if (!sl || nspots == 0) {
        return OSH_EINVAL;
    }

    *sl = calloc(nspots, sizeof(struct beam_spot));
    if (!*sl) {
        return OSH_ENOMEM;
    }

    return OSH_OK;
}

void osh_beam_spots_free(struct beam_spot *sl) {
    free(sl);
}

int osh_beam_shared_init(struct beam_shared *shared) {
    if (!shared) {
        return OSH_EINVAL;
    }
    memset(shared, 0, sizeof *shared);

    shared->sad[0] = 0.0;
    shared->sad[1] = 0.0;
    shared->focus = 0.0;
    shared->theta = 0.0; /* rad — along +Z by default */
    shared->phi = 0.0;   /* rad */
    shared->use_div = 0;
    shared->use_sad = 0;
    shared->sad_was_set = 0;

    return OSH_OK;
}
