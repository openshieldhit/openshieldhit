#include <string.h>

#include "osh_rc.h"
#include "osh_beam_spots.h"


int osh_beam_spots_init(struct beam_spot **sl, size_t nspots) {
    if (!sl || nspots == 0) return OSH_EINVAL;

    *sl = calloc(1, sizeof(struct beam_spot));
    if (!*sl) return OSH_ENOMEM;

    return OSH_OK;
}


int osh_beam_spots_free(struct beam_spot *sl) {
    if (!sl) return OSH_EINVAL;
    free(sl);
    return OSH_OK;
}


int osh_beam_shared_init(struct beam_shared *shared) {
    if (!shared) return OSH_EINVAL;
    memset(shared, 0, sizeof *shared);
    shared->tcut[0] = 0.0;
    shared->tcut[1] = 2500.0; /* MeV */
    shared->pcut[0] = 0.0;
    shared->pcut[1] = 2500.0; /* MeV/c */
    shared->sad[0] = 0.0;
    shared->sad[1] = 0.0; /* cm */
    shared->focus = 0.0; /* cm */
    shared->theta = 0.0; /* rad */
    shared->phi = 0.0;   /* rad */
    shared->use_pmax = 0;
    shared->use_psigma = 0;
    shared->use_div = 0;
    shared->use_sad = 0;
    return OSH_OK;

}


int osh_beam_shared_free(struct beam_shared *shared) {
    if (!shared) return OSH_EINVAL;
    free(shared);
    return OSH_OK;
}