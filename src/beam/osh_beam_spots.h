#ifndef _OSH_BEAM_SPOTLIST
#define _OSH_BEAM_SPOTLIST

#include "osh_beam.h"

/* get an empty beam spot list */
int osh_beam_spots_init(struct beam_spot **sl, size_t nspots);
int osh_beam_spots_free(struct beam_spot *sl);

int osh_beam_shared_init(struct beam_shared *shared);
int osh_beam_shared_free(struct beam_shared *shared);


#endif /* _OSH_BEAM_SPOTLIST */