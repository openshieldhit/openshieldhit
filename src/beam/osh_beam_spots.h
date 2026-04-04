#ifndef _OSH_BEAM_SPOTLIST
#define _OSH_BEAM_SPOTLIST

#include "osh_beam.h"

/* get an empty beam spot list */
int osh_beam_spots_init(struct beam_spot **sl, size_t nspots);
int osh_beam_spots_free(struct beam_spot *sl);

/* Load an external spot list and replace beam->spots with the loaded rows.
 * Optional columns inherit defaults from the already-parsed beam->spots[0]
 * template. The loaded list replaces the template spot; it is not appended. */
int osh_beam_spotlist_load(struct beam_workspace *beam);

int osh_beam_shared_init(struct beam_shared *shared);
int osh_beam_shared_free(struct beam_shared *shared);

#endif /* _OSH_BEAM_SPOTLIST */
