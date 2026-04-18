#ifndef OSH_BEAM_SPOTS_H
#define OSH_BEAM_SPOTS_H

#include "osh_beam.h"

/**
 * @brief Allocate an array of nspots zero-initialised beam_spot structs.
 *
 * @param[out] sl      Receives the allocated array; caller must free with
 *                     osh_beam_spots_free().
 * @param[in]  nspots  Number of spots to allocate; must be > 0.
 *
 * @returns OSH_OK on success, OSH_EINVAL if sl is NULL or nspots is 0,
 *          OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_beam_spots_init(struct beam_spot **sl, size_t nspots);

/**
 * @brief Release a spot array allocated by osh_beam_spots_init().
 *
 * @param[in] sl  Spot array to free; may be NULL.
 */
void osh_beam_spots_free(struct beam_spot *sl);

/**
 * @brief Zero-initialise a beam_shared struct and set safe default values.
 *
 * @param[in,out] shared  Struct to initialise; must not be NULL.
 *
 * @returns OSH_OK on success, OSH_EINVAL if shared is NULL.
 */
int osh_beam_shared_init(struct beam_shared *shared);

#endif /* OSH_BEAM_SPOTS_H */
