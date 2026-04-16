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
int osh_beam_spots_init(struct beam_spot **sl, size_t nspots);

/**
 * @brief Release a spot array allocated by osh_beam_spots_init().
 *
 * @param[in] sl  Spot array to free; must not be NULL.
 *
 * @returns OSH_OK on success, OSH_EINVAL if sl is NULL.
 */
int osh_beam_spots_free(struct beam_spot *sl);

/**
 * @brief Load an external spot-list file and replace beam->spots.
 *
 * @details
 * Reads a whitespace-delimited numeric file (5, 6, 7, 9, or 11 columns per
 * row). Comment lines (starting with a character in OSH_READLINE_COMMENT) and
 * blank lines are skipped. The file is read twice: once to detect the column
 * layout and count rows, and once to fill the allocated spot array.
 *
 * Optional columns inherit their defaults from the template spot already
 * stored in beam->spots[0] by the beam.dat parser. The loaded list replaces
 * the template; it is not appended.
 *
 * @param[in,out] beam           Workspace with spots[0] initialised as a template.
 * @param[in]     spotlist_path  Absolute or already-resolved spot-list path.
 *
 * @returns OSH_OK on success, OSH_EIO on file errors, OSH_EINVAL on format
 *          errors, OSH_ENOMEM on allocation failure.
 */
int osh_beam_spotlist_load(struct beam_workspace *beam, char const *spotlist_path);

/**
 * @brief Zero-initialise a beam_shared struct and set safe default values.
 *
 * @param[in,out] shared  Struct to initialise; must not be NULL.
 *
 * @returns OSH_OK on success, OSH_EINVAL if shared is NULL.
 */
int osh_beam_shared_init(struct beam_shared *shared);

#endif /* OSH_BEAM_SPOTS_H */
