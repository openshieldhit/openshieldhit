#ifndef OSH_APP_OSH_BEAM_SPOTLIST_H
#define OSH_APP_OSH_BEAM_SPOTLIST_H

#include <stddef.h>

#include "openshieldhit/beam.h"
#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"

/**
 * @brief Import a legacy USECBEAM spot list from text into cold beam spots.
 *
 * @details
 * This is an app-layer helper for the `USECBEAM` parser path. The file format
 * is whitespace-delimited numeric text with one spot per row and one shared
 * column layout across the file. Missing optional columns inherit defaults
 * from @p template_spot.
 *
 * The returned spot array is owned by the caller and should be released with
 * `free()` after it has been copied into a beam workspace.
 */
enum osh_status osh_beam_spotlist_import(char const *path,
                                         struct osh_diag_sink const *diag,
                                         struct osh_beam_spot const *template_spot,
                                         struct osh_beam_spot **spots_out,
                                         size_t *nspots_out);

#endif /* OSH_APP_OSH_BEAM_SPOTLIST_H */
