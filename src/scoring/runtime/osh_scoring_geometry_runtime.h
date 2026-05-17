#ifndef OSH_SCORING_GEOMETRY_RUNTIME_H
#define OSH_SCORING_GEOMETRY_RUNTIME_H

#include <stddef.h>

#include "common/raytrace/osh_raytrace.h"
#include "scoring/runtime/osh_scoring_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief One compiled axis of a scoring geometry.
 */
struct osh_scoring_axis_runtime {
    double lo;     /* Lower bound [cm]. */
    double hi;     /* Upper bound [cm]. */
    int nbins;     /* Number of bins. */
    char label[4]; /* "X", "Y", "Z", or "R" — null-terminated. */
};

/**
 * @brief A contiguous span of pages sharing the same geometry and score kind.
 *
 * @details
 * Pages are laid out so all pages of the same geometry are contiguous.
 * Groups subdivide that span by score kind so the hot path can loop over
 * pages of one kind without branching on kind per step.
 */
struct osh_scoring_geometry_score_group {
    size_t first_page;                      /* Index of first page in this group. */
    size_t npages;                          /* Number of pages in this group. */
    enum osh_scoring_score_kind score_kind; /* Shared score kind for all pages. */
};

/**
 * @brief Compiled scoring geometry, ready for hot-path bin lookup.
 *
 * @details
 * Holds the resolved geometry kind, axis bounds, bin counts, and the
 * contiguous page range owned by this geometry in the flat page array.
 * Groups further subdivide the page span by score kind.
 */
struct osh_scoring_geometry_runtime {
    char *kind;                                      /* Geometry type keyword (lowercase, owned). */
    char *name;                                      /* User-visible name (owned). */
    struct osh_scoring_axis_runtime *axes;           /* Axis array (owned). */
    struct osh_scoring_geometry_score_group *groups; /* Score-kind groups (owned). */
    size_t naxes;                                    /* Number of axes. */
    double rot_theta_deg;                            /* Rotation polar angle [deg]. */
    double rot_phi_deg;                              /* Rotation azimuth angle [deg]. */
    size_t nbins;                                    /* Total number of bins (product of axes). */
    size_t first_page;                               /* Index of first page in the flat page array. */
    size_t npages;                                   /* Number of pages owned by this geometry. */
    size_t ngroups;                                  /* Number of score-kind groups. */
    enum osh_scoring_geo_kind geo_kind;              /* Resolved geometry kind enum. */
    int zone_start;                                  /* First zone (Zone geometry only). */
    int zone_stop;                                   /* Last zone  (Zone geometry only). */
    char has_rotation;                               /* Non-zero when rot_* fields are set. */
    struct osh_raytrace_grid vox_grid;               /* DicomCT/DicomRTDOSE: grid geometry (origin, spacing, n). */
    char *rtdose_template_path; /* Non-NULL when FileFormat RTDOSE; owned; path to RTDOSE template. */
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_GEOMETRY_RUNTIME_H */
