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
 * A group is not a set of spatial bins or ray crossings.  It is a run of output
 * pages that all refer to the same geometry and the same score kind, e.g. several
 * DOSE pages with different filters, Settings overrides, or differential axes.
 *
 * Pages are laid out so all pages of the same geometry are contiguous.  Groups
 * subdivide that span by score kind so the hot path can select the estimator once
 * and then loop only over compatible pages.
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
 *
 * ### Axis content by geometry kind
 *
 * @c axes[] carries axis descriptors identified by their @c label field.
 * The array order reflects declaration order in @c detect.dat and must
 * not be assumed to be fixed.  Always look up by label:
 *
 * | geo_kind              | axes present           | labels         |
 * |-----------------------|------------------------|----------------|
 * | OSH_SCORING_GEO_MESH  | 3 (X, Y, Z)            | "X", "Y", "Z" |
 * | OSH_SCORING_GEO_CYL   | 2 (radial + axial)     | "R", "Z"      |
 * | OSH_SCORING_GEO_ZONE  | 0 (zone list instead)  | —             |
 *
 * For @c CYL the @c axes array contains exactly one entry labelled @c "R"
 * and one labelled @c "Z"; which comes first depends on the input file.
 * Code that needs nr or nz must scan by label, not by index.
 */
struct osh_scoring_geometry_runtime {
    char *kind;                                      /* Geometry type keyword (lowercase, owned). */
    char *name;                                      /* User-visible name (owned). */
    struct osh_scoring_axis_runtime *axes;           /* Axis array (owned); see label-lookup note above. */
    struct osh_scoring_geometry_score_group *groups; /* Score-kind groups (owned). */
    size_t naxes;                                    /* Number of axes. */
    double t[16];                                    /* Universe→local affine transform (row-major 4×4). */
    size_t nbins;                                    /* Total number of bins (product of axes). */
    size_t first_page;                               /* Index of first page in the flat page array. */
    size_t npages;                                   /* Number of pages owned by this geometry. */
    size_t ngroups;                                  /* Number of score-kind groups. */
    enum osh_scoring_geo_kind geo_kind;              /* Resolved geometry kind enum. */
    size_t *zone_indices;                            /* Internal GEMCA zone indices, one per Zone bin (owned). */
    double *zone_vol_inv;                            /* Zone geometry only: 1/volume [1/cm3], length nzone_indices. */
    size_t nzone_indices;                            /* Number of explicit Zone bins. */
    char has_rotation;                               /* Non-zero when t[] is valid and axes are in local frame. */
    char *rtdose_template_path; /* Non-NULL when FileFormat RTDOSE; owned; path to RTDOSE template. */
    double *cyl_vol_inv;        /* [cyl_nr] 1/V per R-bin, precomputed at compile; NULL for non-CYL */
    size_t cyl_nr;              /* R bin count (nr); 0 for non-CYL */
    double *bin_vol_inv;        /* [nbins] 1/volume per spatial bin, precomputed at compile. The single
                                   geometry-agnostic volume source for volume-normalised estimators; applied
                                   in postprocess (not at score time). Uniform for Mesh, per-R for Cyl,
                                   per-zone for Zone. NULL when the geometry has no meaningful volume. */
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_GEOMETRY_RUNTIME_H */
