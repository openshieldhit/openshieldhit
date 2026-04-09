#ifndef OSH_SCORING_GEOMETRY_RUNTIME_H
#define OSH_SCORING_GEOMETRY_RUNTIME_H

#include <stddef.h>

#include "scoring/runtime/osh_scoring_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_scoring_axis_runtime {
    double lo;     /* Lower bound [cm]. */
    double hi;     /* Upper bound [cm]. */
    int nbins;     /* Number of bins; negative means unset. */
    char label[4]; /* "X", "Y", "Z", "R" — null-terminated. */
};

struct osh_scoring_geometry_score_group {
    size_t first_page;
    size_t npages;
    enum osh_scoring_score_kind score_kind;
};

struct osh_scoring_geometry_runtime {
    char *kind;
    char *name;
    struct osh_scoring_axis_runtime *axes;
    struct osh_scoring_geometry_score_group *groups;
    size_t naxes;
    double rot_theta_deg;
    double rot_phi_deg;
    size_t nbins;
    size_t first_page;
    size_t npages;
    size_t ngroups;
    enum osh_scoring_geo_kind geo_kind;
    int zone_start;
    int zone_stop;
    char has_rotation;
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_GEOMETRY_RUNTIME_H */
