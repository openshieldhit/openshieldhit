#ifndef OSH_BEAM_PREPARED_H
#define OSH_BEAM_PREPARED_H

#include <stddef.h>

struct osh_beam_prepared {
    double *cum_wt; /* SOBP cumulative weights, length nspots */
    double *tm;     /* affine transforms, 16 doubles per spot */
    double wt_sum;
    double emax;
    double pmax;
    double tcut_lo; /* TCUT0 sampling-truncation lower bound, absolute MeV; 0 with tcut_hi == 0 means unset */
    double tcut_hi; /* TCUT0 sampling-truncation upper bound, absolute MeV; <= 0 means "no truncation" */
    size_t nspots;
};

#endif /* OSH_BEAM_PREPARED_H */
