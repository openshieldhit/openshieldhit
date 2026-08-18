#ifndef OSH_BEAM_PREPARED_H
#define OSH_BEAM_PREPARED_H

#include <stddef.h>

#include "random/osh_rng_gauss_trunc.h"

struct osh_beam_prepared {
    double *cum_wt; /* SOBP cumulative weights, length nspots */
    double *tm;     /* affine transforms, 16 doubles per spot */
    /* Inverse-CDF constants for the TCUT0-truncated energy draw, one entry per
     * spot, or NULL when TCUT0 is unset.  Per spot rather than per run because
     * the window (tcut_lo, tcut_hi) is shared but t0/tsigma are not, so the
     * standardised cuts -- and hence the interval constants -- differ per spot. */
    struct osh_gauss_trunc *etrunc;
    double wt_sum;
    double emax;
    double pmax;
    double tcut_lo; /* TCUT0 sampling-truncation lower bound, absolute MeV; 0 with tcut_hi == 0 means unset */
    double tcut_hi; /* TCUT0 sampling-truncation upper bound, absolute MeV; <= 0 means "no truncation" */
    size_t nspots;
};

#endif /* OSH_BEAM_PREPARED_H */
