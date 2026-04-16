#ifndef OSH_BEAM_PREPARED_H
#define OSH_BEAM_PREPARED_H

#include <stddef.h>

struct osh_beam_prepared {
    double *cum_wt; /* SOBP cumulative weights, length nspots */
    double *tm;     /* affine transforms, 16 doubles per spot */
    double wt_sum;
    double emax;
    double pmax;
    size_t nspots;
};

#endif /* OSH_BEAM_PREPARED_H */
