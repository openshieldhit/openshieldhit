#include "common/osh_interpolate.h"

#include <math.h>

#include "common/osh_logger.h"

double osh_interpolate_flin(float x, float const *xx, float const *ff, unsigned int n, int mode) {
    float x0 = 0, x1 = 0;
    float y0 = 0, y1 = 0;
    long int i = 0;
    float y = 0;

    if (x < xx[0]) { /* exceeds lower bound */
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            osh_error("osh_interpolate: lower bound exceeded.");
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            /* The case x == xx[0] or x == xx[n-1] may seem to be lost here, but it is not:
               in this case, OSH_INTERPOLATE_OOB_NEAREST will give the same as OSH_INTERPOLATE_OOB_EXTRAPOL.
             */
            return ff[0];
        default:
            /* by default, assume OSH_INTERPOLATE_OOB_EXTRAPOL */
            break;
        }
        x0 = xx[0];
        x1 = xx[1];
        y0 = ff[0];
        y1 = ff[1];
    } else if (x > xx[n - 1]) { /* exceeds upper bound */
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            osh_error("osh_interpolate: upper bound exceeded.");
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            /* See comment above */
            return ff[n - 1];
        default:
            /* by default, do an extrapolate */
            break;
        }
        x0 = xx[n - 2];
        x1 = xx[n - 1];
        y0 = ff[n - 2];
        y1 = ff[n - 1];
    } else {
        /* binary search here, return index i to lower element*/
        i = osh_binary_search_f(x, xx, n);
        if (i < 0) {
            osh_error("osh_interpolate: bad binary_search. This should not happen.");
            return NAN;
        }
        /* the highest possible number i can achieve is len - 2. */
        x0 = xx[i];
        x1 = xx[i + 1];
        y0 = ff[i];
        y1 = ff[i + 1];
    }
    y = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    return y;
}

double osh_interpolate_dlin(double x, double const *xx, double const *ff, unsigned int n, int mode) {
    double x0 = 0, x1 = 0;
    double y0 = 0, y1 = 0;
    long int i = 0;
    double y = 0;

    if (x < xx[0]) { /* exceeds lower bound */
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            osh_error("osh_interpolate: lower bound exceeded.");
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            /* The case x == xx[0] or x == xx[n-1] may seem to be lost here, but it is not:
               in this case, OSH_INTERPOLATE_OOB_NEAREST will give the same as OSH_INTERPOLATE_OOB_EXTRAPOL.
             */
            return ff[0];
        default:
            /* by default, assume OSH_INTERPOLATE_OOB_EXTRAPOL */
            break;
        }
        x0 = xx[0];
        x1 = xx[1];
        y0 = ff[0];
        y1 = ff[1];
    } else if (x > xx[n - 1]) { /* exceeds upper bound */
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            osh_error("osh_interpolate: upper bound exceeded.");
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            /* See comment above */
            return ff[n - 1];
        default:
            /* by default, do an extrapolate */
            break;
        }
        x0 = xx[n - 2];
        x1 = xx[n - 1];
        y0 = ff[n - 2];
        y1 = ff[n - 1];
    } else {
        /* binary search here, return index i to lower element*/
        i = osh_binary_search_d(x, xx, n);
        if (i < 0) {
            osh_error("osh_interpolate: bad binary_search. This should not happen.");
            return NAN;
        }
        /* the highest possible number i can achieve is len - 2. */
        x0 = xx[i];
        x1 = xx[i + 1];
        y0 = ff[i];
        y1 = ff[i + 1];
    }
    y = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    return y;
}

double osh_interpolate_dloglog(double xin, double const *xx, float const *ff, unsigned int n, int mode) {
    long int i;
    double x0, x1, y0, y1, log_x, log_x0, log_x1, t;

    if (xin < xx[0]) {
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            osh_error("osh_interpolate_dloglog: lower bound exceeded.");
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            return (double) ff[0];
        default: /* OSH_INTERPOLATE_OOB_EXTRAPOL: use first interval */
            x0 = xx[0];
            x1 = xx[1];
            y0 = (double) ff[0];
            y1 = (double) ff[1];
            break;
        }
    } else if (xin > xx[n - 1]) {
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            osh_error("osh_interpolate_dloglog: upper bound exceeded.");
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            return (double) ff[n - 1];
        default: /* OSH_INTERPOLATE_OOB_EXTRAPOL: use last interval */
            x0 = xx[n - 2];
            x1 = xx[n - 1];
            y0 = (double) ff[n - 2];
            y1 = (double) ff[n - 1];
            break;
        }
    } else {
        i = osh_binary_search_d(xin, xx, n);
        if (i < 0) {
            osh_error("osh_interpolate_dloglog: bad binary_search.");
            return NAN;
        }
        x0 = xx[i];
        x1 = xx[i + 1];
        y0 = (double) ff[i];
        y1 = (double) ff[i + 1];
    }

    /* Fall back to linear interpolation if endpoints are non-positive. */
    if (x0 <= 0.0 || x1 <= 0.0 || y0 <= 0.0 || y1 <= 0.0) {
        t = (xin - x0) / (x1 - x0);
        return y0 + t * (y1 - y0);
    }

    log_x = log(xin);
    log_x0 = log(x0);
    log_x1 = log(x1);
    t = (log_x - log_x0) / (log_x1 - log_x0);
    return exp(log(y0) + t * (log(y1) - log(y0)));
}

long int osh_binary_search_f(float x, float const *xx, unsigned long int len) {
    unsigned int n0, n1, n2; /* lower, middle, upper index */

    n0 = 0;
    n2 = len - 1;
    n1 = (n0 + n2) >> 1; /* gcc -O2: "i >> 1" produces less instructions than "i / 2" */

    while (n0 < n2) {

        if (x > xx[n1]) {
            n0 = n1;
        } else {
            n2 = n1;
        }

        n1 = (n0 + n2) >> 1;

        /* check if interval was found, if yes, exit */
        if ((n2 - n0) <= 1) {
            /* often last step is repeated even if element is found. Not sure if it can be sped up */
            return (long int) n0;
        }
    }
    return -1;
}

long int osh_binary_search_d(double x, double const *xx, unsigned long int len) {
    unsigned int n0, n1, n2; /* lower, middle, upper index */

    n0 = 0;
    n2 = len - 1;
    n1 = (n0 + n2) >> 1; /* gcc -O2: "i >> 1" produces less instructions than "i / 2" */

    while (n0 < n2) {

        if (x > xx[n1]) {
            n0 = n1;
        } else {
            n2 = n1;
        }

        n1 = (n0 + n2) >> 1;

        /* check if interval was found, if yes, exit */
        if ((n2 - n0) <= 1) {
            /* often last step is repeated even if element is found. Not sure if it can be sped up */
            return (long int) n0;
        }
    }
    return -1;
}

long int osh_binary_search_i2(int16_t x, int16_t const *xx, unsigned long int len) {
    unsigned long int n0, n1, n2;

    n0 = 0;
    n2 = len - 1;
    n1 = (n0 + n2) >> 1;

    while (n0 < n2) {

        if (x > xx[n1]) {
            n0 = n1;
        } else {
            n2 = n1;
        }

        n1 = (n0 + n2) >> 1;

        if ((n2 - n0) <= 1) {
            return (long int) n0;
        }
    }
    return -1;
}

long int osh_binary_search_upper_d(double x, double const *xx, unsigned long int len) {
    unsigned int n0, n1, n2; /* lower, middle, upper index */

    if (x <= xx[0])
        return 0;

    n0 = 0;
    n2 = len - 1;
    n1 = n2 >> 1; /* gcc -O2: "i >> 1" produces less instructions than "i / 2" */

    while (n0 < n2) {

        if (x > xx[n1]) {
            n0 = n1;
        } else {
            n2 = n1;
        }

        n1 = (n0 + n2) >> 1;

        /* check if interval was found, if yes, exit */
        if ((n2 - n0) <= 1) {
            /* often last step is repeated even if element is found. Not sure if it can be sped up */
            return (long int) n2;
        }
    }
    return -1;
}

long int osh_binary_search_nurep_d(double x, double const *xx, unsigned long int n) {

    unsigned long int ju, jm, jl;
    unsigned long j = 0;
    int ascnd;

    jl = 0;     /* Initialize lower ...*/
    ju = n + 1; /* ... and upper limits. */

    ascnd = (xx[n] >= xx[1]);

    while (ju - jl > 1) {    /* If we are not yet done,*/
        jm = (ju + jl) >> 1; /* compute a midpoint, */

        if ((x >= xx[jm]) == ascnd) /* paranthesis not in original code. This is tested and OK. */
            jl = jm;                /* and replace either the lower limit*/
        else
            ju = jm; /* or the upper limit, as appropriate. */
    } /* Repeat until the test condition is satisfied. */
    if (x == xx[1])
        j = 1; /* Then set the output */
    else if (x == xx[n])
        j = n - 1;
    else
        j = jl;

    return j;
}
