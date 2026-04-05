#include "common/osh_interpolate.h"

#include <math.h>

#include "common/osh_logger.h"

/**
 * @brief Linear interpolation in a given table. Single precision version.
 *
 * @details finds y for f(x) where x is given and f() is a table given as a set of *xx and *yy values.
 *
 * @param[in] xin - given x value
 * @param[in] *xx - sorted array holding x values, increasing and guaranteed no double occurences.
 * @param[in] *ff - corresponding array holding y values, i.e. *yy = f(*xx)
 * @param[in] n - size of *xx and *ff array in number of elements.
 * @param[in] mode - set behaviour if out-of-bounds, see OSH_INTERPOLATE_OOB_* defines. Default is extrapolation.
 *
 * @returns y = f(x)
 *
 * @warning *xx array must be increasing monotonically without double occurences.
 *
 * @author Niels Bassler
 */
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

/**
 * @brief Linear interpolation in a given table. Double precision version.
 *
 * @details finds y for f(x) where x is given and f() is a table given as a set of *xx and *yy values.
 *
 * @param[in] xin - given x value
 * @param[in] *xx - sorted array holding x values, increasing and guaranteed no double occurences.
 * @param[in] *ff - corresponding array holding y values, i.e. *yy = f(*xx)
 * @param[in] n - size of *xx and *ff array in number of elements.
 * @param[in] mode - set behaviour if out-of-bounds, see OSH_INTERPOLATE_OOB_* defines. Default is extrapolation.
 *
 * @returns y = f(x)
 *
 * @warning *xx array must be increasing monotonically without double occurences.
 *
 * @author Niels Bassler
 */
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

/**
 * @brief Binary search algorithm for moderately fast lookup. Single precision version.
 *
 * @param[in] x - value to look up in *xx
 * @param[in] *xx - sorted array
 * @param[in] len - size of *xx which is equal to number of available elements.
 *
 * @returns  lower index pointing to closest value. I.e. value is between xx[index] and xx[index + 1].
 *           If x <= xx[0] then index = 0 is returned.
 *           If x >= xx[len-1] then index = len-2 is returned.
 *           -1 if trouble.
 *
 * @author Niels Bassler
 */
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

/**
 * @brief Binary search algorithm for moderately fast lookup. Double precision version.
 *
 * @param[in] x - value to look up in *xx
 * @param[in] *xx - sorted array
 * @param[in] len - size of *xx which is equal to number of available elements.
 *
 * @returns  lower index pointing to closest value. I.e. value is between xx[index] and xx[index + 1].
 *           If x <= xx[0] then index = 0 is returned.
 *           If x >= xx[len-1] then index = len-2 is returned, e.g.:  [0|1|2|3|4]  -> 6 planes, last index is 4.
 *           -1 if trouble.
 *
 * @author Niels Bassler
 */
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

/**
 * @brief Binary search algorithm for moderately fast lookup. Double precision version.
 *
 * @param[in] x - value to look up in *xx
 * @param[in] *xx - sorted array
 * @param[in] len - size of *xx which is equal to number of available elements.
 *
 * @returns  upper index pointing to closest value. I.e. value is between xx[index] and xx[index + 1].
 *           If x <= xx[0] then index = 0 is returned.
 *           If x >= xx[len-1] then index = len-1 is returned
 *           -1 if trouble.
 *
 * @author Niels Bassler
 */
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

/**
 * @brief Binary search algorithm NumericalRecipies in C version. Double precision version.
 *
 * @details Note that array size (len) is number of elements (n) + 1.
 *
 * @param[in] x - value to look up in *xx
 * @param[in] *xx - sorted array, first element is stored in xx[1], last element in xx[n].
                    (Note than array size len = n + 1.
 * @param[in] n - number of available elements which is one less than the array size, i.e. xx[1...n].
 *
 * @returns  lower index pointing to closest value. I.e. value is between xx[index] and xx[index + 1].
 *           If x < xx[1] then index = 0 is returned.
 *           If x = xx[1] then index = 1 is returned.
 *           If x = xx[len-1] then index = n - 1 is returned.
 *           If x > xx[len-1] then index = n     is returned.
 *
 * @author NumericalRecipies, adapted by Niels Bassler
 */
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
