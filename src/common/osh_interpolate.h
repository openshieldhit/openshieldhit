#ifndef OSH_INTERPOLATE_H
#define OSH_INTERPOLATE_H

/**
 * @file osh_interpolate.h
 *
 * @brief Small interpolation and binary-search helpers for tabulated physics data.
 *
 * @details
 * The public API here is intentionally low-level: caller-owned sorted grids go
 * in, interpolated values or bracketing indices come out. These helpers are
 * used both during setup-time table generation and for a few runtime lookups.
 */

#include <stdint.h>

/* Define behaviour if out of bounds */
#define OSH_INTERPOLATE_OOB_ZERO 0     /* will return 0.0 if out of bounds */
#define OSH_INTERPOLATE_OOB_NEAREST 1  /* will return nearest endpoint value f(x) */
#define OSH_INTERPOLATE_OOB_EXTRAPOL 2 /* will return extrapolated value */
#define OSH_INTERPOLATE_OOB_ERR 255    /* will log an error and return NAN */

/**
 * @brief Linear interpolation in a single-precision table.
 *
 * @param[in] xin   Query x value.
 * @param[in] xx    Strictly increasing x-grid.
 * @param[in] ff    Corresponding y-values.
 * @param[in] len   Number of points in @p xx and @p ff.
 * @param[in] mode  Out-of-bounds behaviour, see OSH_INTERPOLATE_OOB_*.
 *
 * @returns Interpolated y = f(xin).
 *
 * @warning The @p xx array must be strictly increasing without duplicate values.
 */
double osh_interpolate_flin(float xin, float const *xx, float const *ff, unsigned int len, int mode);

/**
 * @brief Linear interpolation in a double-precision table.
 *
 * @param[in] xin   Query x value.
 * @param[in] xx    Strictly increasing x-grid.
 * @param[in] ff    Corresponding y-values.
 * @param[in] len   Number of points in @p xx and @p ff.
 * @param[in] mode  Out-of-bounds behaviour, see OSH_INTERPOLATE_OOB_*.
 *
 * @returns Interpolated y = f(xin).
 *
 * @warning The @p xx array must be strictly increasing without duplicate values.
 */
double osh_interpolate_dlin(double xin, double const *xx, double const *ff, unsigned int len, int mode);

/**
 * @brief Log-log interpolation with a double-precision x-grid and float y-values.
 *
 * @details
 * Interpolates in log(x)–log(y) space, which is appropriate for quantities
 * that are well approximated by a power law between tabulated points, such as
 * mass stopping powers from LOADDEDX source tables.
 *
 * Both the x and y values at the bracket endpoints must be strictly positive
 * for the log-log path to be taken.  If either is non-positive the function
 * falls back to linear interpolation between those two points.
 *
 * Out-of-bounds behaviour is controlled by @p mode exactly as in the linear
 * variants.  The typical usage for resampling stopping-power tables is
 * OSH_INTERPOLATE_OOB_NEAREST (flat extrapolation).
 *
 * @param[in] xin  Query x value.
 * @param[in] xx   Strictly increasing x-grid (double).
 * @param[in] ff   Corresponding y-values (float).
 * @param[in] len  Number of points.
 * @param[in] mode OOB behaviour: OSH_INTERPOLATE_OOB_*.
 *
 * @returns Interpolated y value.
 */
double osh_interpolate_dloglog(double xin, double const *xx, float const *ff, unsigned int len, int mode);

/**
 * @brief Binary search on a strictly increasing float array.
 *
 * @param[in] x    Query value.
 * @param[in] xx   Sorted array.
 * @param[in] len  Number of points.
 *
 * @returns Lower bracketing index such that @p x lies between
 *          `xx[index]` and `xx[index + 1]`. If `x <= xx[0]`, returns `0`.
 *          If `x >= xx[len - 1]`, returns `len - 2`. Returns `-1` on error.
 */
long int osh_binary_search_f(float x, float const *xx, unsigned long int len);

/**
 * @brief Binary search on a strictly increasing double array.
 *
 * @param[in] x    Query value.
 * @param[in] xx   Sorted array.
 * @param[in] len  Number of points.
 *
 * @returns Lower bracketing index such that @p x lies between
 *          `xx[index]` and `xx[index + 1]`. If `x <= xx[0]`, returns `0`.
 *          If `x >= xx[len - 1]`, returns `len - 2`. Returns `-1` on error.
 */
long int osh_binary_search_d(double x, double const *xx, unsigned long int len);

/**
 * @brief Binary search returning the upper bracketing index.
 *
 * @param[in] x    Query value.
 * @param[in] xx   Sorted array.
 * @param[in] len  Number of points.
 *
 * @returns Upper bracketing index. If `x <= xx[0]`, returns `0`.
 *          If `x >= xx[len - 1]`, returns `len - 1`. Returns `-1` on error.
 */
long int osh_binary_search_upper_d(double x, double const *xx, unsigned long int len);

/**
 * @brief Binary search on a sorted int16 array, returning the lower bracket.
 *
 * @param[in] x    Query value.
 * @param[in] xx   Sorted array.
 * @param[in] len  Number of points.
 *
 * @returns Lower bracketing index. If `x <= xx[0]`, returns `0`.
 *          If `x >= xx[len - 1]`, returns `len - 2`. Returns `-1` on error.
 */
long int osh_binary_search_i2(int16_t x, int16_t const *xx, unsigned long int len);

/**
 * @brief Numerical Recipes style binary search on a 1-based double array.
 *
 * @details
 * The array layout follows the classic Numerical Recipes convention:
 * valid data are stored in `xx[1] ... xx[n]`, so the backing array has size
 * `n + 1`.
 *
 * @param[in] x   Query value.
 * @param[in] xx  Sorted 1-based array.
 * @param[in] n   Number of valid elements.
 *
 * @returns Lower bracketing index in the Numerical Recipes convention.
 */
long int osh_binary_search_nurep_d(double x, double const *xx, unsigned long int n);

#endif /* OSH_INTERPOLATE_H */
