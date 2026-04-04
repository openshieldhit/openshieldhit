#ifndef OSH_INTERPOLATE_H
#define OSH_INTERPOLATE_H

#include <stdint.h>

/* Define behaviour if out of bounds */
#define OSH_INTERPOLATE_OOB_ZERO 0     /* will return 0.0 if out of bounds */
#define OSH_INTERPOLATE_OOB_NEAREST 1  /* will return lowest or highest bound */
#define OSH_INTERPOLATE_OOB_EXTRAPOL 2 /* will return extrapolated value */
#define OSH_INTERPOLATE_OOB_ERR 255    /* will throw an error and exit code */

double osh_interpolate_flin(float xin, float const *xx, float const *ff, unsigned int len, int mode);
double osh_interpolate_dlin(double xin, double const *xx, double const *ff, unsigned int len, int mode);

long int osh_binary_search_f(float x, float const *xx, unsigned long int len);
long int osh_binary_search_d(double x, double const *xx, unsigned long int len);
long int osh_binary_search_upper_d(double x, double const *xx, unsigned long int len);
long int osh_binary_search_i2(int16_t x, int16_t const *xx, unsigned long int len);

long int osh_binary_search_nurep_d(double x, double const *xx, unsigned long int n);

#endif /* OSH_INTERPOLATE_H */
