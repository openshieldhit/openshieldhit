#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "osh_vect.h"

void test_dot_product() {
    double a[OSH_VECT_DIM] = {1.0, 2.0, 3.0};
    double b[OSH_VECT_DIM] = {4.0, -5.0, 6.0};
    double result = osh_vect_dot(a, b);
    assert(fabs(result - (1*4 + 2*(-5) + 3*6)) < OSH_VECT_EPS);
}


void test_cross_product() {
    double i[OSH_VECT_DIM] = {1.0, 0.0, 0.0};
    double j[OSH_VECT_DIM] = {0.0, 1.0, 0.0};
    double k[OSH_VECT_DIM];
    osh_vect_cross(i, j, k);
    assert(fabs(k[0] - 0.0) < OSH_VECT_EPS);
    assert(fabs(k[1] - 0.0) < OSH_VECT_EPS);
    assert(fabs(k[2] - 1.0) < OSH_VECT_EPS);
}


void test_norm() {
    double len2_after;
    double u[OSH_VECT_DIM] = {3.0, 4.0, 0.0};
    double v[OSH_VECT_DIM];

    double len2_before = osh_vect_len2(u);
    assert(fabs(len2_before - 25.0) < OSH_VECT_EPS);  /* 3**2 + 4**2 = 25 */

    osh_vect_norm2(u, v);  /* normalize u into v */

    len2_after = osh_vect_len2(v);
    assert(fabs(len2_after - 1.0) < OSH_VECT_EPS);  /* should be unit vector */
}


int main(void) {
    printf("Running osh_vect tests...\n");

    test_dot_product();
    test_cross_product();
    test_norm();

    printf("All tests passed.\n");
    return 0;
}

