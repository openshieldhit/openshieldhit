#include <math.h>
#include <stdio.h>

#include "common/osh_vect.h"
#include "test_assert.h"

void test_dot_product(void) {
    double a[OSH_VECT_DIM] = {1.0, 2.0, 3.0};
    double b[OSH_VECT_DIM] = {4.0, -5.0, 6.0};
    ASSERT_TRUE(fabs(osh_vect_dot(a, b) - (1 * 4 + 2 * (-5) + 3 * 6)) < OSH_VECT_EPS);
}

void test_cross_product(void) {
    double i[OSH_VECT_DIM] = {1.0, 0.0, 0.0};
    double j[OSH_VECT_DIM] = {0.0, 1.0, 0.0};
    double k[OSH_VECT_DIM];
    osh_vect_cross(i, j, k);
    ASSERT_TRUE(fabs(k[0] - 0.0) < OSH_VECT_EPS);
    ASSERT_TRUE(fabs(k[1] - 0.0) < OSH_VECT_EPS);
    ASSERT_TRUE(fabs(k[2] - 1.0) < OSH_VECT_EPS);
}

void test_norm(void) {
    double u[OSH_VECT_DIM] = {3.0, 4.0, 0.0};
    double v[OSH_VECT_DIM];
    ASSERT_TRUE(fabs(osh_vect_len2(u) - 25.0) < OSH_VECT_EPS); /* 3**2 + 4**2 = 25 */

    osh_vect_norm2(u, v); /* normalize u into v */

    ASSERT_TRUE(fabs(osh_vect_len2(v) - 1.0) < OSH_VECT_EPS); /* should be unit vector */
}

void test_affine_bzalign_transform(void) {
    double p_local[3] = {1.0, 2.0, 3.0};
    double v_local[3] = {0.0, 0.0, 1.0};
    double origin_local[3] = {4.0, 5.0, 6.0};
    double zdir_world[3] = {0.0, 0.0, 1.0};
    double tm[16];
    double p_world[3];
    double v_world[3];

    osh_vect_setup_tmatrix_bzalign_affine(origin_local, zdir_world, tm);
    osh_vect_trans_point_affine(p_local, p_world, tm);
    osh_vect_trans_vector_affine(v_local, v_world, tm);

    ASSERT_TRUE(fabs(p_world[0] - 5.0) < OSH_VECT_EPS);
    ASSERT_TRUE(fabs(p_world[1] - 7.0) < OSH_VECT_EPS);
    ASSERT_TRUE(fabs(p_world[2] - 9.0) < OSH_VECT_EPS);

    ASSERT_TRUE(fabs(v_world[0] - 0.0) < OSH_VECT_EPS);
    ASSERT_TRUE(fabs(v_world[1] - 0.0) < OSH_VECT_EPS);
    ASSERT_TRUE(fabs(v_world[2] - 1.0) < OSH_VECT_EPS);
}

int main(void) {
    printf("Running osh_vect tests...\n");

    test_dot_product();
    test_cross_product();
    test_norm();
    test_affine_bzalign_transform();

    printf("All tests passed.\n");
    return 0;
}
