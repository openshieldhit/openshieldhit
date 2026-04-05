#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/osh_interpolate.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_interpolate_dlin(void) {
    double xx[3] = {0.0, 1.0, 2.0};
    double yy[3] = {0.0, 10.0, 20.0};
    double y;

    y = osh_interpolate_dlin(0.5, xx, yy, 3u, OSH_INTERPOLATE_OOB_EXTRAPOL);
    ASSERT_TRUE(fabs(y - 5.0) < 1e-12);

    y = osh_interpolate_dlin(-1.0, xx, yy, 3u, OSH_INTERPOLATE_OOB_ZERO);
    ASSERT_TRUE(fabs(y) < 1e-12);

    y = osh_interpolate_dlin(3.0, xx, yy, 3u, OSH_INTERPOLATE_OOB_NEAREST);
    ASSERT_TRUE(fabs(y - 20.0) < 1e-12);

    y = osh_interpolate_dlin(-3.0, xx, yy, 3u, OSH_INTERPOLATE_OOB_NEAREST);
    ASSERT_TRUE(fabs(y - 0.0) < 1e-12);
}

static void test_interpolate_flin(void) {
    float xx[3] = {0.0f, 1.0f, 2.0f};
    float yy[3] = {0.0f, 10.0f, 20.0f};
    double y;

    y = osh_interpolate_flin(3.0f, xx, yy, 3u, OSH_INTERPOLATE_OOB_NEAREST);
    ASSERT_TRUE(fabs(y - 20.0) < 1e-6);

    y = osh_interpolate_flin(-3.0f, xx, yy, 3u, OSH_INTERPOLATE_OOB_NEAREST);
    ASSERT_TRUE(fabs(y - 0.0) < 1e-6);
}

static void test_binary_search_upper_d(void) {
    double xx[4] = {1.0, 3.0, 10.0, 20.0};

    ASSERT_TRUE(osh_binary_search_upper_d(0.0, xx, 4u) == 0);
    ASSERT_TRUE(osh_binary_search_upper_d(1.0, xx, 4u) == 0);
    ASSERT_TRUE(osh_binary_search_upper_d(1.5, xx, 4u) == 1);
    ASSERT_TRUE(osh_binary_search_upper_d(3.0, xx, 4u) == 1);
    ASSERT_TRUE(osh_binary_search_upper_d(9.9, xx, 4u) == 2);
    ASSERT_TRUE(osh_binary_search_upper_d(20.0, xx, 4u) == 3);
    ASSERT_TRUE(osh_binary_search_upper_d(25.0, xx, 4u) == 3);
}

static void test_binary_search_i2(void) {
    int16_t xx[5] = {-1000, -100, 0, 500, 1600};

    ASSERT_TRUE(osh_binary_search_i2(-1500, xx, 5u) == 0);
    ASSERT_TRUE(osh_binary_search_i2(-1000, xx, 5u) == 0);
    ASSERT_TRUE(osh_binary_search_i2(-50, xx, 5u) == 1);
    ASSERT_TRUE(osh_binary_search_i2(0, xx, 5u) == 1);
    ASSERT_TRUE(osh_binary_search_i2(1200, xx, 5u) == 3);
    ASSERT_TRUE(osh_binary_search_i2(1600, xx, 5u) == 3);
}

int main(void) {
    test_interpolate_dlin();
    test_interpolate_flin();
    test_binary_search_upper_d();
    test_binary_search_i2();
    return 0;
}
